"""Project commit boundaries under the existing flowcell execution lease."""

import hashlib
import os
import shutil
import uuid

from contextlib import contextmanager, nullcontext
from pathlib import Path

from bcl2fastq_pipeline import (
    analysis_snapshots,
    notification_delivery,
    notifications,
    processing_times,
    project_corrections,
    project_scheduling,
)
from bcl2fastq_pipeline.state import StateConflictError, project_analysis_completed_at, utcnow


def _copy(source, target):
    temporary = target.with_name(f".{target.name}.tmp")
    shutil.copy2(source, temporary)
    with temporary.open("rb") as handle:
        os.fsync(handle.fileno())
    os.replace(temporary, target)
    descriptor = os.open(target.parent, os.O_RDONLY)
    try:
        os.fsync(descriptor)
    finally:
        os.close(descriptor)


def completed_deliveries(state):
    return {
        project
        for project, progress in state.get("project_progress", {}).items()
        if progress.get("finalization", {}).get("status") == "completed"
    }


def recovery_projects(state, stage):
    if stage not in {"reporting", "finalization"} or state["status"] not in {
        "failed",
        "interrupted",
    }:
        return set()
    return completed_deliveries(state)


class ProjectLifecycle:
    def __init__(self, cfg, store):
        self.cfg = cfg
        self.store = store
        self.run_id = cfg.run.run_id
        self.delivery_mode = False
        self.analysis_timer = None

    def plan_analysis(self, projects):
        def update(state):
            state["projects"] = sorted(set(state["projects"]) | set(projects))
            return state

        self.store.mutate(self.run_id, update)

    @property
    def analysis_committed(self):
        selected = self.cfg.run.selected_projects
        return any(
            progress.get("analysis", {}).get("status") == "completed"
            for project, progress in self.store.read(self.run_id)
            .get("project_progress", {})
            .items()
            if selected is None or project in selected
        )

    def start_analysis(self, project, identity, workflow):
        schedule = self.cfg.run.project_schedule or {}
        self.store.start_project_analysis(
            self.run_id,
            project,
            uuid.uuid4().hex,
            identity,
            workflow,
            self.cfg.output_path,
            schedule.get("projects", {}).get(project),
        )

    def _payload(self, project, kind):
        from bcl2fastq_pipeline import workflow_config  # noqa: PLC0415

        cfg = self.cfg
        payload = notifications.make_payload(
            cfg,
            kind,
            projects=[project],
            processing_timing=processing_times.snapshot(self.store.read(self.run_id)),
        )
        context = workflow_config.delivery_project_config(cfg, project)
        custom = dict(context.custom if context else cfg.run.custom)
        payload.update(
            custom=custom,
            libprep=context.libprep if context else cfg.run.libprep,
            user=custom.get("User", cfg.run.user),
        )
        if context:
            payload.update(
                sample_sheet=str(context.sample_sheet),
                sample_submission_form=str(context.submission_form),
            )
        payload["project_contexts"] = {
            project: {
                "libprep": payload["libprep"],
                "workflow": workflow_config.project_workflow(cfg, project),
                "custom": custom,
                "sample_sheet": payload["sample_sheet"],
                "submission_form": payload["sample_submission_form"],
            }
        }
        return payload

    def complete_analysis(self, project):
        from bcl2fastq_pipeline import analysis_qc  # noqa: PLC0415

        cfg = self.cfg
        progress = self.store.read(self.run_id)["project_progress"][project]["analysis"]
        generation = progress["generation"]
        destination = (
            Path(self.store.read(self.run_id)["output_path"])
            / "provenance"
            / "project-analysis"
            / generation
            / project
        )
        for directory in (
            destination.parent.parent.parent,
            destination.parent.parent,
            destination.parent,
        ):
            if directory.is_symlink():
                raise StateConflictError(
                    "Project analysis provenance directories must not be symlinks"
                )
        destination.mkdir(parents=True, exist_ok=False)
        payload = self._payload(project, "processed")
        date = self.run_id.split("_", 1)[0]
        sources = [cfg.output_path / f"multiqc_{project}_{date}.html"]
        sources.extend(
            path
            for path in (
                cfg.output_path / f"all_samples_web_summary_{project}_{date}.html",
                cfg.output_path / f"configmaker-analysis-{project}.json",
                cfg.output_path / f"{project}_samplesheet.tsv",
            )
            if path.is_file()
        )
        for key in ("sample_sheet", "sample_submission_form"):
            if payload.get(key):
                source = Path(payload[key])
                sources.append(source)
                payload[key] = str(destination / source.name)
        for source in sources:
            _copy(source, destination / source.name)
        for directory in (
            destination.parent,
            destination.parent.parent,
            destination.parent.parent.parent,
        ):
            descriptor = os.open(directory, os.O_RDONLY)
            try:
                os.fsync(descriptor)
            finally:
                os.close(descriptor)
        payload.update(
            output_path=str(destination),
            generation=generation,
            attempt=progress["attempt"],
            analysis_qc=analysis_qc.collect(cfg, [project]),
        )
        payload["artifact_hashes"] = {
            str(path): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in destination.iterdir()
        }
        context = payload["project_contexts"][project]
        context.update(
            sample_sheet=payload["sample_sheet"], submission_form=payload["sample_submission_form"]
        )
        self.store.complete_project_analysis(self.run_id, project, generation, payload)
        with self.analysis_timer.exclude() if self.analysis_timer is not None else nullcontext():
            notification_delivery.deliver_pending(cfg, self.store, self.run_id, kind="processed")

    def _hydrate_analysis(self, project):
        from bcl2fastq_pipeline import workflow_config  # noqa: PLC0415

        state = self.store.read(self.run_id)
        if (
            state.get("project_progress", {}).get(project, {}).get("analysis", {}).get("status")
            == "completed"
        ):
            return
        completed = project_analysis_completed_at(state, project)
        if not completed:
            raise StateConflictError(f"No completed analysis for project {project}")
        metadata = state["stages"]["analysis"]["metadata"]
        analysis = {
            "status": "completed",
            "generation": hashlib.sha256(
                f"{self.run_id}:{project}:{completed}".encode()
            ).hexdigest(),
            "attempt": state["attempt"],
            "completed_at": completed,
            "legacy": True,
            "workdir": metadata.get("workdirs", {}).get(project),
            "workflow": workflow_config.project_workflow(self.cfg, project),
        }

        def update(current):
            current.setdefault("project_progress", {}).setdefault(project, {})["analysis"] = (
                analysis
            )
            return current

        self.store.mutate(self.run_id, update)

    @contextmanager
    def _measure_project(self, project, step):
        with processing_times.measure(self.store, self.run_id, step):
            yield

        def accumulated(state):
            record = state["processing_timings"][step][-1]
            state["project_progress"][project].setdefault("timings", {})[step] = record[
                "duration_seconds"
            ]
            record["duration_seconds"] = sum(
                progress.get("timings", {}).get(step, 0)
                for progress in state["project_progress"].values()
            )
            return state

        self.store.mutate(self.run_id, accumulated)

    def finalize(self):
        from bcl2fastq_pipeline import afterFastq  # noqa: PLC0415

        self.delivery_mode = True
        cfg = self.cfg
        projects = afterFastq.get_project_names(afterFastq.get_project_dirs(cfg))
        for project in project_scheduling.project_order(cfg, projects):
            if project in completed_deliveries(self.store.read(self.run_id)):
                self._inventory(project)
                continue
            self._hydrate_analysis(project)

            def starting(state):
                state["project_progress"][project]["finalization"] = {
                    "status": "running",
                    "started_at": utcnow(),
                }
                return state

            self.store.mutate(self.run_id, starting)
            with self._measure_project(project, "archiving"):
                afterFastq.archive_worker(cfg, projects=[project])
            date = self.run_id.split("_", 1)[0]
            with self._measure_project(project, "archive_checksums"):
                for name in (f"{project}_{date}.7za", f"QC_{project}_{date}.7za"):
                    archive = cfg.output_path / name
                    if not archive.is_file() or archive.stat().st_size == 0:
                        raise RuntimeError(f"Project archive is missing or empty: {archive}")
                    afterFastq.md5sum_archive(archive)
            state = self.store.read(self.run_id)
            snapshots = analysis_snapshots.prepare(state, [project])
            analysis = state["project_progress"][project]["analysis"]
            if snapshots[project]["status"] != "available" and not analysis.get("legacy"):
                raise RuntimeError(f"Successful analysis snapshot unavailable for {project}")
            project_corrections.publish_project(self.store, self.run_id, project)
            payload = self._payload(project, "finalized")
            payload.update(
                output_path=state["output_path"],
                generation=analysis["generation"],
                attempt=state["attempt"],
            )
            self.store.complete_project_finalization(
                self.run_id, project, payload, snapshots[project]
            )
            project_corrections.recover(self.store, self.run_id)
            analysis_snapshots.recover(self.store, self.run_id)
            self._inventory(project)
            notification_delivery.deliver_pending(cfg, self.store, self.run_id, kind="finalized")

    def _inventory(self, project):
        from flowcell_manager.flowcell_manager import (  # noqa: PLC0415
            add_flowcell,
            list_flowcell_all,
        )

        state = self.store.read(self.run_id)
        timestamp = state["project_completion"][project]["completed_at"]
        inventory = list_flowcell_all(state["output_path"])
        if (
            (inventory["project"] == project)
            & (inventory["timestamp"] == timestamp)
            & (inventory["archived"] == "0")
        ).any():
            return
        add_flowcell(project=project, path=state["output_path"], timestamp=timestamp)

    def complete_flowcell(self):
        from bcl2fastq_pipeline import afterFastq  # noqa: PLC0415

        projects = afterFastq.get_project_names(afterFastq.get_project_dirs(self.cfg))
        if set(projects) - completed_deliveries(self.store.read(self.run_id)):
            raise StateConflictError("Project deliveries remain incomplete")
        return self.store.complete_run(self.run_id, sorted(projects), project_boundaries=True)
