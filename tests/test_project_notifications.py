"""A selected project's notification uses its own metadata."""

from bcl2fastq_pipeline.config import ProjectConfig

from bcl2fastq_pipeline import notifications
from tests.support.fixtures import configured_bfq


def test_selected_project_notification_does_not_inherit_primary_project_metadata(tmp_path):
    cfg, _source, output = configured_bfq(tmp_path)
    first, second = "GCF-2099-001", "GCF-2099-002"
    cfg.run.user = "primary-owner"
    for project, user, kit, workflow in (
        (first, "primary-owner", "Kit A", "default"),
        (second, "selected-owner", "Kit B", "wgs"),
    ):
        context = ProjectConfig(
            project,
            output / f"SampleSheet_{project}.csv",
            output / f"{project}_Sample-Submission-Form.xlsx",
            {"User": user, "Libprep": kit},
            (151, 151),
        )
        context.pipeline = workflow
        cfg.run.project_configs[project] = context
    payload = notifications.make_payload(cfg, "processed", projects=[second])
    assert set(payload["project_contexts"]) == {second}
    rendered = notifications._run_summary(payload)
    assert "selected-owner" in rendered and "Kit B" in rendered and "wgs" in rendered
    assert "primary-owner" not in rendered and "Kit A" not in rendered
