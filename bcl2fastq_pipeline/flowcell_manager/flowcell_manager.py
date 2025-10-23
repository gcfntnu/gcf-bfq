#!/usr/bin/env python

import argparse
import datetime
import shutil
import subprocess

from pathlib import Path

import pandas as pd

from bcl2fastq_pipeline.config import PipelineConfig

pd.set_option("display.max_rows", 5000)
pd.set_option("display.max_columns", 6)


def get_cfg():
    """Return the active PipelineConfig, loading a static one if needed."""
    try:
        return PipelineConfig.get()
    except RuntimeError:
        # Not initialized yet → load static config only
        return PipelineConfig.load("/config/bcl2fastq.ini")


def add_flowcell(**args):
    cfg = get_cfg()
    row_list = [
        {
            "project": args["project"],
            "flowcell_path": args["path"],
            "timestamp": args["timestamp"],
            "archived": 0,
        }
    ]
    df = pd.DataFrame(row_list)
    flowcells_processed = pd.read_csv(cfg.static.paths.manager_dir / "flowcells.processed")
    flowcells_processed = pd.concat([flowcells_processed, df], sort=True)
    flowcells_processed.to_csv(
        cfg.static.paths.manager_dir / "flowcells.processed",
        index=False,
        columns=["project", "flowcell_path", "timestamp", "archived"],
    )


def archive_flowcell(**args):
    force = args.get("force", False)
    flowcell = Path(args["flowcell"])
    cfg = get_cfg()
    flowcells_processed = pd.read_csv(cfg.static.paths.manager_dir / "flowcells.processed")
    fc_for_deletion = flowcells_processed.loc[flowcells_processed["flowcell_path"] == str(flowcell)]
    if fc_for_deletion.empty:
        print("No such flowcell in inventory!")
        return
    confirm = "yes"
    if not force:
        print("Please confirm deletion of the following flowcell and the contained projects.\n")
        print(fc_for_deletion)
        confirm = input("Delete? (yes/no): ").lower()
    if confirm == "yes":
        deletions = [(flowcell / pid) for pid in fc_for_deletion["project"]]
        deletions += list(flowcell.rglob("*.bam"))
        deletions += list(flowcell.rglob("*.bam.bai"))
        deletions += list(flowcell.glob("*.fastq.gz"))
        deletions += list(flowcell.glob("*.7za"))
        print(f"DELETING: {len(deletions)} items from {flowcell}")
        for i, item in enumerate(deletions, 1):
            print(f"[{i}/{len(deletions)}] Deleting {item}")

        for d in deletions:
            if d.is_dir():
                shutil.rmtree(d)
            else:
                d.unlink()

        flowcells_processed.loc[
            flowcells_processed["flowcell_path"] == str(flowcell), "archived"
        ] = datetime.datetime.now()
        flowcells_processed.to_csv(
            cfg.static.paths.manager_dir / "flowcells.processed",
            index=False,
            columns=["project", "flowcell_path", "timestamp", "archived"],
        )
    else:
        print("Skipping...")


def rerun_flowcell(**args):
    force = args.get("force", False)
    flowcell = args["flowcell"]
    cfg = get_cfg()
    flowcells_processed = pd.read_csv(cfg.static.paths.manager_dir / "flowcells.processed")
    fc_for_deletion = flowcells_processed.loc[flowcells_processed["flowcell_path"] == flowcell]
    if fc_for_deletion.empty:
        print("No such flowcell in inventory!")
        return
    confirm = "yes"
    if not force:
        print("Please confirm deletion of the following flowcell and the contained projects.\n")
        print(fc_for_deletion)
        confirm = input("Delete? (yes/no): ").lower()
    if confirm == "yes":
        cmd = f"rm -rf {flowcell}"
        print(f"DELETING FLOWCELL: {cmd}")
        subprocess.check_call(cmd, shell=True)
        flowcells_processed = flowcells_processed.loc[
            flowcells_processed["flowcell_path"] != flowcell
        ]
        flowcells_processed.to_csv(
            cfg.static.paths.manager_dir / "flowcells.processed",
            index=False,
            columns=["project", "flowcell_path", "timestamp", "archived"],
        )
    else:
        print("Skipping...")


def list_processed(**args):
    cfg = get_cfg()
    flowcells_processed = pd.read_csv(cfg.static.paths.manager_dir / "flowcells.processed")
    return flowcells_processed.loc[
        (flowcells_processed["timestamp"] != "0") | (flowcells_processed["archived"] != "0")
    ]


def list_all(**args):
    cfg = get_cfg()
    return pd.read_csv(cfg.static.paths.manager_dir / "flowcells.processed")


def list_project(project):
    cfg = get_cfg()
    flowcells_processed = pd.read_csv(cfg.static.paths.manager_dir / "flowcells.processed")
    return flowcells_processed.loc[
        (flowcells_processed["project"] == project) & (flowcells_processed["timestamp"] != "0")
    ]


def list_flowcell(flowcell):
    cfg = get_cfg()
    flowcells_processed = pd.read_csv(cfg.static.paths.manager_dir / "flowcells.processed")
    return flowcells_processed.loc[
        (flowcells_processed["flowcell_path"] == flowcell)
        & (flowcells_processed["timestamp"] != "0")
    ]


def list_flowcell_all(flowcell):
    cfg = get_cfg()
    flowcells_processed = pd.read_csv(cfg.static.paths.manager_dir / "flowcells.processed")
    # USED TO AVOID RUNNING OLD FLOWCELLS
    return flowcells_processed.loc[flowcells_processed["flowcell_path"] == flowcell]


def pretty_print(df):
    print("Project \t Flowcell path \t Timestamp \t Archived")
    for i, row in df.iterrows():
        print(
            "{}\t{}\t{}\t{}".format(
                row["project"], row["flowcell_path"], row["timestamp"], row["archived"]
            )
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers()
    parser_add = subparsers.add_parser("add", help="Add a project to the inventory file.")
    parser_add.set_defaults(func=add_flowcell)
    parser_add.add_argument("project", type=str, help="GCF project number.")
    parser_add.add_argument("path", type=str, help="Flowcell path.")
    parser_add.add_argument(
        "timestamp", type=datetime.datetime.fromisoformat, help="A datetime isoformat string."
    )

    parser_archive = subparsers.add_parser(
        "archive", help="Archive the flowcell by deleting fastq files and .7za archives."
    )
    parser_archive.set_defaults(func=archive_flowcell)
    parser_archive.add_argument("flowcell", type=str, help="Path to flowcell to be archived.")
    parser_archive.add_argument("--force", action="store_true", help="Force archive (no prompt)")

    parser_rerun = subparsers.add_parser(
        "rerun", help="Rerun the flowcell by deleting the output directory."
    )
    parser_rerun.set_defaults(func=rerun_flowcell)
    parser_rerun.add_argument("flowcell", type=str, help="Path to flowcell to be deleted.")
    parser_rerun.add_argument("--force", action="store_true", help="Force archive (no prompt)")

    parser_list = subparsers.add_parser("list", help="List all flowcells.")
    parser_list.set_defaults(func=list_all, print_res=True)

    parser_list_processed = subparsers.add_parser(
        "list-processed",
        help="List only flowcells in the inventory file processed by bfq pipeline.",
    )
    parser_list_processed.set_defaults(func=list_processed, print_res=True)

    args = parser.parse_args()
    if vars(args).get("print_res", False):
        pretty_print(args.func(**vars(args)))
    else:
        args.func(**vars(args))
