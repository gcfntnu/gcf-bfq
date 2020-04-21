#!/usr/bin/env python

import bcl2fastq_pipeline.getConfig
import sys
import pandas as pd
import os
import subprocess
import datetime
import argparse
import glob


pd.set_option('display.max_rows', 5000)
pd.set_option('display.max_columns', 6)

def add_flowcell(**args):
    config = bcl2fastq_pipeline.getConfig.getConfig()
    row_list = [
            {
            'project': args['project'],
            'flowcell_path': args['path'],
            'timestamp': args['timestamp'],
            'archived': 0
            }
            ]
    df = pd.DataFrame(row_list)
    flowcells_processed = pd.read_csv(os.path.join(config.get("FlowCellManager","managerDir"),'flowcells.processed'))
    flowcells_processed = flowcells_processed.append(df,sort=True)
    flowcells_processed.to_csv(
            os.path.join(config.get("FlowCellManager","managerDir"),'flowcells.processed'),
            index=False,
            columns = ['project','flowcell_path','timestamp', 'archived'],
            )


def archive_flowcell(**args):
    force = args.get("force",False)
    flowcell = args['flowcell']
    config = bcl2fastq_pipeline.getConfig.getConfig()
    flowcells_processed = pd.read_csv(os.path.join(config.get("FlowCellManager","managerDir"),'flowcells.processed'))
    fc_for_deletion = flowcells_processed.loc[flowcells_processed['flowcell_path'] == flowcell]
    if fc_for_deletion.empty:
        print("No such flowcell in inventory!")
        return
    confirm = 'yes'
    if not force:
        print("Please confirm deletion of the following flowcell and the contained projects.\n")
        print(fc_for_deletion)
        confirm = input("Delete? (yes/no): ").lower()
    if confirm == 'yes':
        deletions = [os.path.join(flowcell,pid) for pid in fc_for_deletion['project']]
        bam = glob.glob(os.path.join(flowcell, "**", '*.bam'))
        deletions.extend(bam)
        bam_bai = glob.glob(os.path.join(flowcell, "**", '*.bam.bai'))
        deletions.extend(bam_bai)
        deletions.append("{}/*.fastq.gz".format(flowcell))
        deletions.append("{}/*.7za".format(flowcell))
        cmd = "rm -rf {}".format(" ".join(deletions))
        print("DELETING: {}".format(cmd))
        subprocess.check_call(cmd, shell=True)
        flowcells_processed.loc[flowcells_processed['flowcell_path'] == flowcell,'archived'] = datetime.datetime.now()
        flowcells_processed.to_csv(
            os.path.join(config.get("FlowCellManager","managerDir"),'flowcells.processed'),
            index=False,
            columns = ['project','flowcell_path','timestamp','archived'],
            )
    else:
        print("Skipping...")

def rerun_flowcell(**args):
    force = args.get("force",False)
    flowcell = args['flowcell']
    config = bcl2fastq_pipeline.getConfig.getConfig()
    flowcells_processed = pd.read_csv(os.path.join(config.get("FlowCellManager","managerDir"),'flowcells.processed'))
    fc_for_deletion = flowcells_processed.loc[flowcells_processed['flowcell_path'] == flowcell]
    if fc_for_deletion.empty:
        print("No such flowcell in inventory!")
        return
    confirm = 'yes'
    if not force:
        print("Please confirm deletion of the following flowcell and the contained projects.\n")
        print(fc_for_deletion)
        confirm = input("Delete? (yes/no): ").lower()
    if confirm == 'yes':
        cmd = "rm -rf {}".format(flowcell)
        print("DELETING FLOWCELL: {}".format(cmd))
        subprocess.check_call(cmd, shell=True)
        flowcells_processed = flowcells_processed.loc[flowcells_processed['flowcell_path'] != flowcell]
        flowcells_processed.to_csv(
            os.path.join(config.get("FlowCellManager","managerDir"),'flowcells.processed'),
            index=False,
            columns = ['project','flowcell_path','timestamp','archived'],
            )
    else:
        print("Skipping...")

def list_processed(**args):
    config = bcl2fastq_pipeline.getConfig.getConfig()
    flowcells_processed = pd.read_csv(os.path.join(config.get("FlowCellManager","managerDir"),'flowcells.processed'))
    return flowcells_processed.loc[(flowcells_processed['timestamp'] != '0') | (flowcells_processed['archived'] != '0')]

def list_all(**args):
    config = bcl2fastq_pipeline.getConfig.getConfig()
    return pd.read_csv(os.path.join(config.get("FlowCellManager","managerDir"),'flowcells.processed'))

def list_project(project):
    config = bcl2fastq_pipeline.getConfig.getConfig()
    flowcells_processed = pd.read_csv(os.path.join(config.get("FlowCellManager","managerDir"),'flowcells.processed'))
    return flowcells_processed.loc[(flowcells_processed['project'] == project) & (flowcells_processed['timestamp'] != '0')]

def list_flowcell(flowcell):
    config = bcl2fastq_pipeline.getConfig.getConfig()
    flowcells_processed = pd.read_csv(os.path.join(config.get("FlowCellManager","managerDir"),'flowcells.processed'))
    return flowcells_processed.loc[(flowcells_processed['flowcell_path'] == flowcell) & (flowcells_processed['timestamp'] != '0')]

def list_flowcell_all(flowcell):
    #USED TO AVOID RUNNING OLD FLOWCELLS
    config = bcl2fastq_pipeline.getConfig.getConfig()
    flowcells_processed = pd.read_csv(os.path.join(config.get("FlowCellManager","managerDir"),'flowcells.processed'))
    return flowcells_processed.loc[flowcells_processed['flowcell_path'] == flowcell]

def pretty_print(df):
    print("Project \t Flowcell path \t Timestamp \t Archived")
    for i, row in df.iterrows():
        print("{}\t{}\t{}\t{}".format(row['project'], row['flowcell_path'], row['timestamp'], row['archived']))


if __name__=='__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers()
    parser_add = subparsers.add_parser("add",help="Add a project to the inventory file.")
    parser_add.set_defaults(func=add_flowcell)
    parser_add.add_argument("project",type=str,help="GCF project number.")
    parser_add.add_argument("path",type=str,help="Flowcell path.")
    parser_add.add_argument("timestamp",type=datetime.datetime.fromisoformat,help="A datetime isoformat string.")

    parser_archive = subparsers.add_parser("archive",help="Archive the flowcell by deleting fastq files and .7za archives.")
    parser_archive.set_defaults(func=archive_flowcell)
    parser_archive.add_argument("flowcell",type=str,help="Path to flowcell to be archived.")
    parser_archive.add_argument("--force",action="store_true",help="Force archive (no prompt)")

    parser_rerun = subparsers.add_parser("rerun",help="Rerun the flowcell by deleting the output directory.")
    parser_rerun.set_defaults(func=rerun_flowcell)
    parser_rerun.add_argument("flowcell",type=str,help="Path to flowcell to be deleted.")
    parser_rerun.add_argument("--force",action="store_true",help="Force archive (no prompt)")

    parser_list = subparsers.add_parser("list",help="List all flowcells.")
    parser_list.set_defaults(func=list_all,print_res=True)

    parser_list_processed = subparsers.add_parser("list-processed",help="List only flowcells in the inventory file processed by bfq pipeline.")
    parser_list_processed.set_defaults(func=list_processed,print_res=True)

    args = parser.parse_args()
    if vars(args).get("print_res",False):
        pretty_print(args.func(**vars(args)))
    else:
        args.func(**vars(args))

