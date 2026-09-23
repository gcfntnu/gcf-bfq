import multiprocessing as mp
import os
import shutil
import subprocess
import tempfile

from pathlib import Path

import pandas as pd

forward = pd.read_csv("qiaseq_primers_fwd.csv", index_col=0)
reverse = pd.read_csv("qiaseq_primers_rev.csv", index_col=0)


def remove_region_marker(path):
    """Apply the former sed replacement without interpreting the path as shell input."""
    path = Path(path)
    with path.open() as source, tempfile.NamedTemporaryFile(
        mode="w", dir=path.parent, delete=False
    ) as destination:
        for line in source:
            destination.write(line.replace(":region=no_adapter", ""))
        temporary_path = Path(destination.name)
    temporary_path.replace(path)


def compress_fastqs(input_paths, output_path):
    """Stream several FASTQs through pigz into one gzip file."""
    output_path = Path(output_path)
    with output_path.open("wb") as output_fh:
        process = subprocess.Popen(
            ["pigz", "-6", "-p", "8"],
            stdin=subprocess.PIPE,
            stdout=output_fh,
        )
        try:
            for input_path in input_paths:
                with Path(input_path).open("rb") as input_fh:
                    shutil.copyfileobj(input_fh, process.stdin)
            process.stdin.close()
            returncode = process.wait()
        except Exception:
            process.terminate()
            process.wait()
            raise
    if returncode:
        raise subprocess.CalledProcessError(returncode, process.args)


def cutadapt_worker(fname):
    fname = Path(fname)
    sample = fname.name.removesuffix("_R1.fastq.gz")
    regions = ["unknown"]
    temporary_inputs = set()

    for region, row in forward.iterrows():
        unknown_r1 = Path(f"{sample}_unknown_R1.fastq")
        unknown_r2 = Path(f"{sample}_unknown_R2.fastq")
        if unknown_r1.exists():
            input_r1 = Path(f"input_{unknown_r1.name}")
            input_r2 = Path(f"input_{unknown_r2.name}")
            unknown_r1.replace(input_r1)
            unknown_r2.replace(input_r2)
            remove_region_marker(input_r1)
            remove_region_marker(input_r2)
            temporary_inputs.update((input_r1, input_r2))
        else:
            input_r1 = fname
            input_r2 = fname.with_name(fname.name.replace("R1.fastq", "R2.fastq"))

        cmd = [
            "cutadapt",
            "-g",
            f"{region}={row['primer']}",
            "-G",
            f"{region}={reverse.loc[region, 'primer']}",
            "--pair-adapters",
            "--no-indels",
            "-e",
            "0.1",
            "--untrimmed-output",
            str(unknown_r1),
            "--untrimmed-paired-output",
            str(unknown_r2),
            "--suffix",
            ":region={name}",
            "-o",
            f"{sample}_{{name}}_R1.fastq",
            "-p",
            f"{sample}_{{name}}_R2.fastq",
            str(input_r1),
            str(input_r2),
        ]
        regions.append(region)
        with Path("log", f"{sample}_qiaseq_demultiplex.log").open("ab") as log_fh:
            subprocess.check_call(cmd, stdout=log_fh)

    for input_path in temporary_inputs:
        input_path.unlink(missing_ok=True)

    r1_paths = [Path(f"{sample}_{region}_R1.fastq") for region in regions]
    r1_paths = [path for path in r1_paths if path.exists()]
    r2_paths = [path.with_name(path.name.replace("R1.fastq", "R2.fastq")) for path in r1_paths]

    compress_fastqs(r1_paths, f"{sample}_R1.fastq.gz")
    compress_fastqs(r2_paths, f"{sample}_R2.fastq.gz")

    for path in [*r1_paths, *r2_paths]:
        path.unlink(missing_ok=True)
    print(f"finished sample: {sample}")


def main():
    os.makedirs("log", exist_ok=True)
    r1_fastqs = list(Path("data").glob("*R1.fastq.gz"))
    with mp.Pool(32) as pool:
        pool.map(cutadapt_worker, r1_fastqs)


if __name__ == "__main__":
    main()
