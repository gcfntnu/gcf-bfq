import os
import subprocess
import multiprocessing as mp
import glob
import pandas as pd

forward = pd.read_csv("qiaseq_primers_fwd.csv", index_col=0)
reverse = pd.read_csv("qiaseq_primers_rev.csv", index_col=0)

def cutadapt_worker(fname):
    sample = os.path.basename(fname).replace("_R1.fastq.gz","")
    for i, r in forward.iterrows():
        if os.path.exists("{}_unknown_R1.fastq".format(sample)):
            fname = "{}_unknown_R1.fastq".format(sample)
            cp_cmd = "mv -f {} input_{}".format(fname, fname)
            subprocess.check_call(cp_cmd, shell=True)
            cp_cmd = "mv -f {} input_{}".format(fname.replace("R1.fastq","R2.fastq"), fname.replace("R1.fastq","R2.fastq"))
            subprocess.check_call(cp_cmd, shell=True)
            sed_cmd = "sed -i -e s/:region=no_adapter//g input_{}".format(fname)
            subprocess.check_call(sed_cmd, shell=True)
            sed_cmd = "sed -i -e s/:region=no_adapter//g input_{}".format(fname.replace("R1.fastq","R2.fastq"))
            subprocess.check_call(sed_cmd, shell=True)


        cmd = "cutadapt -g {region}={fwd_primer} -G {region}={rev_primer} --pair-adapters --no-indels -e 0.1 --untrimmed-output {unknown_r1} --untrimmed-paired-output {unknown_r2} --suffix ':region={{name}}' -o {sample}_{{name}}_R1.fastq -p {sample}_{{name}}_R2.fastq {r1} {r2} >> log/{sample}_qiaseq_demultiplex.log".format(
            sample = sample,
            unknown_r1 = "{}_unknown_R1.fastq".format(sample),
            unknown_r2 = "{}_unknown_R2.fastq".format(sample),
            region = i,
            fwd_primer = r['primer'],
            rev_primer = reverse.loc[i, 'primer'],
            r1 = ("input_" + fname) if "unknown_R1.fastq" in fname else fname,
            r2 = ("input_" + fname.replace("R1.fastq", "R2.fastq")) if "unknown_R1.fastq" in fname else fname.replace("R1.fastq", "R2.fastq")
            )
        subprocess.check_call(cmd, shell=True)

    rm_cmd = "rm input_{}_*fastq".format(sample)
    subprocess.check_call(rm_cmd, shell=True)




r1 = glob.glob(os.path.join("data","*R1.fastq.gz"))

p = mp.Pool(4)
p.map(cutadapt_worker, r1)
p.close()
p.join()
