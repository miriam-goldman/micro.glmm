import os
import sys
import pandas as pd
import re
from pathlib import Path
from collections import Counter, OrderedDict


def parse_trim_summary_paired(f):
    for line in f.readlines():
        if line.startswith('Input Read'):
            vals = re.findall('\D+\: (\d+)', line)
            keys = ('input', 'both_kept','fwd_only','rev_only','dropped')
            return(OrderedDict(zip(keys, vals)))


def parse_trim_summary_unpaired(f):
    for line in f:
        if line.startswith('Input Read'):
            vals = re.findall('\D+\: (\d+)', line)
            keys = ('input', 'both_kept', 'dropped')
            return(OrderedDict(zip(keys, vals)))


def parse_bbduk_summary(f, paired_end):
    for line in f.readlines():
        if line.startswith("Low entropy discards"):
            vals = re.findall('(\d+)', line)[0]
            keys = ('low_entropy')
            if paired_end:
                vals = float(vals)/2
            ret = OrderedDict()
            ret['low_entropy'] = vals
            return ret


def summarize_bbduk(dfile, paired_end):
    dname = os.path.basename(dfile).split('.log')[0]
    with open(dfile) as df:
        bbduk_data = parse_bbduk_summary(df, paired_end)
    bbduk_data['sample'] = dname
    bbduk_data.move_to_end('sample', last=False)
    #sys.stderr.write("trim data: {}\n".format(trim_data))
    return pd.DataFrame(bbduk_data, index=[dname])


def summarize_trim(tfile, paired_end):
    tname = os.path.basename(tfile).split('.out')[0]
    with open(tfile) as tf:
        if paired_end:
            trim_data = parse_trim_summary_paired(tf)
        else:
            trim_data = parse_trim_summary_unpaired(tf)
    trim_data['sample'] = tname
    trim_data.move_to_end('sample', last=False)
    #sys.stderr.write("trim data: {}\n".format(trim_data))
    return pd.DataFrame(trim_data, index=[tname])


def merge_decontam(files):
    df = pd.read_csv(files[0], sep='\t')
    for file in files[1:]:
        df_next = pd.read_csv(file, sep='\t')
        df = pd.merge(df, df_next, on=df.columns[0])
    return df


paired_end = snakemake.config["all"]["paired_end"]
trim_list = [summarize_trim(tf, paired_end) for tf in snakemake.input.trim_files]
trim_df = pd.concat(trim_list)

bbduk_list = [summarize_bbduk(df, paired_end) for df in snakemake.input.bbduk_files]
bbduk_df = pd.concat(bbduk_list)

decontam_df = merge_decontam(snakemake.input.decontam_files)

df = trim_df.merge(bbduk_df, on='sample').merge(decontam_df, on='sample')

for col in df.columns[1:]:
    df[col] = df[col].astype(int)
df['non_host'] = df['both_kept'] - df['host']

df.to_csv(snakemake.output[0], sep="\t", index=False)
