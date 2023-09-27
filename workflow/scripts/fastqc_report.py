import os
import sys
import re
import pandas
from io import StringIO


def parse_fastqc_quality(filename):
    with open(filename) as f:
        report = f.read()
    tableString = re.search(
        '\>\>Per base sequence quality.*?\n(.*?)\n\>\>END_MODULE',
        report, re.DOTALL).group(1)

    f_s = StringIO(tableString)
    df = pandas.read_csv(
        f_s, sep='\t', usecols=['#Base', 'Mean'], index_col='#Base')
    sample_name = os.path.basename(filename.split('_fastqc')[0])
    df.columns=[sample_name]
    f_s.close()
    return(df)


quality_list = [parse_fastqc_quality(file) for file in snakemake.input.files]
quality_table = pandas.concat(quality_list, axis=1).transpose()
quality_table.to_csv(snakemake.output[0], sep="\t", index_label="Samples")
