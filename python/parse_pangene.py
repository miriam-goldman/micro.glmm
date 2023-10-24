# Parge MIDASDB gene_info.txt to extract core/accessory gene
# Chunyu Zhao 2023-06-20

import os
import sys
import argparse
import pandas as pd
import numpy as np


"""
Input:
    {midasdb}/pangenomes/{species_id}/gene_info.txt
    {midasdb}/pangenomes/{species_id}/genes.len
    {midasdb}/genomes.tsv
Output: three TSV files
    - centroids_to_repgenes: by_col to only representative genes
    - centroids_prevalence: by_col, prevalence, gene_length
    - centroid_matrix: by_col - by - genome presence/absence matrix
Methods:
    prevalence := occurrence / total_genomes
Example Usage:
    python parse_pangene.py --species_id 102478 --midasdb /path/to/local-midasddb --genome_list /path/to/list_of_genomes --outdir /path/to/output/directory
"""

def main():
    p = argparse.ArgumentParser(prog="python parse_pangene.py", description='Compute pan-gene prevalaence for given species.')
    p.add_argument(
        "--species_id", type=str, required=True,
        help="species_id")
    p.add_argument(
        "--midasdb", required=True,
        help="Path to MIDASDB")
    p.add_argument(
        "--genome_list", required=True,
        help="list of candidate genomes")
    p.add_argument(
        "--output_dir", required=True,
        help="path to base output directory")
    p.add_argument(
        '--by_col', type=str, default="centroid_95",
        help=f"the centroid column to aggregate occurrence on.")
    p.add_argument(
        '--core_cutoff', type=float, default="0.9",
        help=f"The cutoff for define core gene based on prevalence.")

    args = p.parse_args()

    species_id = args.species_id
    midasdb = args.midasdb
    by_col = args.by_col
    genome_file = args.genome_list
    outdir = args.output_dir
    core_cutoff = args.core_cutoff

    pangenome_db = f"{midasdb}/pangenomes/{species_id}"

    outdir = f"{outdir}/{species_id}"
    os.makedirs(f"{outdir}", exist_ok=True)

    # read list of genomes
    genome_df = pd.read_csv(genome_file, delimiter='\t', names=['genome_id'])
    total_genomes = len(genome_df)

    # read MIDASDB TOC and get representative genome of given species_id
    toc = pd.read_csv(f"{midasdb}/genomes.tsv", delimiter="\t")
    my_rep_genome = toc[(toc['species'] == int(species_id)) & (toc['genome_is_representative'] == 1)]
    my_rep_genome = my_rep_genome.iloc[0,2]

    is_rep_in = genome_df['genome_id'].isin([my_rep_genome]).any()
    if not is_rep_in:
        print("Representative genome is not in the given list of genomes.")

    # read gene_info
    gene_info = pd.read_csv(f"{pangenome_db}/gene_info.txt", delimiter="\t")

    # add genome_id column
    gene_info['genome_id'] = gene_info['gene_id'].replace(to_replace=r'_(.*)', value='', regex=True)
    gene_info['genome_id'] = gene_info['genome_id'].replace('UHGG', 'GUT_GENOME', regex=True)

    # read genes.len
    gene_len_df = pd.read_csv(f"{pangenome_db}/genes.len", delimiter="\t", names=['gene_id', "genome_id", "gene_length"])

    # write rep-gene to centroid-by-col to centroid_to_repgenes
    gene_to_rep = gene_info[gene_info['genome_id'] == my_rep_genome]
    gene_to_rep[['gene_id', by_col]].to_csv(f"{outdir}/{by_col}_to_repgene.tsv", sep='\t', index=False)

    # subset genes_info by only list_of_genomes
    gene_info = pd.merge(gene_info, genome_df, on='genome_id')
    gene_info = gene_info[[by_col, 'genome_id']].drop_duplicates()

    # compute prevalence
    prevalence = gene_info.groupby(by_col)['genome_id'].count().reset_index(name='prevalence')
    prevalence['prevalence'] = prevalence['prevalence'] / total_genomes
    # add centorid_length directly from the gene_id
    prevalence = pd.merge(prevalence, gene_len_df, left_on=[by_col], right_on=['gene_id'])
    prevalence = prevalence[[by_col, 'prevalence', 'gene_length']]
    # write all genes' prevalence to centroid_prevalence.tsv
    prevalence.columns = [by_col, "centroid_prevalence", "centroid_length"]
    prevalence.to_csv(f"{outdir}/{by_col}_prevalence.tsv", sep='\t', index=False)

    # filter out core genes
    prevalence = prevalence[prevalence['centroid_prevalence'] <=  core_cutoff]
    # subset gene_info by only centroids meeting the above non-core filters
    gene_info = pd.merge(gene_info, prevalence, on=by_col)
    gene_info = gene_info[[by_col, 'genome_id']]
    gene_info.to_csv(f"{outdir}/{by_col}_matrix.tsv", sep='\t', index=False)

    # convert long dataframe to wide matrix
    #pivot_df = gene_info.groupby([by_col, 'genome_id']).size().unstack(fill_value=0).astype(int)
    #pivot_df.to_csv(f"{outdir}/{by_col}_matrix.tsv", sep='\t', index=False)

main()
