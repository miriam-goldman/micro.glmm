# Parge MIDASDB gene_info.txt to extract core/accessory gene
# Chunyu Zhao 2023-02-16
# updated by Miriam 2023-10-16

import os
import sys
import argparse
from collections import defaultdict


"""
Input:
    s3://microbiome-pollardlab/uhgg_v1/pangenomes/{species_id}/gene_info.txt.lz4: one gene_id per row
    s3://microbiome-pollardlab/uhgg_v1/pangenomes/{species_id}/genes.len.lz4: one gene_id per row; gene_id, genome_id, gene_len
    UHGG (genomes.tsv)
Output: per-species gene_occurrence.tsv file
Methods:
    Given a list of genomes for given species, we first extracted all the
    genes from those genomes; Then counted the corresponding
    centroid_95 (by default) of these genes, andcompute the genome occurrence;
    Then for each ceontrid_95, we computed the prevalence, as the occurrence / total_genomes.

Example Usage:
    python parse_gene.py --species_id 102478 --pandb_dir /Users/chunyu.zhao/Desktop/2023-02-16-parse_pangene/midasdb
"""

def read_species_tsv(species_file):
    list_of_species = list()
    list_of_rep= list()
    with open(species_file) as instream:
        for line in instream:
            list_of_species.append(line.rstrip().split('\t')[1])
            list_of_rep.append("GUT_GENOME"+line.rstrip().split('\t')[2][10:16])
    return dict(zip(list_of_species,list_of_rep))

def read_genelen(glen_file):
    glen_dict = defaultdict(int)
    all_glen_dict = defaultdict(int)
    with open(glen_file) as instream:
        for line in instream:
            gene_id, genome_id, gene_len = line.rstrip().split('\t')
            if genome_id not rep_genome:
                glen_dict[gene_id] = gene_len
            all_glen_dict[gene_id] = gene_len
    return glen_dict, all_glen_dict


def read_geneinfo(ginfo_file, list_of_genes, by_col, rep_genome):
    centroid_counter = defaultdict(list)
    centroid_to_rep = defaultdict()
    with open(ginfo_file) as instream:
        header = instream.readline().strip('\n').split('\t')
        index_col = header.index(by_col)
        for line in instream:
            line = line.rstrip().split('\t') #<- read in each line of the gene_info.txt

            gene_id = line[0]
            centroid_id = line[index_col]
            genome_id = gene_id.split("_")[0]
            if gene_id in list_of_genes:
                centroid_counter[centroid_id].append(genome_id)

            # The FULL map for gene from representative genome to centroid_95
            if genome_id == rep_genome:
                centroid_to_rep[gene_id] = centroid_id
    return centroid_counter, centroid_to_rep


def compute_prevalence(centroid_counter, total_genomes, all_glen_dict, cutoff):
    centroid_occurence = defaultdict(dict)
    #centroid_matrix = defaultdict(lambda: defaultdict(int))
    centroid_matrix = defaultdict(dict)
    for cg, cl in centroid_counter.items():
        prev = len(set(cl)) / total_genomes
        centroid_occurence[cg]["gene_occurence"] = prev
        centroid_occurence[cg]["gene_length"] = all_glen_dict[cg]
        # Only write strictly non-core genes to martrix
        if prev < cutoff:
            set_of_genomes = set(cl)
            for _ in set_of_genomes:
                centroid_matrix[cg][_] = 1
    return centroid_occurence, centroid_matrix


def write_tofile(out_file, centroid_occurence, by_col):
    with open(out_file, "w") as ostream:
        ostream.write('\t'.join([by_col, "centroid_prevalence", "centroid_ength"]) + '\n')
        for centroid, cd in centroid_occurence.items():
            o, l = cd.values()
            o = format(o, ".3f")
            ostream.write('\t'.join([centroid, o, l]) + '\n')


def write_matfile(mat_file, centroid_matrix, by_col):
    with open(mat_file, "w") as ostream:
        ostream.write('\t'.join([by_col, "genome_id"]) + '\n')
        for k, d in centroid_matrix.items():
            gn = list(d.keys())[0]
            ostream.write('\t'.join([k, gn]) + '\n')


def write_repgene(rep_file, centroid_to_rep, by_col):
    with open(rep_file, "w") as ostream:
        ostream.write('\t'.join(["rep_gene_id", by_col]) + '\n')
        for g, c in centroid_to_rep.items():
            ostream.write('\t'.join([g, c]) + '\n')


def main():
    p = argparse.ArgumentParser(prog="python parse_gene.py", description='compute gene occurrence for given species.')
    p.add_argument(
        "--species_id", type=str, required=True,
        help="species_id")
    p.add_argument(
        "--pandb_dir", required=True,
        help="Path to MIDASDB pangenome base directory")
    p.add_argument(
        '--by_col', type=str, default="centroid_95",
        help=f"the centroid column to aggregate occurrence on.")
    p.add_argument(
        '--core_cutoff', type=float, default="0.9",
        help=f"The cutoff for define core gene based on prevalence.")
    p.add_argument(
        '--genome_info_file', type=str,
        help=f"File that contains the reference noted for each species")
    args = p.parse_args()
    by_col = args.by_col
    ginfo_file = f"{args.pandb_dir}/{args.species_id}/gene_info.txt"
    glen_file = f"{args.pandb_dir}/{args.species_id}/genes.len"
    out_file = f"{args.pandb_dir}/{args.species_id}/centroid_prevalence.tsv"
    mat_file = f"{args.pandb_dir}/{args.species_id}/centroid_matrix.tsv"
    rep_file = f"{args.pandb_dir}/{args.species_id}/centroid_to_repgenes.tsv"
    toc=read_species_tsv(args.genome_info_file)
    rep_genome = toc[args.species_id]
    core_cutoff = args.core_cutoff
    print(f"The representative genome is {rep_genome}")

    list_of_genomes = read_genomes(genome_file)
    total_genomes = len(set(list_of_genomes))
    print(f"Total number of genomes: {total_genomes}")

    glen_dict, all_glen_dict = read_genelen(glen_file, list_of_genomes)
    # We only care about genes from the specified list of genomes
    list_of_genes = set(glen_dict.keys())
    total_genes = len(list_of_genes)
    print(f"Total number of pan genes: {total_genes}")

    centroid_counter, centroid_to_rep = read_geneinfo(ginfo_file, list_of_genes, by_col, rep_genome)
    centroid_occurence, centroid_matrix = compute_prevalence(centroid_counter, total_genomes, all_glen_dict, core_cutoff)


    write_tofile(out_file, centroid_occurence, by_col)
    write_matfile(mat_file, centroid_matrix, by_col)
    write_repgene(rep_file, centroid_to_rep, by_col)

main()
