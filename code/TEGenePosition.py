#### te.pipeline 2025-09-19
#### TE-Gene position
#### Based on the TE annotation file and the gene annotation file
import argparse
import os
import re

import pandas as pd

def read_gtf(file, feature_type=None):
    df = pd.read_csv(file, sep='\t', header=None,
                     names=['chromosome', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
    if feature_type:
        df = df[df['feature'] == feature_type]

    # Grouped by chromosomes
    grouped = df.groupby('chromosome')
    return {chrom: grouped.get_group(chrom) for chrom in grouped.groups}

def check_te_in_gene(te, gene, upstream_downstream=40_000):
    gene_start, gene_end, gene_strand = gene['start'], gene['end'], gene['strand']
    te_start, te_end = te['start'], te['end']

    if te_end < gene_start - upstream_downstream:
        return None
    elif te_start >= gene_start and te_end <= gene_end:
        return 0
    elif te_start < gene_start and te_end >= gene_start - upstream_downstream:
        distance = int(gene_start - te_end)
        if distance < 0:
            distance = 0
        if gene_strand == '+':
            return distance
        else:
            return -distance
    elif te_start <= gene_end + upstream_downstream and te_end > gene_end:
        distance = int(te_start - gene_end)
        if distance < 0:
            distance = 0
        if gene_strand == '+':
            return -distance
        else:
            return distance
    else:
        return None


def analyze_te_insertions(gene_file, te_file, feature_type, upstream_downstream=40_000):
    genes_by_chrom = read_gtf(gene_file, feature_type=feature_type)
    tes_by_chrom = read_gtf(te_file)

    results = []

    for chrom, genes in genes_by_chrom.items():
        if chrom in tes_by_chrom:
            tes = tes_by_chrom[chrom]
            te_index = 0

            for _, gene in genes.iterrows():
                gene_name = str(gene['attribute'])

                # Eliminate TE whose starting position is located outside the upstream range of gene
                while te_index < len(tes) and tes.iloc[te_index]['end'] < gene['start'] - upstream_downstream:
                    te_index += 1

                # Start checking
                current_te_index = te_index

                while current_te_index < len(tes):
                    te = tes.iloc[current_te_index]
                    te_name = str(te['attribute'])

                    # starting position of TE is beyond the downstream region of gene, stop the inspection
                    if te['start'] > gene['end'] + upstream_downstream:
                        break

                    # position of TE and gene
                    position = check_te_in_gene(te, gene, upstream_downstream)
                    if position != 'None':
                        results.append({
                            'Gene_name': gene_name,
                            'Gene_start': int(gene['start']),
                            'Gene_end': int(gene['end']),
                            'Gene_strand': str(gene['strand']),
                            'TE_name': te_name,
                            'TE_start': int(te['start']),
                            'TE_end': int(te['end']),
                            'Chromosome': te['chromosome'],
                            'Position': str(position)
                        })

                    current_te_index += 1

    return pd.DataFrame(results)

if __name__ == '__main__':
    tool_name = 'TEGenePosition'
    version_num = '0.0.1'
    default_feature_type = 'gene'
    default_up_down_distance = 40_000

    # 1.parse args
    describe_info = '########################## ' + tool_name + ', version ' + str(version_num) + ' ##########################'
    parser = argparse.ArgumentParser(description=describe_info)
    parser.add_argument('--gene_path', required=True, metavar='gene_path',
                        help='Input the path of gene gtf.')
    parser.add_argument('--TE_gtf_path', required=True, metavar='TE_gtf_path',
                        help='Input the path of TE gtf.')
    parser.add_argument('--output_dir', required=False, metavar='output_dir',
                        help='Input the path of output directory')
    parser.add_argument('--feature_type', required=False, metavar='feature_type',
                        help='Input the types of genes (gene/mRNA). default=' + default_feature_type)
    parser.add_argument('--up_down_distance', required=False, metavar='up_down_distance',
                        help='Gene-centered upstream and downstream distance. default=' + str(default_up_down_distance))


    args = parser.parse_args()

    gene_path = args.gene_path
    TE_gtf_path = args.TE_gtf_path
    output_dir = args.output_dir
    feature_type = args.feature_type
    up_down_distance = args.up_down_distance

    if not os.path.isabs(gene_path):
        gene_path = os.path.abspath(gene_path)
    if not os.path.isabs(TE_gtf_path):
        TE_gtf_path = os.path.abspath(TE_gtf_path)

    if output_dir is None:
        output_dir = os.getcwd()
    else:
        output_dir = output_dir

    if not os.path.isabs(output_dir):
        output_dir = os.path.abspath(output_dir)

    if feature_type is None:
        feature_type = default_feature_type
    else:
        feature_type = feature_type

    if up_down_distance is None:
        up_down_distance = default_up_down_distance
    else:
        up_down_distance = int(up_down_distance)

    results = analyze_te_insertions(gene_path, TE_gtf_path, feature_type, up_down_distance)
    output_file = output_dir + '/TE_Gene_position.tsv'
    results.to_csv(output_file, sep='\t', index=False)