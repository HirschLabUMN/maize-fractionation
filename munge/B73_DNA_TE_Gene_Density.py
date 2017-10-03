#!/usr/bin/env python
"""Calculate gene desnity, DNA TE density, and RNA TE density in sliding windows
across the B73 genome. This will be used to make a big summary plot for later.
Takes four arguments:
    1) Genes GFF
    2) DNA TE GFF
    3) RNA TE GFF
    4) Reference .fai file (from samtools faidx)
The default window size is 1Mb, and the default window step is 500kb.
"""

import sys


def get_chrom_len(fai):
    """Read the fasta index (.fai. from samtools faidx) to get the lengths of
    the chromosomes, in bases."""
    chrom_lengths = {}
    with open(fai, 'r') as f:
        for line in f:
            tmp = line.strip().split('\t')
            chrom_lengths[tmp[0]] = int(tmp[1])
    return chrom_lengths


def parse_gff(gff):
    """Parse a GFF and return a data structure of the following format:
    {
        chr: [gene1_start, gene2_start, ...],
        chr: ...
    }

    We only care about the start positions, since we going to be counting in
    intervals, and not doing any sophisticated operations on them."""
    gff_data = {}
    with open(gff, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                tmp = line.strip().split('\t')
                if tmp[0] in gff_data:
                    gff_data[tmp[0]].append(int(tmp[3]))
                else:
                    gff_data[tmp[0]] = [int(tmp[3])]
    return gff_data


def generate_windows(start, stop, size, step):
    """Build a generator to return endpoints for a sliding window."""
    curr_pos = start
    while curr_pos < stop:
        if curr_pos + size < stop:
            yield (curr_pos, curr_pos + size)
            curr_pos += step
        else:
            yield (curr_pos, stop)
            curr_pos += stop


def main(genes, dna_te, rna_te, fai, winsize=1000000, step=500000):
    """Main function."""
    c_lengths = get_chrom_len(fai)
    gene_pos = parse_gff(genes)
    dna_te_pos = parse_gff(dna_te)
    rna_te_pos = parse_gff(rna_te)
    #   Print a header
    print 'Chromosome\tMidpoint\tGenesPerMb\tDNATEsPerMb\tRNATEsPerMb'
    #   Sort chromosomes, and iterate
    for c in sorted(c_lengths):
        #   Generate windows
        windows = generate_windows(1, c_lengths[c], winsize, step)
        for w in windows:
            #   Count how many features are in the window
            gene_in_win = [
                True
                if x >= w[0] and x <= w[1]
                else False
                for x
                in gene_pos[c]
                ]
            #   Same for DNA TEs
            dna_te_in_win = [
                True
                if x >= w[0] and x <= w[1]
                else False
                for x
                in dna_te_pos[c]
                ]
            rna_te_in_win = [
                True
                if x >= w[0] and x <= w[1]
                else False
                for x
                in rna_te_pos[c]
                ]
            #   Build the string to print
            toprint = '\t'.join([
                c,
                str((w[0]+w[1])/2.0),
                str((float(sum(gene_in_win)) / (w[1] - w[0])) * winsize),
                str((float(sum(dna_te_in_win)) / (w[1] - w[0])) * winsize),
                str((float(sum(rna_te_in_win)) / (w[1] - w[0])) * winsize)
                ]
                )
            print toprint
    return


main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
