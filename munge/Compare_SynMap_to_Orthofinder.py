#!/usr/bin/env python
"""Compare the acnestral locus assignment from two approaches. The first was
multiple sequence alignment of gene groups identifed through Maize-ancestral
BLAST searches through SynMap in CoGe. The second was orthologous group
identification with publicly available grass species, then multiple
sequence alignment of maize-ancestral genes. Because the SynMap searches were
done with an older version of the Sorghum assembly, we also need to translate
the Sb... numbers to Sobic... numbers. Takes four arguments:
    1) Sorghum name translation from Phytozome
    2) List of tandem duplicate IDs and genes in them
    3) SynMap assignment table
    4) Orthofinder assignment table
"""

import sys
import pprint


def parse_translation(key_file):
    """Parse the Phytozome name file, and store it in a dictionary."""
    sb_key = {}
    with open(key_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            else:
                tmp = line.strip().split('\t')
                # If the number of elements in the line is not 4, then the gene
                # is unmapped in the new assembly. We cannot do anything about
                # this.
                if len(tmp) != 4:
                    continue
                else:
                    # We will key on the old ID
                    sb_key[tmp[3]] = tmp[2]
    return sb_key


def parse_clusters(cluster_key):
    """Make a dictionary that is keyed on the cluster ID, and a list of genes
    as the value."""
    cluster = {}
    with open(cluster_key, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            cluster[tmp[0]] = tmp[1].split(',')
    return cluster


def parse_synmap(synmap, t_key):
    """Parse the SynMap summary table and keep it in a dictionary. Be sure to
    convert all old Sorghum gene IDs into new ones, if possible. If there
    multiple genes in the "Best_Gene" column, ignore it, since it would not have
    been assigned anyway."""
    syn_matches = {}
    with open(synmap, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split('\t')
                # slice off the first four characters from the cluster ID, since
                # the SynMap table has 'Anc_' in front of it
                cluster = tmp[0][4:]
                out_gene = tmp[4]
                mz_genes = tmp[5].split(',')
                if out_gene.startswith('S'):
                    if out_gene in t_key:
                        real_out = t_key[out_gene]
                    else:
                        real_out = out_gene
                else:
                    real_out = out_gene
                syn_matches[cluster] = [(m, real_out) for m in mz_genes]
    return syn_matches


def parse_orthofinder(orthofinder):
    """Parse the Orthofinder summary table and keep it in a dictionary. Keep
    track of the rice genes, too."""
    ortho_matches = {}
    with open(orthofinder, 'r') as f:
        for index, line in enumerate(f):
            if index == 0:
                continue
            else:
                tmp = line.strip().split('\t')
                cluster = tmp[0]
                mz_gene = tmp[2]
                sb_gene = tmp[3]
                os_gene = tmp[5]
                # Some clusters have multiple assignments - keep track of these!
                if cluster in ortho_matches:
                    ortho_matches[cluster].append((mz_gene, sb_gene, os_gene))
                else:
                    ortho_matches[cluster] = [(mz_gene, sb_gene, os_gene)]
    return ortho_matches


def compare_tables(ck, sm, og):
    """Compare the SynMap and Orthofinder outputs, and print which match and
    which mismatch. Also summarize which relationships could not be tested due
    to gene IDs not mapping."""
    # We keep track of a big matrix to print out. The first row is the header.
    to_print = []
    for c in sorted(ck):
        if c in sm and c in og:
            # If the cluster is both in SynMap BLAST and Orthofinder, then we
            # pull out the relationships from each method
            synmap = sm[c]
            ortho = og[c]
            # Set some default values for Sb and Os genes
            sb_gene = 'NA'
            os_gene = 'NA'
            # We want to have pairwise comparisons, so we decompose the
            # potentially multi-element groupings into pairwise ones
            ortho_sb = [(s[0], s[1]) for s in ortho]
            ortho_os = [(s[0], s[2]) for s in ortho]
            # We now have a list of pairwise relationships between maize
            # genes and ancestral genes from Sorghum and/or rice. We want to
            # identify overlapping relationships.
            rep_v_sb = list(set(synmap) & set(ortho_sb))
            rep_v_os = list(set(synmap) & set(ortho_os))
            # Determine the representative gene. If rep_v_sb exists, then there
            # is a Sorghum match. If rep_v_os exists, then there is a rice
            # match. We prefer Sorghum to rice.
            if rep_v_sb:
                rep_gene = rep_v_sb[0][0]
                syn = False
            elif rep_v_os:
                rep_gene = rep_v_os[0][0]
                syn = False
            else:
                # In the case of a mismatch, we take the acenstral gene from the
                # SynMap BLAST
                syn = True
                rep_gene = synmap[0][0]
                # Check the outgroup gene. If it starts with 'Sobic' then it's
                # Sorghum.
                anc_gene = synmap[0][1]
                if anc_gene.startswith('Sobic'):
                    sb_gene = anc_gene
                elif anc_gene.startswith('OS'):
                    os_gene = anc_gene
            # We also want to keep track of the Sorghum gene and the rice gene
            # that were chosen by Orthofinder
            if rep_v_sb:
                sb_gene = rep_v_sb[0][1]
            if rep_v_os:
                os_gene = rep_v_os[0][1]
            # Append a row to the matrix to print, only if the rep gene
            # is not NA
            if sb_gene == 'NA' and os_gene == 'NA':
                continue
            elif syn:
                to_print.append([c, ','.join(ck[c]), rep_gene, sb_gene,
                                 os_gene, "SynMap_Chosen"])
            elif not syn:
                to_print.append([c, ','.join(ck[c]), rep_gene, sb_gene,
                                 os_gene, "Both"])
            # If both sb_gene and os_gene are NA, then we have a mapping
            # failure
            # Then, for the rest of the Otherfinder clusters, we append them
            # to the table to print
            for other_relationship in og[c]:
                if other_relationship[0] == rep_gene:
                    continue
                else:
                    to_print.append([c, ','.join(ck[c]), other_relationship[0],
                                     other_relationship[1],
                                     other_relationship[2], "Other_Orthologue"])
        if c in sm and c not in og:
            # If we get here, then the Relationship was found in SynMap BLAST,
            # and not in Orthofinder. These will also be the only relationships
            # that have Sorghum assembly mapping failures.
            synmap = sm[c]
            # First, check the format of the ancestral gene ID. If it is
            # a Sb0... ID, then we cannot use it.
            anc = synmap[0][1]
            if anc.startswith('Sb'):
                continue
            else:
                # Then, check the length of the SynMap BLAST list. if it is
                # greater than 1, then there is more than one dS relationship
                # identified, and we cannot choose a representative gene.
                if len(synmap) > 1:
                    mz = sorted([s[0] for s in synmap])[0]
                    # The ancestral locus should be the same for all elements of
                    # this list, so we just take the first one
                    anc = synmap[0][1]
                    # If the ancestral locus starts with Sobic, then it is a
                    # Sorghum gene. Else, it is a rice gene.
                    if anc.startswith('Sobic'):
                        sb_gene = anc
                        os_gene = 'NA'
                    else:
                        sb_gene = 'NA'
                        os_gene = anc
                    # We will take the leftomost gene, per convention
                    to_print.append([c, ','.join(ck[c]), mz, sb_gene, os_gene, 
                                     'BLAST_Leftmost'])
                else:
                    # For the cases where SynMap BLAST identifies one winner
                    # and does not have a mapping failure, just print it out
                    mz = synmap[0][0]
                    anc = synmap[0][1]
                    if anc.startswith('Sobic'):
                        sb_gene = anc
                        os_gene = 'NA'
                    else:
                        sb_gene = 'NA'
                        os_gene = anc
                    to_print.append([c, ','.join(ck[c]), mz, sb_gene, os_gene,
                                     'BLAST'])
        if c not in sm and c in og:
            # If we get here, then the cluster has Orthofinder matches, but does
            # not have SynMap BLAST matches. This is easy, we just print it out
            # as usual
            for rel in og[c]:
                to_print.append([c, ','.join(ck[c]), rel[0], rel[1], rel[2],
                                 'Orthofinder_Only'])
    return to_print


def main(key, cluster_key, synmap, orthofinder):
    """Main function."""
    sb_key = parse_translation(key)
    clust_key = parse_clusters(cluster_key)
    synmap_table = parse_synmap(synmap, sb_key)
    ortho_table = parse_orthofinder(orthofinder)
    s_o_comp = compare_tables(clust_key, synmap_table, ortho_table)
    # Print out the summary table
    print '\t'.join([
        'Cluster_ID',
        'Cluster_Genes',
        'Maize_Gene',
        'Sorghum_Gene',
        'Orzya_Gene',
        'Note'])
    for row in s_o_comp:
        print '\t'.join(row)
    return


if len(sys.argv) != 5:
    print """
Compare the acnestral locus assignment from two approaches. The first was
multiple sequence alignment of gene groups identifed through Maize-ancestral
BLAST searches through SynMap in CoGe. The second was orthologous group
identification with publicly available grass species, then multiple
sequence alignment of maize-ancestral genes. Because the SynMap searches were
done with an older version of the Sorghum assembly, we also need to translate
the Sb... numbers to Sobic... numbers.  Takes four arguments:
    1) Sorghum name translation from Phytozome
    2) List of tandem duplicate IDs and genes in them
    3) SynMap assignment table
    4) Orthofinder assignment table"""
    exit(1)

main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
