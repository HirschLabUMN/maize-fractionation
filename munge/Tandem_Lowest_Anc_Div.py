#!/usr/bin/env python
"""Choose the tandem duplicate that has the lowest dS from the outgroup as the
most "represenntative" for the false tandem duplicate clusters. Takes two
arguments:
    1) "False" tandem clusters file
    2) dN, dS, divergence summary table from yn00, with ancestral
"""

import sys
import pprint

try:
    tandemfile = sys.argv[1]
    divfile = sys.argv[2]
except IndexError:
    print 'Incorrect number of arguments supplied. The arguments should be'
    print '    Tandem Cluster File'
    print '    yn00 summary file, with ancestral seqs'
    exit(1)

tandems = []
with open(tandemfile, 'r') as f:
    for line in f:
        # We prepend 'Anc_' because that's how the clusters with ancestral
        # sequence are named.
        tandems.append('Anc_' + line.strip())

# Then, store a dictionary of the tandem clusters. There will be several
# filters applied to the clusters before they make it intto this dictionary:
#   1) They must be in the supplied clusters file
#   2) Only comparisons to the ancestral gene are listed
#   3) dS between the gene and the ancestral gene must be below 2.5
filtered_clusters = {}
gene_membership = {}
with open(divfile, 'r') as f:
    for index, line in enumerate(f):
        if index == 0:
            continue
        else:
            tmp = line.strip().split()
            cluster_id = tmp[0]
            genes = tmp[1]
            gene_membership[cluster_id] = genes.split(',')[1:]
            og = tmp[1].split(',')[0]
            gene_a = tmp[2]
            gene_b = tmp[3]
            dS = tmp[9]
            dSSE = tmp[10]
            # If the cluster ID is not one we wnated to analyze, or the second
            # gene in the comparison was maize, or dS is too high, we skip it
            if (cluster_id not in tandems) or (gene_b.startswith('Zm')) or (float(dS) + float(dSSE) > 2.5) or (dS == 'nan'):
                continue
            else:
                if cluster_id in filtered_clusters:
                    filtered_clusters[cluster_id][gene_a] = (og, float(dS), float(dSSE))
                else:
                    filtered_clusters[cluster_id] = {gene_a: (og, float(dS), float(dSSE))}


# Print a header
print '\t'.join([
    'Cluster_ID',
    'Cluster_Size',
    'Num_Lowest_Div',
    'Cluster_Genes',
    'Outgroup_Gene',
    'Best_Gene',
    'Outgroup_dS',
    'Outgroup_dS_SE'])
# Then, we want to apply some filters to the remaining clusters. 
for clust_id in sorted(filtered_clusters.keys()):
    # First, find out how many genes this cluster has
    numgenes = len(filtered_clusters[clust_id])
    num_tot = len(gene_membership[clust_id])
    # If it has just one gene, return that
    if numgenes == 1:
        bestgene = filtered_clusters[clust_id].keys()[0]
        bestdS = filtered_clusters[clust_id][bestgene][1]
        bestdSSE = filtered_clusters[clust_id][bestgene][2]
        outgroup = filtered_clusters[clust_id][bestgene][0]
        numbest = '1'
    else:
        # Else we have to choose the best one. We do this by the gene with the
        # lowest dS +/- the standard error. First, calculate the intervals for
        # each dS+-SE.
        intervals = []
        for d in filtered_clusters[clust_id].iteritems():
            # Append a tuple that has the form:
            #   (Maize_Gene, Outgroup_Gene, dS-SE, dS+SE)
            intervals.append((d[0], d[1][0], d[1][1]-d[1][2], d[1][1]+d[1][2]))
        # Then, sort the intervals by the lower bound of dS
        intervals = sorted(intervals, key=lambda x: x[2])
        # then, identify overlapping intervals with the first one
        overlapping = []
        for i in range(1, len(intervals)):
            base = intervals[0]
            comp = intervals[i]
            # If the comparison interval start is greater than the start of the
            # base interval, and less than the end, it is overlapping
            if comp[2] >= base[2] and comp[2] <= base[3]:
                overlapping.append(comp)
        # Then, if there are no overlapping genes, then choose the gene with the lowest mean
        if not overlapping:
            bestgene = intervals[0][0]
            outgroup = intervals[0][1]
            bestdS = (intervals[0][2]+intervals[0][3])/2
            bestdSSE = (intervals[0][3]-intervals[0][2])/2
            numbest = '1'
        else:
            bestgene = ','.join([g[0] for g in overlapping+[intervals[0]]])
            outgroup = filtered_clusters[clust_id][bestgene.split(',')[0]][0]
            bestdS = 'NA'
            bestdSSE = 'NA'
            numbest = str(len(overlapping) + 1)
    print '\t'.join([
        clust_id,
        str(num_tot),
        numbest,
        ','.join(gene_membership[clust_id]),
        outgroup,
        bestgene,
        str(bestdS),
        str(bestdSSE)])

