#!/usr/bin/env python
"""Reads the list of tandem duplicates and tries to determine if one of them
is likely a pseudogene based on sequence similarity to known TEs and the number
of exons. Takes three arguments:
    1) Tandem duplicates file, comma-seprated, one cluster per line
    2) GFF of the representative transcripts for each gene
    3) FASTA with the longest CDS of each gene
    4) Path to TE BLAST database"""

import sys
import subprocess
import tempfile

try:
    from Bio.Blast.Applications import NcbiblastnCommandline
    from Bio.Blast import NCBIXML
except ImportError:
    print "You need Biopython to run this script."
    exit(1)


def parse_gff(g):
    """Read GFF data. For each gene, store its length and the number of exons
    in its representative transcript. This will return a data structure of the
    following form:
    {
        gene_id: {
            Length: INT,
            Chromosome: STR,
            Start: INT,
            Stop: INT,
            Exons: [EX1, EX2, ..., ]
            },
        gene_id: {
            Length: INT,
            Chromosome: STR,
            Start: INT,
            Stop: INT,
            Exons: [EX1, EX2, ..., ]
            },
        ...
    }"""
    gff_dat = {}
    with open(g, 'r') as f:
        for line in f:
            #   Skip header lines
            if line.startswith('#'):
                continue
            else:
                tmp = line.strip().split('\t')
                #   We want to save information from genes
                if tmp[2] == 'gene':
                    metadata = tmp[-1].split(';')
                    for m in metadata:
                        if m.startswith('ID'):
                            gid = m[3:]
                            break
                    if gid not in gff_dat:
                        gff_dat[gid] = {
                            'Length': int(tmp[4]) - int(tmp[3]),
                            'Chromosome': tmp[0],
                            'Exons': []}
                elif tmp[2] == 'exon':
                    metadata = tmp[-1].split(';')
                    for m in metadata:
                        if m.startswith('Parent'):
                            parent = m[7:].split('_')[0]
                        if m.startswith('ID') or m.startswith('exon_id'):
                            eid = m[3:]
                    if parent in gff_dat:
                        gff_dat[parent]['Exons'].append(eid)
                elif tmp[2] == 'CDS':
                    metadata = tmp[-1].split(';')
                    for m in metadata:
                        if m.startswith('Parent'):
                            parent = m[7:].split('_')[0]
    return gff_dat


def te_check(geneid, cds, te_db, evalue='1e-10'):
    """BLAST the gene ID against the TE database, and determine if any of the
    sequences produce significant hits. Return e-values of the hits if it does,
    return False othersise."""
    #   Open a temporary file to hold the genomic sequence of the gene
    dup_temp = tempfile.NamedTemporaryFile(
        mode='w+t',
        prefix='Fractionation_BlastSearch_',
        suffix='.fasta')
    blast_temp = tempfile.NamedTemporaryFile(
        mode='w+t',
        prefix='Fractionation_BlastSearch_',
        suffix='.xml')
    #   Build the command line for sequence extraction
    seq_extract = [
        'samtools',
        'faidx',
        cds,
        geneid]
    #   Run the command, save the stdout for writing into the tempfile
    p = subprocess.Popen(
        seq_extract,
        shell=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    out, err = p.communicate()
    dup_temp.write(out)
    #   Then run BLAST
    cline = NcbiblastnCommandline(
        query=dup_temp.name,
        out=blast_temp.name,
        db=te_db,
        evalue=evalue,
        perc_identity=75,
        qcov_hsp_perc=50,
        outfmt=5)
    cline()
    #   Then parse the BLAST output
    blast_temp.seek(0)
    try:
        blast_records = NCBIXML.parse(blast_temp)
        blast_records = list(blast_records)
    except ValueError:
        return False
    evals = []
    names = []
    #   Else, save the e-values, and the sequences that match
    for rec in blast_records:
        if len(rec.alignments) == 0:
            return False
        for aln in rec.alignments:
            for hsp in aln.hsps:
                evals.append(hsp.expect)
                names.append(aln.title)
    #   Close the temp file handles to clean up
    dup_temp.close()
    blast_temp.close()
    return (evals, names)


def main(tandem, gff_file, cds, te_db):
    """Main function."""
    #   Print a header
    print '\t'.join(
        [
            'GenesInCluster',
            'Lengths',
            'NumExons',
            'OneExon',
            'WithTE',
            'OneExonTE',
            'NumExonTE',
            'NumExonWOTE',
            'LengthsWithTE',
            'LengthsWOTE',
            'TEsRepresented',
            'Filtered'
        ])
    gff_data = parse_gff(gff_file)
    #   Start reading through the tandem duplicates file
    with open(tandem, 'r') as f:
        for line in f:
            tmp = sorted(line.strip().split(','))
            num_exons = [len(gff_data[g]['Exons']) for g in tmp]
            lengths = ','.join([str(gff_data[g]['Length']) for g in tmp])
            #   Get the genes that have only one exon
            one_exon = [x[0] for x in zip(tmp, num_exons) if x[1] == 1]
            #   Get which genes have TE sequence similarity
            te_genes = []
            te_names = []
            for geneid in tmp:
                te_sim = te_check(geneid, cds, te_db)
                if te_sim:
                    te_genes.append(geneid)
                    n = [s.split(' ')[1] for s in te_sim[1]]
                    te_names += n
            #   Get only the unique TE names
            te_names = list(set(te_names))
            #   Now, print out the table!
            l_withte = [str(gff_data[g]['Length']) for g in te_genes]
            l_wote = [str(gff_data[g]['Length']) for g in tmp if g not in te_genes]
            if len(te_names) == 0:
                te_pres = 'NA'
            else:
                te_pres = ','.join(te_names)
            if len(l_withte) == 0:
                lenwithte = 'NA'
            else:
                lenwithte = ','.join(l_withte)
            if len(l_wote) == 0:
                lenwote = 'NA'
            else:
                lenwote = ','.join(l_wote)
            cluster = ','.join(tmp)
            if len(one_exon) == 0:
                oexon = 'NA'
            else:
                oexon = ','.join(one_exon)
            if len(te_genes) == 0:
                withte = 'NA'
            else:
                withte = ','.join(te_genes)
            representative = [gid for gid in tmp if gid not in one_exon + te_genes]
            if len(representative) == 0:
                rep = 'NA'
            else:
                rep = ','.join(representative)
            nexons = ','.join([str(s) for s in num_exons])
            exonte = list(set(one_exon).intersection(set(te_genes)))
            if len(exonte) == 0:
                exonte = 'NA'
            else:
                exonte = ','.join(exonte)
            nexonswote = [str(x[1]) for x in zip(tmp, num_exons) if x[0] not in te_genes]
            nexonswte = [str(len(gff_data[g]['Exons'])) for g in te_genes]
            if len(nexonswte) == 0:
                numexon_wte = 'NA'
            else:
                numexon_wte = ','.join(nexonswte)
            if len(nexonswote) == 0:
                numexon_wote = 'NA'
            else:
                numexon_wote = ','.join(nexonswote)
            print '\t'.join(
                [
                    cluster,
                    lengths,
                    nexons,
                    oexon,
                    withte,
                    exonte,
                    numexon_wte,
                    numexon_wote,
                    lenwithte,
                    lenwote,
                    te_pres,
                    rep
                ])
    return


main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
