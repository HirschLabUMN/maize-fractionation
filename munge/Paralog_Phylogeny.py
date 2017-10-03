#!/usr/bin/env python
"""Script to build alignments and phylogenies with putative duplicated loci
in the maize genome. Borrows heavily from BAD_Mutations. Uses TBLASTX to
identify putative homologues in publicly available databases of sequences from
Phytozome and Ensembl Plants, aligns them with PASTA, and estimates a
phylogeny from the sequences."""

try:
    import argparse
except ImportError:
    print "This program requires the 'argparse' module."
    exit(1)

try:
    from Bio import SeqIO
    from Bio import SeqRecord
    from Bio.Seq import Seq
    from Bio.Blast.Applications import NcbitblastxCommandline
    from Bio.Blast import NCBIXML
except ImportError:
    print "This program requires the Biopython library."
    exit(1)

#   Import the standard library modules here
import tempfile
import subprocess
import os
import time
import re
import shutil
import sys

#   A usage message
def usage():
    print """
Usage: Paralog_Phylogeny.py -b BASE -f FOCAL -d DUP [ -d DUP ...] -o OUT

will generate a multiple sequence alignment and maximum likelihood phylogenetic
tree based on best TBLASTX matches to the gene sequence given in the FOCAL
FASTA file. This script searches the BASE directory for BLAST databases, then
uses TBLASTX to search the FOCAL gene against each species. The best matches are
then translated into amino acid sequences and aligned. The sequences given
in the DUP FASTA files are appended to the list of orthologous sequences, but
not used in searching. Phylogenetic relationships are then estimated, and the
alignment is back-translated into nucleotides. The tree and mulitple sequence
alignment are copied to the OUT directory.

To set up the BASE directory, use the 'fetch' subcommand from the
BAD_Mutations program, available at https://github.com/MorrellLAB/BAD_Mutations.
There are instructions for the 'fetch' subcommand in the BAD_Mutations user
manual.

Requires the following packages:
    PASTA (for alignment)
    NCBI BLAST+ executables
    Biopython
    Python 2.7.x
"""
    exit(1)


#   Parse the arguments
def parse_args():
    """Create an ArgumentParser object to validate and store user input."""
    parser = argparse.ArgumentParser(
        description='Align putative paralogs to phylogentic samples.',
        add_help=True)
    parser.add_argument(
        '--base',
        '-b',
        required=False,
        help='Base directory for species databases. Defaults to .',
        default=os.getcwd())
    parser.add_argument(
        '--focal',
        '-f',
        required=True,
        help='Focal gene FASTA file. This is the one that will be used to BLAST'
             ' against the other species.')
    parser.add_argument(
        '--duplicated',
        '-d',
        action='append',
        default=None,
        help='FASTA treated as duplicate. Will be concatenated to alignment '
             'after focal gene BLAST. May be specified multiple times.')
    parser.add_argument(
        '--output',
        '-o',
        default=os.getcwd(),
        help='Output directory for MSA and tree. Defaults to .')
    args = parser.parse_args()
    return args


#   Validate the arguments
def validate_args(a):
    """Validate arugments. Die if any of the files are unreadable."""
    #   See if the base directory exists.
    if not os.path.isdir(a['base']):
        print 'Supplied base directory is not readable, or does not exist.'
        exit(1)
    #   and the destination directory
    if not os.path.isdir(a['base']):
        print 'Supplied output directory is not readable, or does not exist.'
        exit(1)
    #   Then, try the focal gene
    try:
        open(a['focal'])
    except IOError:
        print 'Supplied focal gene file is not readable, or does not exist.'
        exit(1)
    #   Then, try the duplicates
    for d in a['duplicated']:
        try:
            open(d)
        except IOError:
            print ('One of the duplicated gene files is not readable, or does '
                  'not exist.')
            exit(1)
    return a


#   Get a list of the BLAST databases in the base directory
def blast_dbs(base):
    """Finds all files in the base directory that correspond to BLAST
    databases. Lifted from BAD_Mutations."""
    #   First, find all files ending in '.fa' - these are BLAST databases
    command = ' '.join(['find', base, '-name', '"*.fa"'])
    raw_files = subprocess.check_output(command, shell=True)
    file_list = raw_files.strip().split('\n')
    return file_list


#   Run BLAST with the focal gene and return the hits
def run_blast(db, focal):
    """For each database, run TBLASTX against it and pull the best hit."""
    #   Create a temp file
    temp_output = tempfile.NamedTemporaryFile(
        mode='w+t',
        prefix='BLASTSearch_',
        suffix='_BLASTout.xml')
    #   Start building a command line
    cline = NcbitblastxCommandline(
        query=focal,
        out=temp_output.name,
        db=db,
        evalue=0.05,
        outfmt=5,
        max_target_seqs=1)
    #   And then execute it
    cline()
    return temp_output


#   Extract the best hit from the results handle
def get_hit(results):
    """Parse the XML BlAST report, extract the sequence ID of the best hit."""
    #   Define a helper function for extracting the best hit out of the ugly
    #   BlastRecord object
    def best_hit(brecord):
        for aln in brecord.alignments:
            for hsp in aln.hsps:
                frames = [str(f) for f in hsp.frame]
                if hsp.expect <= 0.05:
                    return aln.title
        else:
            return None

    #   Seek to the beginning of the handle for reading.
    results.seek(0)
    #   Parse the XML
    blast_records = list(NCBIXML.parse(results))
    #   If the list is empty, then we have no hits
    if len(blast_records) > 0:
        #   For each record
        for rec in blast_records:
            best = best_hit(rec)
            if best:
                break
    else:
        return None
    #   Close the temporary file to clean up; it's automatically deleted
    results.close()
    return best


#   Define a function to get the sequence of the homologue
def extract_hit_seq(db, hit_id):
    """Uses the 'blastdbcmd' command to extract the FASTA sequence of the best
    BLAST hit."""
    #   Run blastdbcmd to get the sequence ID from the database
    cmd = ['blastdbcmd', '-entry', hit_id, '-db', db]
    p = subprocess.Popen(
        cmd,
        shell=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    out, err = p.communicate()
    return (out, err)


#   Translate the sequences into amino acid sequences
def prepare_sequences(nuc):
    """Prepares the CDS sequences for alignment with Pasta. Checks if any
    sequences are not multiples of 3, and appends N if not. Translates
    the nucleotide sequences to amino acids for protein alignment with
    Pasta, then removes the trailing stop codon, if it is present."""
    #   First, read in the sequences and check their length
    input_seqs = list(SeqIO.parse(nuc.name, 'fasta'))
    #   Start an empty dictionary for holding the nucleotide sequences
    orig_seqs = {}
    #   Start accumulating translated sequences to write into the
    #   alignment input file.
    tl_seqs = []
    for i in input_seqs:
        #   If the length is not a multiple of 3, then we have to add Ns to
        #   make it so.
        if len(i) % 3 != 0:
            to_add = 3 - (len(i) % 3)
            #   Tack on the original sequence, with some appended Ns so we
            #   can recreate the nucleotide alignment later.
            orig_seqs[i.id] = str(i.seq) + to_add*'N'
            #   Create the new sequence. Pasta chokes on ambiguous
            #   amino acids that aren't X, so we replace them all with X.
            fixed_seq = Seq(str(i.seq) + to_add*'N').translate()
            sub_seq = Seq(
                re.sub(
                    'B|Z|J|O|U',
                    'X',
                    str(fixed_seq),
                    re.I))
            new_seq = SeqRecord.SeqRecord(
                sub_seq,
                id=i.id,
                description='')
        else:
            orig_seqs[i.id] = str(i.seq)
            sub_seq = Seq(
                re.sub(
                    'B|Z|J|O|U',
                    'X',
                    str(i.seq.translate()),
                    re.I))
            new_seq = SeqRecord.SeqRecord(
                sub_seq,
                id=i.id.replace(':', '_'),
                description='')
        tl_seqs.append(new_seq)
    #   Then, we have to iterate through the translated sequences and
    #   check for sequences ending in stop codons. Pasta hates these, so
    #   we will prune them. We also check for those with internal stop
    #   codons, and skip those.
    fixed_tl_seqs = []
    for i in tl_seqs:
        #   If we find an internal stop codon
        if re.match(r'.+\*[^$]', str(i.seq)):
            continue
        else:
            #   Otherwise, just strip the stop codon off the end and save
            if i.seq.endswith('*'):
                fixed_tl_seqs.append(
                    SeqRecord.SeqRecord(i.seq[:-1], id=i.id, name='', description='')
                    )
            else:
                fixed_tl_seqs.append(i)
    #   Then, we open another temporary file to hold our amino acid
    #   sequences.
    prot_in = tempfile.NamedTemporaryFile(
        mode='w+t',
        prefix='PastaInput_',
        suffix='.fasta')
    #   And write the protein sequences into it
    SeqIO.write(fixed_tl_seqs, prot_in, 'fasta')
    prot_in.flush()
    return (orig_seqs, prot_in)


#   Then align the sequences, return the names of the alignment and tree files
def pasta_align(prot):
    """Align the amino acid sequences with Pasta."""
    #   Get the temp out directory, because Pasta writes a lot of stuff there.
    pasta_out = tempfile.gettempdir()
    #   We make a job name from the time in microseconds
    #   This should be good enough...
    pasta_job = 'pastajob_' + '%.6f' % time.time()
    #   Create the command line
    cmd = [
        'run_pasta.py',
        '-d',
        'protein',
        '--no-return-final-tree-and-alignment',
        '--num-cpus=1',
        '--job='+pasta_job,
        '--temporaries=' + pasta_out,
        '-i',
        prot.name,
        '--mask-gappy-sites=12',
        '--iter-limit=5',
        '-o',
        pasta_out]
    #   Then, we'll execute it
    p = subprocess.Popen(
        cmd,
        shell=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE)
    out, err = p.communicate()
    #   Then, build the output name
    #   The structure of the Pasta output files is
    #       Jobname.marker001.Unaligned_name.aln
    #       Jobname.tre
    aln_out = ''.join([
        pasta_out,
        os.path.sep,
        pasta_job,
        '.marker001.',
        prot.name.split(os.path.sep)[-1].replace('.fasta', ''),
        '.aln'])
    tree_out = ''.join([
        pasta_out,
        os.path.sep,
        pasta_job,
        '.tre'])
    return (aln_out, tree_out)


#   Back-translate the alignment into nucleotides
def back_translate(p_aln, orig):
    """Back-translates from amino acid to nucleotide, using the original
    input sequences as a guide to avoid ambiguity. Assumes that a non-gap
    character in the amino acid alignment will be faithfully represented
    by a triplet in the source sequence, and will not check identity of
    translated codons."""
    #   Iterate through the aligned protein file, and start rebuilding the
    #   original nucleotide sequence.
    aln_prot = SeqIO.parse(p_aln, 'fasta')
    bt_seqs = []
    for rec in aln_prot:
        #   Check if we have name mismatch. This *shouldn't* happen, but
        #   this should stop some errors
        if rec.id in orig:
            in_seq = orig[rec.id]
            rebuilt_seq = ''
            pos = 0
            #   Now, iterate through the amino acid sequence and rebuild
            #   the nucleotide sequence.
            for aa in rec.seq:
                #   If we have a gap, put in three gaps
                if aa == '-':
                    rebuilt_seq += '---'
                else:
                    rebuilt_seq += in_seq[pos:pos+3]
                    pos += 3
        #   Then, put in a new SeqRecord with the sequence and the name, so
        #   we can write a fasta file.
        bt_seqs.append(SeqRecord.SeqRecord(
            Seq(rebuilt_seq),
            id=rec.id,
            name='',
            description=''))
    #   And create a new temporary file for them
    final_seqs = tempfile.NamedTemporaryFile(
        mode='w+t',
        prefix='BackTranslated_',
        suffix='.fasta',
        delete=False)
    SeqIO.write(bt_seqs, final_seqs, 'fasta')
    final_seqs.flush()
    return final_seqs


def main():
    """Main function."""
    #   Check the number of arguments first
    if len(sys.argv) == 1:
        usage()
    #   vars() returns a dictionary. For some reason, argparse doesn't want to
    #   return just a dictionary.
    args = vars(parse_args())
    v_args = validate_args(args)
    bdbs = blast_dbs(args['base'])
    #   Create a named temporary file for storing the homologue sequences
    homologues = tempfile.NamedTemporaryFile(
        mode='w+t',
        prefix='BestHits_',
        suffix='_SeqList.fasta')
    for db in bdbs:
        b_results = run_blast(db, args['focal'])
        hit = get_hit(b_results)
        if hit:
            o, e = extract_hit_seq(db, hit.split(' ')[0])
            db_seq = ''.join(o.split('\n')[1:])
            #   Then, search for the ambiguous nucleotides, skip if they are
            #   found.
            if re.search('S|W|R|K|Y|M|V|B|D|H|N', db_seq, re.I):
                continue
            else:
                spname = os.path.basename(db)
                #   Then split on . and take the first part
                spname = '>' + spname.split('.')[0]
                #   Remove any weird characters in the species name
                weird = ['/', '\\', ':', '+', '.']
                for c in weird:
                    spname = spname.replace(c, '')
                fasta = spname + '\n' + db_seq + '\n'
                homologues.write(fasta)
                homologues.flush()

    #   Next, tack on the query sequence, and the duplicate sequences
    for query in [v_args['focal']] + v_args['duplicated']:
        s = SeqIO.read(query, 'fasta')
        q_fasta = '>' + s.id + '\n' + s.seq + '\n'
        homologues.write(str(q_fasta))
        homologues.flush()

    #   Then, prepare the sequences for PASTA alignment
    orig_seqs, prot_seqs = prepare_sequences(homologues)
    #   And align them
    alignment, tree = pasta_align(prot_seqs)
    #   Backtranslate
    nuc_aln = back_translate(alignment, orig_seqs)
    #   Then, copy the files to the destination directory
    paln_name = os.path.basename(args['focal']).replace('.fasta', '_prot_MSA.fasta')
    naln_name = os.path.basename(args['focal']).replace('.fasta', '_nuc_MSA.fasta')
    tree_name = os.path.basename(args['focal']).replace('.fasta', '.tree')
    full_paln_name = os.path.join(args['output'], paln_name)
    full_naln_name = os.path.join(args['output'], naln_name)
    full_tree_name = os.path.join(args['output'], tree_name)
    open(full_paln_name, 'w').close()
    open(full_naln_name, 'w').close()
    open(full_tree_name, 'w').close()
    shutil.copy2(alignment, full_paln_name)
    shutil.copy2(nuc_aln.name, full_naln_name)
    shutil.copy2(tree, full_tree_name)

    #   Clean up
    prot_seqs.close()
    nuc_aln.close()

main()
