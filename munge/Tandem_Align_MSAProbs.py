#!/usr/bin/env python
"""Take aligned tandem cluster sequences and re-align them with MSAprobs. This
will only run on Linux (at least, unless a MacOS MSAProbs binary is available).
Takes two arguments:
    1) Directory ofsequence alignments
    2) Output directory
"""

import sys
import os
import subprocess
import tempfile
try:
    from Bio import SeqIO
    from Bio import Seq
    from Bio import SeqRecord
except ImportError:
    print 'This script requires Biopython.'
    exit(1)


def get_alignments(aln_dir, ext='.fasta'):
    """Given a directory, find all files with the given extension, and return
    as a list."""
    full_dir = os.path.abspath(os.path.expanduser(aln_dir))
    dir_contents = os.listdir(full_dir)
    alns = [f for f in dir_contents if f.endswith(ext)]
    # Generate full paths for the alignments
    full_alns = [os.path.join(full_dir, f) for f in alns]
    return full_alns


def translate_dealign(aln):
    """De-align the aligned sequences (remove gaps), and translate them to
    amino acid. Write the de-aligned amino acid sequences to a temporary file,
    to align later."""
    to_realign = []
    for s in SeqIO.parse(open(aln), 'fasta'):
        nogaps = ''.join([b for b in s.seq if b != '-'])
        # Translate it
        nogaps_trans = str(Seq.Seq(nogaps).translate())
        # Remove any stop codons
        nogaps_trans = nogaps_trans.split('*')[0]
        to_realign.append(
            SeqRecord.SeqRecord(
                Seq.Seq(nogaps_trans),
                id=s.id,
                name='',
                description=''
                )
            )
    # Open a temporary file
    t = tempfile.NamedTemporaryFile(
        mode='w+t',
        prefix='msaprobs_realign',
        suffix='.fasta'
        )
    SeqIO.write(to_realign, t.name, 'fasta')
    return t


def run_msaprobs(infile, a_name, outdir):
    """Run MSAprobs to align the tandem sequences. Since MSAprobs is not
    installed on MSI, and it does not compile nicely on MacOS, you must specify
    the path to the msaprobs binary."""
    # Seek to the beginning of the input file so that it can be read
    infile.seek(0)
    # Set full path to msaprobs binary
    msaprobs_path = '/media/sf_VBox_Share/MSAProbs-0.9.7/MSAProbs/msaprobs'
    # Full path to the outputt directory
    out_path = os.path.abspath(os.path.expanduser(outdir))
    # basename for the output file
    out_fname = os.path.basename(a_name).replace('.fasta', '_msaprobs.fasta')
    cmd = [
        msaprobs_path,
        '-o',
        os.path.join(out_path, out_fname),
        '-ir',
        '20',
        infile.name
    ]
    # Run the command
    p = subprocess.Popen(
        cmd,
        shell=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
        )
    out, err = p.communicate()
    # Close the tempfile
    infile.close()
    # return the output name
    return os.path.join(out_path, out_fname)


def backtranslate(amino_acid, nucleotide):
    """Use the nucleotide sequences to convert the aligned amino acids back to
    nucleotides."""
    # Store the non-gapped nucleotide sequences as a dictionary
    nuc_seqs = []
    for s in SeqIO.parse(nucleotide, 'fasta'):
        snogap = ''.join([b for b in s.seq if b != '-'])
        nuc_seqs.append(
            SeqRecord.SeqRecord(
                seq=Seq.Seq(snogap),
                id=s.id,
                name='',
                description=''
                )
            )
    nuc_seqs = SeqIO.to_dict(nuc_seqs)
    # Then, back-translate the aligned sequences
    bt_seq = []
    for prot in SeqIO.parse(amino_acid, 'fasta'):
        codon = 0
        bt = ''
        nuc = nuc_seqs[prot.id]
        for aa in prot:
            if aa == '-':
                bt += '---'
            else:
                bt += nuc[codon*3:(codon*3)+3]
                codon += 1
        bt_seq.append(bt)
    # Write the backtranslated sequences to disk
    bt_name = amino_acid.replace('msaprobs.fasta', 'backtranslated.fasta')
    SeqIO.write(bt_seq, bt_name, 'fasta')
    return


def main(aln_dir, out_dir):
    """Main function."""
    alns = get_alignments(aln_dir)
    for a in alns:
        aa = translate_dealign(a)
        aligned_aa = run_msaprobs(aa, a, out_dir)
        backtranslate(aligned_aa, a)
    return

main(sys.argv[1], sys.argv[2])
