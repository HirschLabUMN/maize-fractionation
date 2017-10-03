#!/usr/bin/env python
"""Calculate dNdS for the alignments of tandem duplicates against Sorghum
bicolor, using the 'yn00' routine of PAML. It uses the Yang and Nielsen (2000)
method for estimating dN/dS, which uses ti:tv ratios and codon usage. Takes two
arguments: 
   1) Path to alignments of tandem duplicates
   2) Output directory"""

import sys
import os
import shutil
import subprocess


def list_alignments(aln_dir, ext='.fasta'):
    """List the full path to files in the alignment directory, so that it is
    easy to iterate through them."""
    dir_contents = os.listdir(aln_dir)
    # Filter to just files with the desired extension
    alignments = [os.path.basename(f) for f in dir_contents if f.endswith(ext)]
    # Append the dirname to the filenames. The os.path.expanduser() function
    # correctly handles the ~ in paths.
    full_aln_dir = os.path.abspath(os.path.expanduser(aln_dir))
    full_path_alignments = [os.path.join(full_aln_dir, f) for f in alignments]
    return full_path_alignments


def generate_ctl(a_name, o_dir):
    """Generate a control file for an analysis, given a sequence name. The
    sample control file from the PAML documentation is:

    seqfile = abglobin.nuc * sequence data file name
    outfile = yn           * main result file
    verbose = 0      * 1: detailed output (list sequences), 0: concise output
    icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below
    weighting = 0  * weighting pathways between codons (0/1)?
    commonf3x4 = 0  * use one set of codon freqs for all pairs (0/1)?"""
    # Open a handle to the control file
    handle = open(os.path.join(o_dir, 'yn00.ctl'), 'w')
    handle.write('seqfile = ' + a_name + '\n')
    handle.write('outfile = ' + a_name + '.yn00_out\n')
    handle.write('verbose = 1\n')
    handle.write('icode = 0\n')
    handle.write('weighting = 0\n')
    handle.write('commonf3x4 = 0\n')
    handle.flush()
    handle.close()
    return


def setup_analysis(aln, o_dir):
    """Given an alignment file and an output directory to run the analysis,
    set up the analysis environment. All the materials for the analysis of one
    gene should be present in its own directory, with a .ctl file."""
    # Get the basename of the alignment, which would be the cluster ID
    aln_basename = os.path.basename(aln)
    # Strip off the extension
    clust_id = aln_basename.split('.')[0]
    # Then make the directory for the analysis
    analysis_dir = os.path.join(
        os.path.abspath(os.path.expanduser(o_dir)),
        clust_id)
    if not os.path.exists(analysis_dir):
        os.makedirs(analysis_dir)
    # Generate the ctl file in the analysis directory
    generate_ctl(aln_basename, analysis_dir)
    # Copy the alignment into the analysis directory
    shutil.copy2(aln, analysis_dir)
    return


def run_analysis(aln, o_dir):
    """Run the analysis."""
    # Get the full path of the analysis directory and cd into it
    aln_id = os.path.basename(aln).split('.')[0]
    analysis_dir = os.path.join(
        os.path.abspath(os.path.expanduser(o_dir)),
        aln_id)
    os.chdir(analysis_dir)
    # Run the anaylsis with the subprocess module, saving stdout and stderr
    cmd = ['yn00', 'yn00.ctl']
    p = subprocess.Popen(
        cmd,
        shell=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
        )
    out, err = p.communicate()
    # Save the output from stdout and stderr
    handle = open(aln_id + '.yn00_stdout', 'w')
    handle.write(out)
    handle.flush()
    handle.close()
    handle = open(aln_id + '.yn00_stderr', 'w')
    handle.write(err)
    handle.flush()
    handle.close()
    return


def main(tandem_dir, dest_dir):
    """Main function."""
    cluster_alignments = list_alignments(tandem_dir)
    # get the full path to the destination directory
    full_dest = os.path.abspath(os.path.expanduser(dest_dir))
    for clust in cluster_alignments:
        setup_analysis(clust, full_dest)
        run_analysis(clust, full_dest)
    return


# execute the main function
main(sys.argv[1], sys.argv[2])
