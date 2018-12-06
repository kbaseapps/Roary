import sys
import os
import json
import subprocess


def run_roary(scratch, gff_folder, params):
    """
    params :
    --------
            scratch    : scratch folder
            gff_folder : file path to folder the input .gff files are located in
            params     : set of parameters for running roary
            -------
                    blast_p_percentage     : minimum percent identity threshold for BLAST P, default 95
                    max_num_clusters       : maximum allowable number of cluster to create, default 50000
                    percent_genes_for_core : percentage of isolates a gene must be in to be core, default 99
    """
    # first verify that there are multiple (at least 2) .gff files in the gff_folder
    gff_folder_dir = os.listdir(gff_folder)
    gff_files = []
    for f in gff_folder_dir:
        if '.gff' in f:
            gff_files.append(os.path.join(gff_folder, f))
    if len(gff_files) < 2:
        raise ValueError("Must Provide at least two Genomes with GFF files.")

    out_dir = os.path.join(scratch, 'outputs')

    # arguments we probably want to include in UI
    blastp_percentage_identity = str(params['blast_p_percentage'])
    no_split_paralogs = False
    max_num_clusters = str(params['max_num_clusters'])
    percent_genes_for_core = str(params['percent_genes_for_core'])

    # internal arguments
    graphing = False
    num_threads = str(4)

    # find out what other arguments we want to use for this from Paramvir
    args = ['roary', '-f', out_dir, '-p', num_threads, '-i',
            blastp_percentage_identity, '-g', max_num_clusters, '-cd', percent_genes_for_core]
    if no_split_paralogs:
        args.append('-s')

    if graphing:
        args.append('-r')
    # finally add files on which to run Roary
    args += gff_files

    # Roary seems to return code 255 with subprocess.check_output,
    proc = subprocess.Popen(args)
    res = proc.wait()

    if res != 0:
        raise RuntimeError("Roary subprocess exited with error code: %f" % res)

    return out_dir
