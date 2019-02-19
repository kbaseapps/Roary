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

    # assert that input paramters are there.
    assert('max_num_clusters' in params), "Maximum number of cluster argument required"
    assert('blast_p_percentage' in params), "minimum BlastP Percentage argument required"
    assert('percent_genes_for_core' in params), "Percent Genes for Core argument required"

    # arguments we probably want to include in UI
    blastp_percentage_identity = str(params['blast_p_percentage'])
    no_split_paralogs = False
    max_num_clusters = str(params['max_num_clusters'])
    percent_genes_for_core = str(params['percent_genes_for_core'])

    # internal arguments
    graphing = True
    num_threads = str(8)

    # run with the following arguments
    args = ['roary', '-f', out_dir, '-p', num_threads, '-i',
            blastp_percentage_identity, '-g', max_num_clusters, '-cd', percent_genes_for_core]

    if no_split_paralogs:
        args.append('-s')

    if graphing:
        args.append('-r')
    # finally add files to run Roary on
    args += gff_files

    print('Running Roary')
    proc = subprocess.Popen(args)
    res = proc.wait()
    print('Roary run complete')
    if res != 0:
        error_message = "Roary subprocess exited with error code: %i" % res
        raise RuntimeError(error_message)

    return out_dir
