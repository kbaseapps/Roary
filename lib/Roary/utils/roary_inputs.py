'''
Download a GenomeSet
'''
import os
import sys
import subprocess
import random

from collections import defaultdict

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.GenomeFileUtilClient import GenomeFileUtil
from installed_clients.AssemblyUtilClient import AssemblyUtil


def download_gffs(cb_url, scratch, genome_set_ref):
    """
    Args:
    cb_url - callback server URL
    scratch - scratch work folder
    genome_set_ref - reference to genome_set object in workspace
    Returns the path to the folder containing .gff files


    we want to first handle a GenomeSet Object "KBaseSearch.GenomeSet" or "KBaseSets.GenomeSet"
    """

    # Get our utilities
    dfu = DataFileUtil(cb_url)
    au = AssemblyUtil(cb_url)
    gfu = GenomeFileUtil(cb_url)

    obj_data = dfu.get_objects({'object_refs': [genome_set_ref]})['data'][0]
    gs_obj = obj_data['data']
    obj_type = obj_data['info'][2]

    if 'KBaseSets.GenomeSet' in obj_type:
        refs = [gsi['ref'] for gsi in gs_obj['items']]
    elif 'KBaseSearch.GenomeSet' in obj_type:
        refs = [gse['ref'] for gse in gs_obj['elements'].values()]
    else:
        raise TypeError(
            'provided input must of type KBaseSets.GenomeSet or KBaseSearch.GenomeSet not ' +
            str *
            (obj_type))

    if len(refs) < 2:
        raise ValueError("Must provide GenomeSet with at least 2 Genomes.")

    # name the output directory
    temp_dir = scratch + '/temp'
    final_dir = scratch + '/gff'

    os.mkdir(final_dir)
    os.mkdir(temp_dir)

    # write file that will help us cat the gff and fasta files
    cat_path = scratch + '/fast_cat.txt'

    with open(cat_path, 'w') as cat_file:
        cat_file.write("##FASTA\n")

    path_to_ref_and_ID_pos_dict = {}
    all_ids = set([])

    for ref in refs:
        gen_obj = dfu.get_objects({'object_refs': [ref]})['data'][0]['data']

        # NO Eukaryotes, NO Fungi,
        # yes bacateria, yes archaea, yes(?) virus
        if gen_obj['domain'] not in ['Bacteria', 'Archaea']:
            raise TypeError(
                'Provided Genomes are not labeled as Bacteria or Archaea. Roary is only equipped to handle Archaea or Bacteria')

        fasta_path = temp_dir + "/" + gen_obj['id'] + ".fa"
        gff_file = gfu.genome_to_gff({'genome_ref': ref, 'target_dir': temp_dir})
        if 'assembly_ref' not in gen_obj.keys():
            raise TypeError("All genomes must contain an 'assembly_ref'")
        else:
            fasta_file = au.get_assembly_as_fasta(
                {'ref': gen_obj['assembly_ref'], 'filename': fasta_path})
            # check that fasta_file exists
            if not os.path.isfile(fasta_file['path']):
                raise ValueError('An input Genome does not have an associated FASTA file.')

        # need to figure out if FASTA is already in gff file
        # not sure if we need to do this step.
        if 'path' in gff_file:
            gff_file_path = gff_file['path']
        elif 'file_path' in gff_file:
            gff_file_path = gff_file['file_path']
        elif 'gff_file' in gff_file:
            gff_file_path = gff_file['gff_file']['path']
        else:
            raise ValueError("No GFF File Path found.")

        assert(os.path.isfile(gff_file_path)), "Could not find input GFF file for object with workspace reference: %s"%ref

        # oki doki, here we wanna make sure that the ID's in the genome object match up with
        # ID's in the gff file. This is importatnt because the pangenome object uses the genome
        # objects (in the pangenomeviewer).

        gff_file_path, ID_to_pos, gffid_to_genid, contains_fasta, all_ids = filter_gff(gff_file_path, gen_obj, all_ids)

        new_file_path = final_dir + "/" + gen_obj['id'] + '.gff'

        if contains_fasta:
            args = ['mv', gff_file_path, new_file_path]
            subprocess.call(args)
        else:
            # NOTE: We have to pipe output of cat call to the new_file_path
            # next we make a new 'gff' file that contains both the gff and fasta information

            args = ['cat', gff_file_path, cat_path, fasta_file['path']]
            catted_files = subprocess.check_output(args)
            with open(new_file_path, 'w') as f:
                f.write(catted_files.decode('utf-8'))

        path_to_ref_and_ID_pos_dict[new_file_path] = (ref, ID_to_pos, gffid_to_genid)

    return final_dir, path_to_ref_and_ID_pos_dict


suffix_list = ['_gene']


def filter_gff(gff_file, genome_obj, all_ids =set([]), overwrite=True):
    # if there is a duplicate ID's toss one of them out (if it is programmed frameshift toss one)
    # Here we throw out the second duplicate ID.


    # ideally we would do the feature mapping to the genes, but there aren't always 'gene' features in the gff files
    # so to get around this we base the gene position in the gff off the CDS position because they are most often
    # consecutively ordered and paired. For the organisms that Roary is supposed to service, there should only be
    # a one to one pairing of 'gene' to 'CDS', so this shouldn't be much of an issue

    features = genome_obj['features']

    # we want to make sure that the feature ID's lineup with the gff ID's
    # and we'll change the feature_pos argument to be the position of the ID
    # in the feature array.
    gen_ids = set([])
    ID_to_pos = {}
    for feat_pos, feature in enumerate(features):
        id_ = feature["id"]
        gen_ids.add(id_)
        ID_to_pos[id_] = feat_pos


    with open(gff_file) as f:
        output = []
        gff_ids = set([])
        contains_fasta = False
        for l in f:
            if "##FASTA" in l:
                # add FASTA sequence to output if it exists and then exit
                contains_fasta = True
                output.append(l)
                output += [j for j in f]
                break
            # if we find a comment add it to output
            if l[0:2] == '##':
                output.append(l)
                continue

            # get ID of feat
            ID = l.split('ID=')[-1].split(';')[0]

            # I wonder if there needs to be some more filtering on the ID
            # looks like there may be multiple IDs that are similar to some extent
            # because 

            # determine type of feature
            feat_type = l.split()[2]

            ids_len = len(gff_ids)
            all_ids_len = len(all_ids)
            gff_ids.add(ID)
            all_ids.add(ID)

            if len(gff_ids) == ids_len:
                # we found a dupliacte and its the second one. don't include it
                continue
            if len(all_ids) == all_ids_len:
                # here we want to add a unique identifier to the end of the ID. we use the gff file name.
                l_before, l_after = l.split(ID)[0], l.split(ID)[1]
                ID = ID + '___' + os.path.basename(gff_file).split('.')[0]
                l = l_before + ID + l_after
            if feat_type == 'CDS':
                output.append(l)

    # we want to make sure that the gen_ids     and gff_ids are the same.
    # diff = gen_ids.symmetric_difference(gff_ids)
    gffid_to_genid = map_gff_ids_to_genome_ids(gff_ids, gen_ids, genome_obj)

    assert(len(output) > 1), "Could not succesfully filter %f. It may be empty or contain no CDS information."%gff_file.split('/')[-1]

    if overwrite:
        # (over)write output to file
        with open(gff_file, 'w') as f:
            for l in output:
                f.write(l)

    return gff_file, ID_to_pos, gffid_to_genid, contains_fasta, all_ids


def mapping_func(gff_id, gen_id):
    '''
    function to map gff IDs to genome IDs.
    '''
    if gff_id in gen_id:
        return True
    if gen_id in gff_id:
        return True
    return False


def filter_gff_id(gff_id):
    if gff_id[-4:] == '.CDS':
        gff_id = gff_id[:-4]
    if gff_id[-5:] == '_gene':
        gff_id = gff_id[:-5]
    return gff_id


def map_gff_ids_to_genome_ids(gff_ids, gen_ids, genome_obj):
    '''
    map gff_ids to their corresponding genome object ids for one genome.

    The reason we do this is that the Pangenome viewer requires the ID's in the pangenome
    object to match with the genome object.

    instead of the symmetric difference, we make sure the genome has all the gff ids

    Params:
        gff_ids: list of gff IDs
        gen_ids: list of genome IDs
    Returns:
        gffid_to_genids: map of gff file ID -> genome object ID
    '''

    # start by mapping from genome_id to gff file
    mapping = defaultdict(lambda:[])
    for gff_id in gff_ids:
        gff_id = filter_gff_id(gff_id)
        contained = False
        for gen_id in gen_ids:
            if mapping_func(gff_id, gen_id):
                mapping[gen_id].append(gff_id)
                contained = True
        if not contained:
            raise ValueError('cannot match gff id %s'%gff_id)

    # mapping
    prob_mapping = {key:mapping[key] for key in mapping if len(mapping[key]) > 1}
    used_set = set([mapping[key][0] for key in mapping if len(mapping[key]) == 1])
    prob_dict = dict(prob_mapping)

    # making sure we have 1 to 1 mapping
    iters = 0
    while len(prob_dict) > 1:
        for gen_id, gff_list in prob_mapping.items():
            for index, gff_id in enumerate(gff_list):
                if gff_id in used_set:
                    gff_list.pop(index)
            if len(gff_list) == 1:
                mapping[gen_id] = gff_list
                used_set.add(gff_list[0])
                prob_dict.pop(gen_id, None)                
            else:
                prob_dict[gen_id] = gff_list
        prob_mapping = dict(prob_dict)
        iters+=1
        if iters >20:
            raise ValueError("Could not resolve mapping of \
            KBaseGenomes.Genome object IDs to GFF file IDs.")

    # remap from gen_id to gff_id
    gffid_to_genid = {mapping[key][0]:key for key in mapping}

    if len(gffid_to_genid) != len(gff_id):
        raise ValueError("Genome object with id %s cannot match all of \
                          its IDs to an ID in its GFF File. "%genome_obj['id'])

    return gffid_to_genid



    # --------------------------------------------------------------------------------------------
    # check_id = 'C6Y50_RS11770'
    # if check_id in gff_ids and check_id in gen_ids:
    #     raise ValueError("%s in both genome ids and gff ids in object %s"%(check_id, genome_obj['id']))
    # if check_id in gff_ids:
    #     raise ValueError("%s is in the gff file ID in file %s"%(check_id, genome_obj['id']))
    # if check_id in gen_ids:
    #     raise ValueError("%s is in the genome ID in file %s"%(check_id, genome_obj['id']))
    # --------------------------------------------------------------------------------------------

    # <<<<<<<OLD ONE>>>>>>>>
    '''
    overlap = gen_ids.intersection(gff_ids)
    gffid_to_genid = {}
    if len(overlap) == len(gff_ids):
        # all the ids overlap
        for o in overlap:
            gff_id = filter_gff_id(o)
            gen_id = o
            gffid_to_genid[gff_id] = gen_id
    else:
        # we should make a mapping from the gff ID's to the genome ID's
        for o in overlap:
            gff_id = filter_gff_id(o)
            gen_id = o
            gffid_to_genid[gff_id] = gen_id

        # get non overlapping ones
        diff = gen_ids - gff_ids
        gff_diff = gff_ids - gen_ids

        mapping = defaultdict(lambda:list)
        used_gff = set([])
        for gen_id in diff:
            # find gff_id that matches with the associated genome id
            for gff_id in gff_diff:
                gff_id = filter_gff_id(gff_id)
                if mapping_func(gff_id, gen_id):
                    # this is where they overlap
                    used_gff.add(gff_id)
                    mapping[gen_id].append(gff_id)

        # we have to make sure each gff_id has a mapping to a genome_id
        if len(gff_diff) == len(used_gff):
            # yay
            pass
        else:
            left_over = gff_diff - used_gff
            for gff_id in left_over:
                gff_id = filter_gff_id(gff_id)
                # now find the gen_id
                in_it = False
                for gen_id in gen_ids:
                    if mapping_func(gff_id, gen_id):
                        mapping[gen_id].append(gff_id)
                        in_it = True
                if not in_it:
                    raise ValueError("GFF ID %s cannot be matched to a genome ID"%gff_id )

        problem_map = {key:mapping[key] for key in mapping if len(mapping[key]) > 1}
        used_set = set([mapping[key][0] for key in mapping if len(mapping[key]) == 1])


        # iteratively find pairings for the key values.
        problem_map_copy = dict(problem_map)
        iters = 0
        while len(problem_map) > 1:
            # 1.) filtering step
            for gen_id, gff_list in problem_map.items():
                for index, item in enumerate(gff_list):
                    if item in used_set:
                        gff_list.pop(index)
                # 2.) add items that have 1 to 1 mapping to used list and remove key, value pair
                if len(gff_list) == 1:
                    mapping[gen_id] = gff_list
                    used_set.add(gff_list[0])
                    problem_map_copy.pop(gen_id, None)
                else:
                    problem_map_copy[gen_id] = gff_list
            problem_map = dict(problem_map_copy)
            iters+=1
            if iters > 20:
                raise ValueError("Could not resolve mapping of KBaseGenomes.Genome object IDs to GFF file IDs.")

        # now we should have a complete 1 to 1 mapping.
        for gen_id in mapping:
            gff_id = mapping[gen_id][0]
            # get rid of any suffixes from the list described above
            if gff_id[-5:] in suffix_list:
                gff_id = gff_id[-5:]
            gffid_to_genid[gff_id] = gen_id



    # Check to see if the gffid_to_genid contains all
    if len(gffid_to_genid) != len(gen_ids) and len(gffid_to_genid) != len(gff_ids):
        raise ValueError("Genome object with id %s cannot match all of its IDs to an ID in its GFF File. "%genome_obj['id'])#,

    return gffid_to_genid
    '''
