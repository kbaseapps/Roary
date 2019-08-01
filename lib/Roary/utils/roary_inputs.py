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


def download_gffs(cb_url, scratch, input_refs):
    """
    Args:
    cb_url - callback server URL
    scratch - scratch work folder
    input_refs - list of references to genome_set or genome objects in workspace
    Returns the path to the folder containing .gff files


    we want to first handle a GenomeSet Object "KBaseSearch.GenomeSet" or "KBaseSets.GenomeSet"
    """

    # Get our utilities
    dfu = DataFileUtil(cb_url)
    au = AssemblyUtil(cb_url)
    gfu = GenomeFileUtil(cb_url)

    obj_data = dfu.get_objects({'object_refs': input_refs})['data']

    refs = []
    for datum in obj_data:
        gs_obj = datum['data']
        obj_type = datum['info'][2]

        if 'KBaseSets.GenomeSet' in obj_type:
            curr_refs = [gsi['ref'] for gsi in gs_obj['items']]
        elif 'KBaseSearch.GenomeSet' in obj_type:
            curr_refs = [gse['ref'] for gse in gs_obj['elements'].values()]
        elif 'KBaseGenome.Genome' in obj_type:
            # not most efficient answer here, but will work for now.
            curr_refs = ['/'.join([datum['info'][6], datum['info'][0], datum['info'][4]])]
        else:
            raise TypeError(
                'provided input(s) must of type KBaseGenomes.Genome, KBaseSets.GenomeSet or '
                ' KBaseSearch.GenomeSet not ' + str(obj_type))
        refs += curr_refs

    if len(refs) < 2:
        raise ValueError("Must provide GenomeSet with at least 2 Genomes, or multiple Genomes")

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
        # NOTE: we are getting rid of this because it is not consistent enough for now...

        # if gen_obj['domain'] not in ['Bacteria', 'Archaea']:
        #     raise TypeError('Provided Genomes are not labeled as Bacteria or Archaea. '
        #                     'Roary is only equipped to handle Archaea or Bacteria')

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

        gen_id_to_pos, contains_fasta, all_ids, gen_ids = filter_gff(gff_file_path, gen_obj, all_ids=all_ids)

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

        path_to_ref_and_ID_pos_dict[new_file_path] = (ref, gen_id_to_pos, gen_ids)

    return final_dir, path_to_ref_and_ID_pos_dict

def toString(s):
    try:
        return str(s, 'utf-8')
    except:
        return s

def filter_gff(gff_file, genome_obj, all_ids=set([]), overwrite=True):
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
    # gen_ids = set([])
    # mrna_id_to_id = {}
    # gff_id_and_type = {}
    gen_ids = []
    gen_id_to_pos = {}
    for feat_pos, feature in enumerate(features):
        if feature.get('id'):
            id_ = feature.get('id')
            id_ = toString(id_)
        else:
            id_ = ""

        if feature.get('cdss'):
            cdss = feature.get('cdss')
            cds = toString(cdss[0])
        else:
            cds = ""

        if feature.get('mrnas'):
            mrnas = feature.get('mrnas')
            mrna = toString(mrnas[0])
        else:
            mrna = ""
        gen_id_to_pos[id_] = feat_pos
        gen_ids.append((id_, cds, mrna))

    gff_ids = set([])
    parent_ids = set([])
    gff_output = []
    length = 0
    with open(gff_file) as f:
        contains_fasta = False
        for l in f:
            if '##FASTA' in l:
                contains_fasta = True
                gff_output.append(l)
                gff_output += [j for j in f]
                break
            if l[:2] =='##':
                gff_output.append(l)
                continue
            feat_type = l.split()[2]
            ID = l.split('ID=')[-1].split(';')[0]
            if 'Parent=' in l:
                parent = l.split("Parent=")[-1].split(";")[0]
            else:
                parent = ''
            if feat_type == "CDS":

                pre_len = len(gff_ids)
                all_pre_len = len(all_ids)

                gff_ids.add(ID)
                all_ids.add(ID)
                if parent:
                    parent_ids.add(parent)
                else:
                    parent_ids.add(ID)

                if len(gff_ids) == pre_len:
                    continue
                elif len(all_ids) == all_pre_len:
                    l_before, l_after = l.split(ID)[0], l.split(ID)[1]
                    ID = ID + "___" + os.path.basename(gff_file).split('.')[0]
                    l = l_before + ID + l_after

                gff_output.append(l)
            length+=1

    # if len(gff_ids) > len(gen_ids):
    #     raise ValueError("More gff ids than there are available genome ids")

    # lets try mapping at end instead
    # gff_id_to_gen_id, gen_ids = map_gff_to_gen(gen_ids, list(parent_ids), gff_file)

    if overwrite:
        with open(gff_file, 'w') as f:
            for l in gff_output:
                f.write(l.rstrip() + '\n')

    return gen_id_to_pos, contains_fasta, all_ids, gen_ids

# def find_pair(gff_id, gen_ids):
#     gff_id = toString(gff_id)
#     for i, gds in enumerate(gen_ids):
#         id_, cds, mrna = gds
#         id_, cds, mrna = toString(id_), toString(cds), toString(mrna)
#         if gff_id == id_:
#             return i
#         if gff_id == cds:
#             return i
#         if gff_id == mrna:
#             return i
#         # last resort check if the gff_id is a substring:
#         if in id_:
#             return i
#     return None


# def map_gff_to_gen(gen_ids, gff_ids, gff_file):
#     gff_id_to_gen_id = {}
#     print('gff_file', gff_file)
#     print('length of gen ids:', len(gen_ids))
#     print('length of gff ids:', len(gff_ids))
#     for gff_id in gff_ids:
#         gff_id = toString(gff_id)
#         i = find_pair(gff_id, gen_ids)
#         if i == None:
#             # if we can't find a pairing, map them together
#             gff_id_to_gen_id[gff_id] = toString(gff_id)
#             # raise ValueError("gff id %s has no matching genome object id"%(gff_id),  gen_ids[:30])
#         else:
#             val = gen_ids.pop(i)
#             gff_id_to_gen_id[gff_id] = toString(val[0])
#     print('gff id to gen id:', len(gff_id_to_gen_id))
#     return gff_id_to_gen_id, gen_ids
