import os
import uuid
import pandas as pd
import numpy as np
import json
from shutil import copyfile

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.WSLargeDataIOClient import WsLargeDataIO

# utils
from .roary_output import format_output_html

def get_col_name_from_path(path):
    return os.path.splitext(path.split('/')[-1])[0]

def toString(s):
    try:
        return str(s, 'utf-8')
    except:
        return s

def filter_gff_id(gff_id):
    if gff_id[-4:] == '.CDS':
        gff_id = gff_id[:-4]
    if gff_id[-6:] == "_CDS_1":
        gff_id = gff_id[:-6]
    return gff_id

def gff_id_error(gff_id, gff_id_to_gen_id, col):
    keys = sorted(list(gff_id_to_gen_id.keys()))
    pre_per = gff_id.split('.')[0]
    key_vals_of_substr = []
    for i in range(1, len(gff_id)):
        if gff_id[:i] in keys:
            key_vals_of_substr.append((gff_id[:i], gff_id_to_gen_id[gff_id[:i]]))

    raise KeyError("gff ID %s from file %s, not in dict: "%(gff_id, col), {key:gff_id_to_gen_id[key] for key in keys[:50]})


def find_pair(gff_id, gen_ids):
    gff_id = toString(gff_id)
    for i, gds in enumerate(gen_ids):
        id_, cds, mrna = gds
        id_, cds, mrna = toString(id_), toString(cds), toString(mrna)
        if gff_id == id_:
            return i
        if gff_id == cds:
            return i
        if gff_id == mrna:
            return i
        # last resort check if the gff_id is a substring:
        if gff_id in id_:
            return i
    return None


def generate_pangenome(gene_pres_abs, path_to_ref_and_ID_pos_dict, pangenome_id, pangenome_name):
    '''
        params:
            gene_pres_abs               : file path to gene_presence_absence.csv output from Roary
            path_to_ref_and_ID_pos_dict : dictionary mapping gff file path to a tuple of (workspace ref, 
                                          {gene ID :-> file position}, {genome object id -> gff id})
            pangenome_id                : pangenome identifier
            pangenome_name              : pangenome display name
        Returns:
            Pangenome                   : KBaseGenomes.Pangenome like object (see type spec) 
                                          https://narrative.kbase.us/#spec/type/KBaseGenomes.Pangenome
    '''
    Pangenome = {}

    Pangenome['genome_refs'] = [tup[0] for tup in path_to_ref_and_ID_pos_dict.values()]
    Pangenome['id'] = pangenome_id
    Pangenome['name'] =  pangenome_name
    Pangenome['type'] = None

    OrthologFamilyList = []

    consistent_cols  = ['Gene','Non-unique Gene name','Annotation','No. isolates',\
                        'No. sequences','Avg sequences per isolate','Genome Fragment',\
                        'Order within Fragment','Accessory Fragment','Accessory Order with Fragment',\
                        'QC','Min group size nuc','Max group size nuc','Avg group size nuc']
    # load gene_pres_abs data and format it into Pangenome dictionary
    df = pd.read_csv(gene_pres_abs)
    # ignore errors to only drop columns that are in the dataframe

    cols = set(df.columns.values) - set(consistent_cols)

    # map the path_to_ref_and_ID_pos_dict to the columns names in the dataframe.
    # and make the mapping between gen_ids and gff_ids
    col_to_ref = {}
    for col in cols:
        start_len = len(col_to_ref)
        for path in path_to_ref_and_ID_pos_dict:
            if col == get_col_name_from_path(path):

                genome_ref, gen_id_to_pos, gen_ids = path_to_ref_and_ID_pos_dict[path]
                # now we would also like to use the information contained here to
                # create an output id to genome_object id mapping
                gff_id_to_gen_id = {}
                gff_ids = list(set(df[col].tolist()))
                for gff_id in gff_ids:
                    if not pd.isnull(gff_id):
                        gff_id = filter_gff_id(toString(gff_id))
                        i = find_pair(gff_id, gen_ids)
                        if i == None:
                            gff_id_to_gen_id[gff_id] = gff_id
                        else:
                            val = gen_ids.pop(i)
                            gff_id_to_gen_id[gff_id] = val[0]

                col_to_ref[col] = gff_id_to_gen_id, genome_ref, gen_id_to_pos

                break
        if len(col_to_ref) == start_len:
            colnames = []
            for p in path_to_ref_and_ID_pos_dict:
                colname = get_col_name_from_path(p)
                colnames.append((p, colname))
            raise ValueError(f"could not find file name match for {col} column"
                             f". should be in: {colnames}")

    df = df.where((pd.notnull(df)), None)

    # now we construct pangenome object
    for i, row in df.iterrows():
        OrthologFamily = {}

        # put in standard arguments for an OrthologFamily as found in Pangenome Spec file.
        OrthologFamily['id'] = row['Gene']
        OrthologFamily['type'] = None
        OrthologFamily['function'] = row['Annotation']  # should be the gene function
        OrthologFamily['md5'] = None
        OrthologFamily['protein_translation'] = None

        orthologs = []

        for col in cols:
            row_gff_id = row[col]
            if not pd.isnull(row_gff_id):
                # find if the gff_id is in fact multiple gff_id's tab delimited
                if '\t' in row_gff_id:
                    gff_ids = row_gff_id.split('\t')
                else:
                    gff_ids = [row_gff_id]

                gff_id_to_gen_id, genome_ref, gen_id_to_pos = col_to_ref[col]
                max_pos = max(list(gen_id_to_pos.values())) + 1
                for gff_id in gff_ids:
                    gff_id = filter_gff_id(toString(gff_id))

                    # check if gff_id in gff_id_to_gen_id
                    gene_id = None
                    if gff_id not in gff_id_to_gen_id:
                        if '___' in gff_id:
                            # chop off extra identifier if it exists
                            gff_id = gff_id.split('___')[0]
                            gff_id = filter_gff_id(gff_id)
                            if gff_id not in gff_id_to_gen_id:
                                # instead of throwing an error, we check if theres a pair in
                                # gen_ids, if not just use the same ID.
                                gene_id = gff_id
                                #gene_id = find_pair(gff_id, remaining_gen_ids)
                                # gff_id_error(gff_id, gff_id_to_gen_id, col)
                            else:
                                gene_id = gff_id_to_gen_id[gff_id]
                        else:
                            gene_id = gff_id
                            # gene_id = find_pair(gff_id, remaining_gen_ids) 
                            # gff_id_error(gff_id, gff_id_to_gen_id, col
                    else:
                        gene_id = gff_id_to_gen_id[gff_id]
                    if gene_id in gen_id_to_pos:
                        feature_pos = gen_id_to_pos[gene_id]
                    else:
                        feature_pos = max_pos
                        max_pos+=1
                    orthologs.append([gene_id, feature_pos, genome_ref])
        OrthologFamily['orthologs'] = orthologs

        OrthologFamilyList.append(OrthologFamily)

    Pangenome['orthologs'] = OrthologFamilyList

    return Pangenome


def upload_pangenome(cb_url, scratch, Pangenome, workspace_name, pangenome_name):
    """
    params:
        cb_url         : callback url
        scratch        : folder path to Pangenome object 
        pangenome      : KBaseGenomes.Pangenome like object
        workspace_name : workspace name
        pangenome_name : Pangenome display name
    Returns:
        pangenome_ref: Pangenome workspace reference
        pangenome_info: info on pangenome object
    """
    dfu = DataFileUtil(cb_url)
    meta = {}
    hidden = 0

    # dump pangenome to scratch for upload
    # data_path = os.path.join(scratch, pangenome_name + '.json')
    # json.dump(pangenome, open(data_path, 'w'))

    if isinstance(workspace_name, int) or workspace_name.isdigit():
        workspace_id = workspace_name
    else:
        workspace_id = dfu.ws_name_to_id(workspace_name)

    save_params = {
        'id': workspace_id,
        'objects': [{
            'type': 'KBaseGenomes.Pangenome',
            'data': Pangenome,
            'name': pangenome_name,
            'meta': meta,
            'hidden': hidden
        }]
    }
    print(f"Saving Pangenome object to {workspace_id}")
    try:
        info = dfu.save_objects(save_params)[0]
    except:
        raise RuntimeError("this is the attempted saved Pangenome: \n"
                            f"{Pangenome}")

    ref = "{}/{}/{}".format(info[6], info[0], info[4])
    print(f"Pangenome saved to {ref}")

    return {
        'pangenome_ref': ref,
        'pangenome_info': info
    }


def roary_report(cb_url, scratch, workspace_name, sum_stats, gene_pres_abs,
                 pangenome_ref, conserved_vs_total_graph, unique_vs_new_graph):
    """
    params:
        cb_url         : callback url
        workspace_name : name of the workspace
        sum_stats      : summary_statistics.txt file output from Roary
        pangenome_ref  : reference to the pangenome object, or None
    Returns:
        report_name : name of report object  
        report_ref  : reference to report object in workspace
    """
    report_name = 'Roary_report_'+str(uuid.uuid4()) 
    dfu = DataFileUtil(cb_url)

    # Convert output files to HTML
    html_output = format_output_html(sum_stats, gene_pres_abs)

    file_dir = os.path.join(scratch, report_name)
    os.mkdir(file_dir)

    copyfile(conserved_vs_total_graph, os.path.join(file_dir, "conserved_vs_total_genes.png"))
    copyfile(unique_vs_new_graph, os.path.join(file_dir, "unique_vs_new_genes.png"))

    html_path = os.path.join(file_dir, 'output.html')
    with open(html_path, 'w') as f:
        f.write(html_output) 


    html_link = {
        'path': file_dir,
        'name':'output.html',
        # 'label':'Summary_Statistics',
        'description':'Roary Gene Statistics html report'
    }

    csv_link = {
        'path':gene_pres_abs,
        'name':'gene_presence_absence.csv',
        'description':"Data Table of Gene Presence, Gene Absence and other information"
    }

    # photo_link_1 = {
    #     'path': conserved_vs_total_graph,
    #     'name':'conserved_vs_total_genes.png',
    #     # 'label':'Conserved_vs_total_genes_graph',
    #     'description':'Graph of conserved genes vs. total genes'
    # }
    # photo_link_2 = {
    #     'path': unique_vs_new_graph,
    #     'name':'unique_vs_new_genes.png',
    #     # 'label':'unique_vs_new_genes_graph',
    #     'description':'Graph of unique genes vs new genes'
    # }

    report_client = KBaseReport(cb_url)
    if pangenome_ref is not None:
        report = report_client.create_extended_report({
            'direct_html_link_index':0,
            'html_links':[html_link],
            'file_links':[csv_link],  #, photo_link_1, photo_link_2],
            'workspace_name': workspace_name,
            'report_object_name': report_name,
            'objects_created': [{
                'ref':pangenome_ref,
                'description':"Pangenome Object"
            }]
        })
    else:
        report = report_client.create_extended_report({
            'direct_html_link_index':0,
            'html_links':[html_link],
            'file_links':[csv_link, photo_link_1, photo_link_2],
            'workspace_name': workspace_name,
            'report_object_name': report_name
        })
    return {
        'report_name':report['name'],
        'report_ref': report['ref']
    }
