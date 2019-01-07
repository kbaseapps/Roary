/*
A KBase module: Roary
*/

module Roary {
    /* roary input
    params:
        workspace_name - name of the narrative workspace 
        ref - KBaseSearch.GenomeSet or KBaseSets.GenomeSet workspace reference
        pangenome_name - [OPTIONAL] name of optional Pangenome output
        blast_p_percentage - minimum percentage identity for blastp 
        max_num_clusters - maximum number of clusters
        percent_genes_for_core - percentage of isolates a gene must be in to be core
    */
    typedef structure {
        string workspace_name;
        string ref;
        string pangenome_name;
        int blast_p_percentage;
        int max_num_clusters;
        int percent_genes_for_core;
    } RoaryParams;

    /* roary output
    params:
        report_name - name of the outgoing Report object
        report_ref - workspace reference of the report object
    */
    typedef structure {
        string report_name;
        string report_ref;
    } RoaryResults;

    /*
        This example function accepts the parameters described in RoaryParams and returns results in RoaryResults inside a KBaseReport.
    */
    funcdef run_Roary(RoaryParams params) returns (RoaryResults output) authentication required;

};
