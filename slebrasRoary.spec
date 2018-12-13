/*
A KBase module: slebrasRoary
TODO change name to Roary
*/

module slebrasRoary {
	/* roary input*/
  /* TODO add comments about each of these params and which ones are required */
    typedef structure {
        string workspace_name;
        string ref;
        string pangenome_name;
        int blast_p_percentage;
        int max_num_clusters;
        int percent_genes_for_core;
    } RoaryParams;

    /* roary output */
    typedef structure {
        string report_name;
        string report_ref;
    } RoaryResults;

    /*
        TODO update this comment
        This example function accepts any number of parameters and returns results in a KBaseReport
    */
    funcdef run_slebrasRoary(RoaryParams params) returns (RoaryResults output) authentication required;

};
