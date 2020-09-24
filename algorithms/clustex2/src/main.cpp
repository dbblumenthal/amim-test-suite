#include "classes.h"

int main(int argc, char *argv[]){
	Parameter new_param = parameter_manipulate( argc, argv);

	ifstream ifs_ppi( (new_param.network_file).c_str() );
	PPI new_net( ifs_ppi );
	new_net.ppi_filter();
	new_net.info();


	ifstream ifs_gene_list( (new_param.gene_file).c_str() );
	Gene_list new_gene_list(ifs_gene_list, new_param);
	new_gene_list.gene_overlap( new_net , new_param);
	//for debug
	//new_gene_list.seeds_print();
	//new_net.info();
	//new_gene_list.info(new_net);
	//
	integer N = ( new_net.ppi_gene() ).size();
	
	MatrixXd diff_kern( N, N );
	
	//input the diffusion kernel file if there are any
	if( new_param.diff_kern_file_input == true )
		new_net.read_diff_kern( diff_kern, new_param);
	if( new_param.density_file_input == true )
		new_gene_list.read_density( new_param);

	if( new_param.gene_density_step == true ){
		new_gene_list.random_walk(new_net, new_param);
		if( new_param.out_middle_results == true)
			new_gene_list.print_density( new_param);
	}

	if ( new_param.diff_kern_step == true ){
		new_net.diff_kern_calc(diff_kern, new_gene_list, new_param);
		if( new_param.out_middle_results == true)
			new_net.print_diff_kern(diff_kern, new_param);
	}

	if( new_param.parameter_calc_step == true ){
		Diff_modules new_clustering;
		new_clustering.simi_thresh_calc(diff_kern, new_net, new_gene_list,new_param);
	}
	
	if( new_param.show_clustering_process == true ){
		Diff_modules new_clustering;
		new_clustering.hma_merge_fix_neighborhood_process(diff_kern, new_gene_list, new_param, new_net );
	}
	
	if( new_param.clustering_step == true ){
		Diff_modules new_clustering;
		if( new_param.simi_thresh_input == false ){
			new_clustering.hma_merge_fix_neighborhood(diff_kern, new_gene_list, new_param, new_net);
		}
		else if( new_param.simi_thresh_input == true ){
			new_clustering.hma_merge_stop(diff_kern, new_gene_list, new_param, new_net);
		}
		//for debug
		//new_clustering.clustering_results_check();
		//
		new_clustering.print_clusters( new_gene_list, new_net, new_param);
	}
	
	


	return 0;
}

