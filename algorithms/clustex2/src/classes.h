#include <map>
#include <string>
#include <set>
#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>

#include "Eigen/Dense"

extern "C" {
#include "include/f2c.h"
#include "include/clapack.h"
};


using namespace Eigen;
using namespace std;

typedef long int integer;
typedef double doublereal;

static omp_lock_t lock;

#define VERSION "2.0.0"
#define SOFTWARE "ClustEx"

class PPI;
class Gene_list;
class Parameter;
class Diff_modules;


class PPI{
public:
	int ppi_size;
	PPI(ifstream &ifs);
	set<string> ppi_gene();
	map<string, map<string, double> > get_co_express();
	map<string, set<string> > get_ppi();
	void ppi_filter();
	void diff_kern_calc( MatrixXd &simi_mat, Gene_list &gene_list, Parameter &params );
	void print_diff_kern(MatrixXd &simi_mat, Parameter &params);// print out the diffusion kernel
	void read_diff_kern(MatrixXd &simi_mat, Parameter&params);//read the diffusion kernel
	void info();
private:
	map<string, map<string,double> > _co_express;
	map<string, set<string> > _ppi_net;
};

class Gene_list{
public:
	Gene_list( ifstream &ifs, Parameter &params);
	set<string> get_seeds();
	void seeds_print();
	map<string, double> get_fold_change();
	void gene_overlap( PPI &network, Parameter &params);
	void random_walk( PPI &network, Parameter &params);
	void info( PPI &network );
	void print_density(Parameter &params);
	void read_density(Parameter &params );
private:
	map<string, double> _fold_change;
	set<string> _seed_gene;
};

class Parameter{
public:
	//threshold
	double restart_prob;
	double diff_beta;
	double seed_gene_thresh;
	double simi_thresh;
	bool simi_thresh_input;
	//double merge_time_cutoff;
	double merge_gene_cutoff;
	int parameter_level;

	//input file
	string gene_file;
	string network_file;
	bool gene_file_input;
	bool network_file_input;
	string density_file;
	string diff_kern_file;
	bool density_file_input;
	bool diff_kern_file_input;
	//step choice
	bool diff_kern_step;
	bool gene_density_step;
	bool parameter_calc_step;
	bool show_clustering_process;
	bool clustering_step;
	//function choice
	bool out_middle_results;
	//job ID
	string job_id;
	bool job_id_input;
	
	
	Parameter();

	void print_parameter();
};


class Diff_modules{
public: 
	void simi_thresh_calc( MatrixXd &simi_mat, PPI &network, Gene_list &gene_list, Parameter &params);
	void hma_origin(  MatrixXd &simi_mat, Gene_list &gene_list, double simi_thresh);// still not defined
	void hma_process( MatrixXd &simi_mat, Gene_list &gene_list, double simi_thresh, ofstream &ofs);
	void hma_merge_stop( MatrixXd &simi_mat, Gene_list &gene_list, Parameter &params, PPI &network);//still not used by simi_thresh_vec
	void hma_merge_fix_neighborhood_process(MatrixXd &simi_mat, Gene_list &gene_list, Parameter &params, PPI &network);
	void hma_merge_fix_neighborhood(MatrixXd &simi_mat, Gene_list &gene_list, Parameter &params, PPI &network);
	
	//hma_tree();
	void print_clusters( Gene_list &gene_list, PPI &network,Parameter &params);
	//for debug
	void clustering_results_check();
	//
	
private:
	map<int, set<string> > clusters;
	map<string, int> gene_of_cluster;
	map<int , vector<int> > clustering_process;// which is used during clustering processing
};



void replace_all( string &str, const string &old_val, const string &new_val);


void string_delimit(string &text, const char sep, vector<string> &words);

typedef pair<string, double> string_double;
inline double cmp( const string_double &x, const string_double &y){
	return x.second > y.second;
}

void help();

Parameter parameter_manipulate(int argc, char *argv[]);

template <typename inputIter , typename elemType>
void mean( inputIter first, inputIter last, elemType &result){
	elemType mean_val = 0;
	int ix=0;
	inputIter first_prime = first;
	for( ; first_prime != last; first_prime++){
		mean_val += *first_prime;
		ix++;
	}
	if( ix == 0 ){
		cout << "Wrong with the input of calculating the mean value!\n";
		cout  << mean_val << "\t";
		for( first_prime = first; first_prime != last; first++ )
			cout << *first_prime << "\t";
		cout << endl;
	}
	result = mean_val / ix;
}
