#include "classes.h"

void string_delimit(string &text, const char sep, vector<string> &words)
{
	/*******************************************************************
	Split a string according to a separator character into a vector
	of strings.
	*******************************************************************/
	words.clear();
	unsigned int i,j;
	i=0;
	while(i<text.size())
	{
		j=text.find_first_of(sep, i);
		if (j>text.size())
		{
			j=text.size();
		}
		if (j!=i)
		{
			words.push_back(text.substr(i, j-i));
		}
		i=j+1;
	}
}


void replace_all(string &str, const string &old_val, const string &new_val){
	/***********************
	function: replace all the &old_val in the &str with the &new_val
	build time: 2013/10/21
	***********************/
	while(true){
		string::size_type pos(0);
		if( (pos = str.find( old_val)) != string::npos)
			str.replace(pos, old_val.length(), new_val);
		else 
			break;
	}
}


PPI::PPI(ifstream &ifs){
	/****
	modification time: 2013/08/15
	function: constructor of PPI network
	requirement:
	(1) input file format: 
	header line: Gene 1 TAB Gene 2 TAB ( TAB weight)
	other lines...
    (2) when the format is wrong, stop clustex2 and call help();
	the function is not completed, leaving part (2) unfinished
	****/
	if( !ifs ){
		cerr << "Did not open up the file containing PPI network, please check again!\n";
		exit(-1);
	}
	string temp_str;
	bool _weighted = false;
	vector<string> header_line;
	getline(ifs, temp_str);
	if(temp_str[temp_str.length()-1] == '\r')
		temp_str = temp_str.substr(0, temp_str.length()-1);
	string_delimit( temp_str, '\t', header_line);
	if( header_line.size() == 3)
		_weighted = true;
	else if( header_line.size() == 2)
		_weighted = false;
	else{
		cerr << "Wrong format of network file!\n";
		exit(-1);
	}
	while( getline(  ifs , temp_str ) ){
		if(temp_str[temp_str.length()-1] == '\r')
			temp_str = temp_str.substr(0, temp_str.length()-1);
		vector<string> str_vec;
		string_delimit(temp_str, '\t', str_vec);
		//neglect self-loop
		if( str_vec[0] == str_vec[1] )
			continue;
		else{
			//first add this edge in the &ppi
			map<string, set<string> >::iterator ppi_ptr;
			ppi_ptr = _ppi_net.find( str_vec[0]);
			if( ppi_ptr != _ppi_net.end() )
				ppi_ptr->second.insert( str_vec[1] );
			else{
				set<string>  temp_set; 
				temp_set.insert( str_vec[1] );
				_ppi_net.insert( map<string, set<string> >::value_type( make_pair( str_vec[0] , temp_set) ) );
			}
			ppi_ptr = _ppi_net.find( str_vec[1] );
			if( ppi_ptr != _ppi_net.end() )
				ppi_ptr->second.insert( str_vec[0] );
			else{
				set<string> temp_set;
				temp_set.insert( str_vec[0] );
				_ppi_net.insert( map<string, set<string> >::value_type( make_pair(str_vec[1] , temp_set) ) );
			}
			//second add the weight in &co_express
			if( _weighted == false ){
				map<string, map<string, double> >::iterator corr_express_iter;
				corr_express_iter = _co_express.find( str_vec[0] );
				if( corr_express_iter != _co_express.end() ){
					//map<string, double> temp_map;
					//temp_map.insert( make_pair( str_vec[1], atof( str_vec[2].c_str() ) ) );
					//_co_express[ corr_express_iter->first ].insert( make_pair( str_vec[1], abs(atof( str_vec[2].c_str() ) ) ) );
					_co_express[ corr_express_iter->first ].insert( map<string, double>::value_type( make_pair( str_vec[1], 1.0 ) ) );
				}
				else{
					map<string, double> temp_map;
					//temp_map.insert( make_pair( str_vec[1], abs(atof( str_vec[2].c_str()) ) ) );
					temp_map.insert( map<string, double>::value_type( make_pair( str_vec[1], 1.0 ) ) );
					_co_express.insert( map<string, map<string, double> >::value_type( make_pair( str_vec[0], temp_map) ) );
				}
				corr_express_iter = _co_express.find( str_vec[1] );
				if( corr_express_iter != _co_express.end() ){
					//map<string, double> temp_map;
					//temp_map.insert( make_pair( str_vec[0], atof(str_vec[2].c_str() ) ) );
					//_co_express[ str_vec[1] ].insert( make_pair( str_vec[0], abs(atof(str_vec[2].c_str() )) ) );
					_co_express[ str_vec[1] ].insert( map<string,double>::value_type( make_pair( str_vec[0], 1 ) ) );
					//corr_express[ str_vec[1] ].insert( temp_map );
				}
				else{
					map<string, double> temp_map;
					//temp_map.insert( make_pair( str_vec[0], abs(atof( str_vec[2].c_str() ) )) );
					temp_map.insert( map<string, double>::value_type( make_pair( str_vec[0], 1 ) ) );
					_co_express.insert( map<string, map<string,double> >::value_type( make_pair( str_vec[1] , temp_map ) ) );
				}
			}
			else{
				map<string, map<string, double> >::iterator corr_express_iter;
				stringstream sstrm;
				sstrm << str_vec[2];
				double temp_val;
				sstrm >> temp_val;
				corr_express_iter = _co_express.find( str_vec[0] );
				if( corr_express_iter != _co_express.end() ){
					_co_express[ corr_express_iter->first ].insert( map<string, double>::value_type( make_pair( str_vec[1], abs(temp_val) ) ) );
				}
				else{
					map<string, double> temp_map;
					temp_map.insert(  map<string, double>::value_type( make_pair( str_vec[1], abs( temp_val ) ) ) );
					_co_express.insert( map<string, map<string,double> >::value_type( make_pair( str_vec[0], temp_map) ) );
				}
				corr_express_iter = _co_express.find( str_vec[1] );
				if( corr_express_iter != _co_express.end() ){
					_co_express[ str_vec[1] ].insert( map<string, double>::value_type( make_pair( str_vec[0], abs(temp_val ) ) ) );
				}
				else{
					map<string, double> temp_map;
					temp_map.insert( map<string, double>::value_type( make_pair( str_vec[0], abs( temp_val ) ) ) );
					_co_express.insert( map<string, map<string,double> >::value_type( make_pair( str_vec[1] , temp_map ) ) );
				}
			}
			
		}
	}
	ppi_size = _ppi_net.size();

}

void PPI::ppi_filter( ){
	/******************
	build time: 2013/08/15
	modification time: 2013/08/19
	function: filter the PPI network, delete isolated genes
	modification: a public ppi_size is added
	******************/
	map<string, set<string> >::iterator _ppi_iter = _ppi_net.begin();
	while( _ppi_iter != _ppi_net.end() ){
		if( _ppi_iter->second.size() == 0){
			map<string, set<string> >::iterator _temp_iter = _ppi_iter;
			_ppi_iter++;
			_co_express.erase( _temp_iter->first );
			_ppi_net.erase( _temp_iter->first );
		}
		else{
			_ppi_iter++;
		}
	}
	ppi_size = _ppi_net.size();
}

set<string> PPI::ppi_gene(){
	/**********
	build time: 2013/08/15
	function: obtain the genes of the ppi network
	**********/
	set<string> _genes;
	for(map<string, set<string> >::iterator temp_iter = _ppi_net.begin(); temp_iter != _ppi_net.end(); temp_iter++)
		_genes.insert( temp_iter->first );
	
	return _genes;
}

map<string, map<string,double> > PPI::get_co_express(){
	/***************
	build time: 2013/08/15
	function: acquire the input expression correlation
	***************/
	return _co_express;
}

map<string, set<string> > PPI::get_ppi(){
	/********************
	build time: 2013/08/16
	function: get the network without expression correlation
	********************/
	return _ppi_net;
}

void PPI::diff_kern_calc(MatrixXd &simi_mat, Gene_list &gene_list, Parameter &params){
	/******************
	build time: 2013/08/16
	modification time: change 0.01 to params.diff_beta
	function: calculate the diffusion kernel of the network
	equation: 
	exp( beta * H) = T * exp( beta * Q) * inv(T)
	T is eigenvector matrix and Q is eigenvalue diagonal matrix
	Q is the negative Laplacian matrix, Q = A - D
	A is the adjacent matrix and D is the degree matrix
	*******************/
	cout << "Start to calculate the diffusion kernel of the gene network (-D step)....";
	map<string, double> fold_change = gene_list.get_fold_change();
	integer N = fold_change.size();
	map<string, double>::iterator fc_iter = fold_change.begin();
	for(int i = 0; i < N; i++){
		map<string, double>::iterator temp_fc_iter = fold_change.begin();
		for(int j = 0; j < N; j++){
			map<string, double>::iterator temp_iter = _co_express[fc_iter->first].find( temp_fc_iter->first);
			if( temp_iter != _co_express[fc_iter->first].end() )
				simi_mat(i,j) = abs(temp_iter->second);
			else
				simi_mat(i,j) = 0;
			temp_fc_iter++;
		}
		fc_iter++;
	}

	VectorXd col_sum = simi_mat.colwise().sum();
	for( int i = 0; i < N; i++)
		simi_mat(i,i) = -col_sum(i);
	//decomposition
	//SelfAdjointEigenSolver<MatrixXd> eigensolver(simi_mat);
	//VectorXd temp_vec = eigensolver.eigenvalues();
	//use lapack
	char JOBZ = 'V';
	char UPLO = 'U';
	integer LDA = N;
	double *eigenvalue = new double[N];
	integer LWORK = 3 * N -1;
	doublereal *WORK = new doublereal[LWORK];
	integer INFO;
	dsyev_( &JOBZ, &UPLO, &N, simi_mat.data(), &LDA, eigenvalue, WORK, &LWORK, &INFO);
	delete WORK;

	if( INFO == 0 ){
		VectorXd temp_vec(N);
		for( int i = 0; i < N; i++)
			temp_vec(i) = exp( params.diff_beta * eigenvalue[i] );
			simi_mat = simi_mat * ( temp_vec.asDiagonal() ) * ( simi_mat.transpose() );
		
		cout << "done!\n";
		delete eigenvalue;
	}
	else{
        cout << "ClustEx 2.0 is unable to calculate your network's diffusion kernel!\n";
		exit(-1);
	}
}

void PPI::print_diff_kern(MatrixXd &simi_mat, Parameter &params){
	/***************
	build time: 2013/08/17
	function: print out the diffusion kernel
	****************/
	string out_file;
	stringstream sstrm;
	string restart_prob_str;
	string diff_beta_str;
	sstrm << params.diff_beta;
	sstrm >> diff_beta_str;
	sstrm.clear();
	sstrm << params.restart_prob;
	sstrm >> restart_prob_str;
	out_file = params.job_id + "_diffusion_kernel_" + diff_beta_str + ".txt";
	ofstream ofs( (out_file).c_str() );
	if( !ofs ){
		cerr << "Cannot open the file to output the diffusion kernel!\n";
	}
	//int N = _ppi_net.size();
	for( int i = 0; i < ppi_size; i++){
		for(int j = 0; j < ppi_size; j++)
			ofs << simi_mat(i,j) << "\t";
		ofs << endl;
	}
	ofs.close();
}

void PPI::read_diff_kern(MatrixXd &simi_mat, Parameter &params){
	/****************
	build time: 2013/08/17
	function: read the diffusion kernel
	****************/

	int N = _ppi_net.size();
	string temp_str;
	int i = 0;
	//int N = _ppi_net.size();
	ifstream ifs( (params.diff_kern_file).c_str() );
	if( !ifs ){
		cerr << "Unable to open the diffusion kernel file!\n";
		help();
		exit(-1);
	}
	while( getline(ifs, temp_str) ){
		if( temp_str[temp_str.length() -1] == '\r' )
			temp_str = temp_str.substr(0, temp_str.length()-1);
		vector<string> str_vec;
		string_delimit( temp_str, '\t', str_vec );
		if( str_vec.size() != N){
			cerr << "Wrong with the diffusion kernel matrix file!\n";
			exit(-1);
		}
		stringstream sstrm;
		for( int j = 0; j < N; j++){
			sstrm.clear();
			sstrm << str_vec[j];
			double temp_val;
			sstrm >> temp_val;
			simi_mat(i,j) = temp_val;
		}
		i++;
	}
	/*debug info
	ofstream ofs("diffusion_kernel.txt");
	for( int i = 0; i < N; i++){
		for( int j = 0; j < N; j++)
			ofs << simi_mat(i,j) << '\t';
		ofs << endl;
	}
	ofs.close();
	*/

	/*
	if( i != ppi_size){
		cerr << "Wrong with the diffusion kernel matrix file!\n";
		exit(-1);
	}
	*/
}

void PPI::info(){
	/***************
	build time: 2013/08/18
	function: display the information of a PPI, including
	(1) how many genes in the map<string, set<string> > _ppi_net
	(2) how many genes in the map<string, map<string, double> > _co_express
	(3) check the above two data structures whether they are consensus
	***************/
	cout << "There are " << _ppi_net.size() << " genes in the network.\n";
	/*
	cout << "There are " << _co_express.size() << " genes in the weighted network.\n";
	for(map<string, set<string>  >::iterator net_iter = _ppi_net.begin(); net_iter != _ppi_net.end(); net_iter++){
		set<string> temp_set_1;
		set<string> temp_set_2;
		map<string, double >::iterator temp_iter = _co_express[net_iter->first].begin();
		for( ; temp_iter != _co_express[net_iter->first].end(); temp_iter++ ){
			if( temp_iter->second < 0 )
				cerr << "A correlation is smaller than 0!\n";
			temp_set_1.insert( temp_iter->first );
		}
		for( set<string>::iterator temp_iter = net_iter->second.begin(); temp_iter != net_iter->second.end(); temp_iter++){
			temp_set_2.insert( *temp_iter );
		}
		if( temp_set_1 != temp_set_2)
			cerr << "The weighted network and the unweighted network are not the same with the gene " << net_iter->first << endl;
	}
	cout << "Done checking the consensus of the weighted network and unweighted network!\n";
	*/
	//how many edges in the network
	int edge_num = 0;
	for( map<string, set<string> >::iterator ppi_iter = _ppi_net.begin(); ppi_iter != _ppi_net.end(); ppi_iter++){
		for (set<string>::iterator set_iter = ppi_iter->second.begin(); set_iter != ppi_iter->second.end(); set_iter++)
			edge_num++;
	}
	cout << "There are " << edge_num/2 << " edges in the network!\n";



}

Gene_list::Gene_list(ifstream &ifs, Parameter &params){
	/*****************
	build time: 2013/08/15
	modification time: 2013/11/15
	modification purpose: if have "\r\n" as the end, delete "\r"
	function: 
	(1) input the seed gene list
	(2) get the initial fold change
	requirement of the input
	(1) header line: Gene (TAB fold change)
	(2) other lines....(format same with header line)
	*****************/
	if( ! ifs ){
		cerr << "Did not open up the gene list file!\n";
		exit(-1);
	}
	string temp_str;
	getline(ifs, temp_str);
	if( temp_str[temp_str.length()-1] == '\r')
		temp_str = temp_str.substr(0, temp_str.length()-1);
	vector<string> header_line;
	string_delimit( temp_str, '\t', header_line);
	bool _weighted;
	if( header_line.size() == 1 )
		_weighted = false;
	else if( header_line.size() == 2)
		_weighted = true;
	else{
		cerr << "The header line is " << temp_str <<  "and its size() is " << header_line.size() << " , the input format of the gene list file is incorrect!\n";;
		exit(-1);
	}

	bool _true_seeds = false;
	if( _weighted == true && params.seed_gene_thresh != 0 )
		_true_seeds = true;

	while( getline(ifs, temp_str) ){
		if( temp_str[temp_str.length()-1] == '\r' )
			temp_str = temp_str.substr(0, temp_str.length()-1);
		vector<string> str_vec;
		string_delimit( temp_str, '\t', str_vec);
		if( _weighted == false ){
			_fold_change.insert(map<string, double>::value_type(make_pair( str_vec[0], 1) ) );
			_seed_gene.insert( str_vec[0] );
		}
		else if( _weighted == true){
			stringstream transfer_str;
			double temp_val;
			transfer_str << str_vec[1];
			transfer_str >> temp_val;
			temp_val = abs(temp_val);
			_fold_change.insert( map<string, double>::value_type( make_pair( str_vec[0], temp_val ) ) );
			if( _true_seeds == false )
				_seed_gene.insert( str_vec[0] );
			else if( _true_seeds == true && temp_val >= params.seed_gene_thresh ){
				_seed_gene.insert( str_vec[0] );
			}
		}
	}

	if( _seed_gene.empty() == true ){
		cerr << "The seed gene threshold may be too large that no genes are considered to be seed genes!\n";
		help();
		exit(-1);
	}
}

set<string> Gene_list::get_seeds( ){
	/************
	build time: 2013/08/15
	function: get the seed genes
	*************/
	return _seed_gene;
}

map<string, double> Gene_list::get_fold_change(){
	/********
	build time: 2013/08/16
	function: get the fold change of each gene
	*********/
	return _fold_change;
}

void Gene_list::gene_overlap(PPI &network, Parameter &params){
	/************
	build time: 2013/08/16
	modification time: 2013/08/19; 2013/10/09; 2013/10/28
	function: 
	(1) delete seed genes that are not in the network
	(2) complete all genes including seed genes and network genes fold change
	(3) if _fold_change is larger than PPI, delete genes that are not in PPI
	(4) if there are no seed genes, then stop running and report the error
	************/
	set<string> gene_of_ppi = network.ppi_gene();
	//(1) delete seed genes that are not in the network
	string str;
    str = params.job_id + "_seed_genes_not_in_network.txt";
	//ofstream ofs("seed genes not in network.txt");
	ofstream ofs( str.c_str() );
    cout << "There are " << _seed_gene.size() << " seed genes in the initial seed gene list!\n";
	set<string>::iterator seed_gene_iter = _seed_gene.begin();
	while( seed_gene_iter != _seed_gene.end() ){
		if( gene_of_ppi.find( *seed_gene_iter ) == gene_of_ppi.end() ){
			set<string>::iterator temp_iter = seed_gene_iter;
			seed_gene_iter++;
			ofs << *temp_iter << endl;
			_fold_change.erase( *temp_iter );
			_seed_gene.erase( temp_iter );
		}
		else
			seed_gene_iter++;
	}
    cout << "Filtered with the gene network, there are " << _seed_gene.size() << " seed genes!\n";
	if( _seed_gene.size() == 0){
		cerr << "There are no seed genes and ClustEx 2.0 stops, please check if your gene name are consensus in both candidate genes list and gene network!\n";
		help();
		exit(-1);
	}
	//(2) complete all genes including seed genes and network genes fold change
	for( set<string>::iterator ppi_gene_iter = gene_of_ppi.begin(); ppi_gene_iter != gene_of_ppi.end(); ppi_gene_iter++){
		if( _fold_change.find( *ppi_gene_iter) == _fold_change.end() ){
			_fold_change.insert( map<string, double>::value_type(make_pair( *ppi_gene_iter, 0) ) );
		}
	}
	map<string, double>::iterator fold_iter = _fold_change.begin();
	while(  fold_iter != _fold_change.end() ){
		if( gene_of_ppi.find( fold_iter->first ) == gene_of_ppi.end() ){
			map<string, double>::iterator temp_iter = fold_iter;
			fold_iter++;
			_fold_change.erase( temp_iter );
		}
		else{
			fold_iter++;
		}
	}
}

void Gene_list::random_walk(PPI &network, Parameter &params){
	/*********************
	build time: 2013/08/16
	function: calculate the final probability distribution
	(1) 1/c * [ E- (1-c)*A ] * p =  p_0
	(2) A is column normalization
	(2) Eigen is used for matrix calculation and so thus clapack
	*********************/
    cout << "Start to calculate the gene importance (-G step)....";
	map<string, map<string,double> > co_express = network.get_co_express();
	integer N = _fold_change.size();
	//complete the initial 
	MatrixXd simi_mat(N,N);
	map<string, double>::iterator fc_iter = _fold_change.begin();
	for(int i = 0; i < N; i++){
		map<string, double>::iterator temp_fc_iter = _fold_change.begin();
		for(int j = 0; j < N; j++){
			map<string, double>::iterator temp_iter = co_express[fc_iter->first].find( temp_fc_iter->first);
			if( temp_iter != co_express[fc_iter->first].end() )
				simi_mat(i,j) = abs(temp_iter->second);
			else
				simi_mat(i,j) = 0;
			temp_fc_iter++;
		}
		fc_iter++;
	}
	
	VectorXd col_sum = simi_mat.colwise().sum();
	for( int i = 0; i < N; i++)
		//for each column
		for( int j = 0; j < N; j++)
			//for each line
			simi_mat(j,i) = simi_mat(j,i)/col_sum(i);
	simi_mat = -(1-params.restart_prob) * simi_mat;
	for( int i = 0; i < N; i++)
		simi_mat(i,i) += 1;
	simi_mat = simi_mat / params.restart_prob;

	VectorXd init_prob(N);
	fc_iter = _fold_change.begin();
	for( int i = 0; i < N; i++){
		init_prob(i) = abs(fc_iter->second);
		fc_iter++;
	}
	init_prob = init_prob / init_prob.sum();

	/*debug info
	cout << "Done preparing all the matrices and start to calculate the final probability distribution!\n";
	*/
	//VectorXd final_prob = simi_mat.colPivHouseholderQr().solve(init_prob);
	integer NRHS = 1;
	// N
	integer LDA = N;
	integer LDB = N;
	integer *IPIV = new integer[N];
	double *final_prob;
	final_prob = new double[N];
	for(int i = 0; i < N; i ++)
		final_prob[i] = init_prob[i];
	integer INFO;
	dgesv_(&N, &NRHS, simi_mat.data(), &LDA, IPIV, final_prob, &LDB, &INFO);
	delete IPIV;

	/*debug info
	cout << "Done calculating!\n";
	*/
	if( INFO == 0){
		cout << "done!\n";
		fc_iter = _fold_change.begin();
		for( int i = 0; i < N; i++){
			if( isnan(final_prob[i]) || isinf(final_prob[i])){
				cerr << "Cannot calculate the gene importance because improper gene-gene correlation matrix!\n";
				delete final_prob;
				exit(-1);
			}
			fc_iter->second = final_prob[i];
			//debug info
			init_prob(i) = final_prob[i];
			//
			fc_iter++;
		}
		//cout << "The sum prob is " << init_prob.sum() << endl;
		double temp_sum_0215 = init_prob.sum();
		if(isnan(temp_sum_0215) || isinf(temp_sum_0215) ){
			cerr << "Cannot calculate the gene importance because improper gene-gene correlation matrix!\n";
			delete final_prob;
			exit(-1);
		}
		delete final_prob;
	}
	else{
		cerr << "Cannot calculate the gene importance because improper gene-gene correlation matrix!\n";
		delete final_prob;
		exit(-1);
	}
	
	//cout << "The sum prob is " << final_prob.sum() << endl;
	//double relative_error = (simi_mat * final_prob - init_prob ).norm() / init_prob.norm();
	//cout << "The relative error is : " << relative_error << endl;

	
	/*for debug
	ofstream ofil("gene_density.txt");
	for(map<string, double>::iterator temp_iter = fold_change.begin(); temp_iter != fold_change.end(); temp_iter++)
		ofil << temp_iter->first << "\t" << temp_iter->second << "\n";
	ofil.close();
	*/
	
}

void Gene_list::seeds_print(){
	/**********
	function: output the seed genes
	build time: 2013/09/28
	***********/
	ofstream ofil("seed_genes_in_network.txt");
	ofil <<  "There are " << _seed_gene.size() << "genes\n";
	set<string>::iterator _iter = _seed_gene.begin();
	for( ; _iter != _seed_gene.end(); _iter++){
		ofil << *_iter << endl;
	}

}


void Gene_list::info( PPI &network ){
	/************
	build time: 2013/08/18
	function: check a Gene_list, including
	(1) _seed_genes exist in _fold_change and the PPI
	(2) _fold_change genes are the same with the PPI
	************/
	bool consensus = true;
	map<string, set<string> > ppi_net = network.get_ppi();
	bool fold_in_ppi = true;
	bool ppi_in_fold = true;
	cout << "Now there are " << _fold_change.size() << " genes which have initial  weights.\n";
	cout << "There are " << _seed_gene.size() << " seed genes.\n";
	for( map<string, double>::iterator fold_iter = _fold_change.begin(); fold_iter != _fold_change.end(); fold_iter++){
		if( ppi_net.find( fold_iter->first) == ppi_net.end() )
			fold_in_ppi = false;
	}
	for( map<string, set<string> >::iterator ppi_iter = ppi_net.begin(); ppi_iter != ppi_net.end(); ppi_iter++){
		if( _fold_change.find( ppi_iter->first) == _fold_change.end() )
			ppi_in_fold = false;
	}
	if( fold_in_ppi == false || ppi_in_fold == false ){
		cerr << "The ppi and the gene list are not consensus!\n";
		consensus = false;
	}
	if( consensus == false )
		cerr << "The genes with initial differential scores are not all in the gene network!\n";
	else
		cout << "The genes with initial differential scores are all in the gene network!\n";
	
	bool consensus_2 = true;
	for( set<string>::iterator temp_iter = _seed_gene.begin(); temp_iter != _seed_gene.end(); temp_iter++){
		if( _fold_change.find( *temp_iter ) == _fold_change.end() ){
			//cerr << "The seed gene list is not consensus with the fold change or the ppi\n";
			consensus_2 = false;
			cerr << "The gene is " << *temp_iter << endl;
		}
	}

	if( consensus_2 == false)
		cerr << "Not all the seed genes are in the network!\n";
	else
		cout << "All the seed genes are in the network!\n";

	if( consensus == true && consensus_2 == true )
		cout << "The seed gene list is part of the network!\n";
	else
		cout << "The seed gene list is not in consistent with the network!\n";

}

void Gene_list::print_density(Parameter &params){
	/***********
	build time: 2013/08/17
	function: print out the gene density score
	***********/
	string out_file;
	stringstream sstrm;
	sstrm << params.restart_prob;
	string restart_prob_str;
	sstrm >> restart_prob_str;
	string diff_beta_str;
	sstrm.clear();
	sstrm << params.diff_beta;
	sstrm >> diff_beta_str;
	out_file = params.job_id + "_gene_importance_" + restart_prob_str  + ".txt";
	ofstream ofs( (out_file).c_str() );
	ofs << "Gene" << '\t' << "density_score" << endl;
	for( map<string, double>::iterator temp_iter = _fold_change.begin(); temp_iter != _fold_change.end(); temp_iter++){
		ofs << temp_iter->first << '\t' << temp_iter->second << endl;
	}
	ofs.close();
}

void Gene_list::read_density(Parameter &params){
	/********
	build time: 2013/08/17
	function: read the gene density score
	requirment:
	(1) the input file MUST have a header line
	*********/

	ifstream ifs( (params.density_file).c_str() );
	if( !ifs){
		cerr << "Cannot open the gene importance file!\n";
		exit(-1);
	}
	string temp_str;
	getline( ifs, temp_str);//the header line
	while( getline( ifs, temp_str) ){
		if( temp_str[temp_str.length()-1] == '\r' )
			temp_str = temp_str.substr( 0, temp_str.length()-1 );
		vector<string> str_vec;
		string_delimit( temp_str, '\t', str_vec );
		stringstream sstrm;
		sstrm << str_vec[1];
		sstrm >> _fold_change[str_vec[0]];
	}
}


void Diff_modules::simi_thresh_calc( MatrixXd &simi_mat, PPI &network, Gene_list &gene_list, Parameter &params){
	/*****************
	build time: 2013/08/16
	modification time: 2013/11/05
	modification purpose: if ClustEx 2.0 calculates the diffusion kernel,it needs to read the diffusion kernel TXT here
	function: try to find out an appropriate set of parameters: similarity threshold and merge time cutoff
	*****************/
    cout << "Start to cluster genes with different neighborhoods (-P step)....";
	map<string,double> _fold_change = gene_list.get_fold_change();
	map<string, int> gene_order;
	int N = _fold_change.size();
	omp_lock_t lock;
	
	/************/
	if( params.diff_kern_step == true ){
		string infile;
		stringstream sstrm;
		string tmp_str;
		sstrm << params.diff_beta;
		sstrm >> tmp_str;
		infile = params.job_id + "_diffusion_kernel_" + tmp_str  + ".txt";
		ifstream ifs( infile.c_str() );
		//int N = _ppi_net.size();
		string temp_str;
		int i = 0;
		while( getline(ifs, temp_str) ){
			vector<string> str_vec;
			string_delimit( temp_str, '\t', str_vec );
			if( str_vec.size() != N){
				cerr << "Wrong with the diffusion kernel matrix file!\n";
				exit(-1);
			}
			//stringstream sstrm;
			for( int j = 0; j < N; j++){
				sstrm.clear();
				sstrm << str_vec[j];
				double temp_val;
				sstrm >> temp_val;
				simi_mat(i,j) = temp_val;
			}
			i++;
		}
	}

	//first load the order of simi_mat corresponding to fold_change
	
	int order = 0;
	for( map<string, double>::iterator gene_fold_iter = _fold_change.begin(); gene_fold_iter != _fold_change.end(); gene_fold_iter++){
		gene_order.insert( map<string,int>::value_type( make_pair( gene_fold_iter->first, order) ) );
		order++;
	}
	//store all the similarities, which will be sorted to find out 100 similarity thresholds
	map<string, set<string> > ppi = network.get_ppi();
	set<double> simi_set;
	for( map<string, set<string> >::iterator ppi_iter = ppi.begin(); ppi_iter != ppi.end(); ppi_iter++){
		set<string>::iterator set_iter = ppi_iter->second.begin();
		for( ; set_iter != ppi_iter->second.end(); set_iter++){
			simi_set.insert( simi_mat( gene_order[ppi_iter->first] , gene_order[*set_iter] ) );
		}
	}

	vector<double> simi_vec;
	for( set<double>::iterator set_iter = simi_set.begin(); set_iter != simi_set.end(); set_iter++)
		simi_vec.push_back( *set_iter );
	sort( simi_vec.begin(), simi_vec.end() );
	integer vec_length = simi_vec.size() -1;
	//debug info
	//cout << "\nThe vector length is " << vec_length << "\n";
	//
	vector<double> simi_thresh_vec;
	// 20 parameters, 0.05 to 0.95; 50 parameters, 0.02 to 0.98; 99 parameters
	if( params.parameter_level == 1){// for webserver
		for( double i = 0.95; i >= 0.695; i = i -0.05 ){
			simi_thresh_vec.push_back( simi_vec[integer( floor( i*vec_length) ) ] );
		}
	}
	else if( params.parameter_level ==2 ){
		for(double i = 0.98; i >= 0.695; i = i - 0.02 )
			simi_thresh_vec.push_back( simi_vec[integer( floor( i*vec_length) ) ] );
	}
	else if( params.parameter_level == 3){
		for( double i = 0.99; i >= 0.695; i = i - 0.01 )
			simi_thresh_vec.push_back( simi_vec[ integer(floor( i * vec_length )) ] );
	}
	else if( params.parameter_level == 4 ){
		for( double i = 0.99;  i>= 0.005; i = i - 0.01 )
			simi_thresh_vec.push_back( simi_vec[ integer(floor( i * vec_length )) ] );
	}
	
	//debug info
	//cout <<"The next vector length is " << simi_thresh_vec.size() << endl;
	//
	
	
	//now use each similarity threshold to clustering a result, and use this result to calculate corresponding parameter
	//string temp_str = params.network_file;
	string temp_str_2;
	temp_str_2 = params.job_id + "_neighborhood.txt";
	//temp_str += temp_str_2;
	//ofstream ofil_temp(temp_str_2.c_str(),ios_base::app);
	ofstream ofil_temp_1(temp_str_2.c_str());
	ofil_temp_1.close();
	ofstream ofil_temp(temp_str_2.c_str(),ios_base::app);
	ofil_temp << "similarity threshold" << '\t' <<"clustered genes" << '\t' << "clustered seed gene numbers" << '\t'
                       << "largest module size" <<'\t' << "seed genes in largest"  << '\t'
					   << "second largest module size" << '\t' << "seed gene in 2nd largest"
					   << '\t' << endl;
    //vector<double> cluster_efficient_vec;
    //map<double, vector<double> > cluster_efficient_map;
    //here is the place need to use OpenMp
    omp_init_lock(&lock);
#pragma omp parallel for
    for( int i = 0 ; i < simi_thresh_vec.size(); i++)
         hma_process(simi_mat, gene_list, simi_thresh_vec[i], ofil_temp);

	ofil_temp.close();
	cout << "done!\n";
    omp_destroy_lock(&lock);
}




void Diff_modules::hma_process(MatrixXd &simi_mat, Gene_list &gene_list, double simi_thresh, ofstream &ofs){
    /**********************************************
    build time: 2013/08/13,08/14
    modification time: 2013/08/16; 2013/09/14;2013/10/2; 2013/10/30;2013/11/15
    function: under each merge time,
    under one specific  similarity threshold, output the clustering process, mainly including under each merge time
    (1) clustered genes and seed genes, or the ratio
    (2) largest module size and seed genes in it, or the ratio
    (3) second largest module size and seed genes in it, or the ratio?
    modification details at 2013/09/14:
    (1) a largest module with its size and a second largest module with its size
    (2) every time a gene is merged, modify the largest module and second largest module
    ************************************************/

	set<string> seed_genes = gene_list.get_seeds();
    map<string,double> fold_change = gene_list.get_fold_change();
    map<string, int> simi_mat_order;
	omp_lock_t lock;
    int order = 0;
    for( map<string, double>::iterator fc_iter = fold_change.begin(); fc_iter != fold_change.end(); fc_iter++ ){
        simi_mat_order.insert( map<string,int>::value_type( make_pair( fc_iter->first, order) ) );
        order++;
    }
    // use HMA method to hierarchically clustering genes
    vector<string_double> fc_vec;
    for( map<string, double>::iterator fc_iter = fold_change.begin(); fc_iter != fold_change.end(); fc_iter++)
        fc_vec.push_back( make_pair( fc_iter->first, fc_iter->second ) );
    sort( fc_vec.begin(), fc_vec.end(), cmp );

    //_clusters_.clear();
    //_gene_of_cluster_.clear();
    //_clustering_process_.clear();
    map<int,set<string> > _clusters_;
    map<string,int> _gene_of_cluster_;
    map<int, vector<int> > _clustering_process_;



    vector<string> clustered_genes;
    // make the first module
    vector<string_double>::iterator fc_vec_iter = fc_vec.begin();
    set<string> temp_1;
    temp_1.insert( fc_vec_iter->first );
    int cluster_num = 1;
    _clusters_.insert( map<int,set<string> >::value_type(make_pair(cluster_num, temp_1)) );
    _gene_of_cluster_.insert( map<string,int>::value_type(make_pair(fc_vec_iter->first, cluster_num)) );
    clustered_genes.push_back( fc_vec_iter->first );

	//make the largest module and 2nd largest module
	int largest_module = cluster_num;
	int largest_module_size = 1;
	int largest_module_seeds = 0;
	int second_largest_module = 0;
	int second_largest_module_size = 0;
	int second_largest_module_seeds = 0;
	//clustered seed genes;
	int clustered_seeds = 0;
	if( seed_genes.find(fc_vec_iter->first) != seed_genes.end() )
	{
		largest_module_seeds++;
		clustered_seeds++;
	}
	//

    fc_vec_iter++;
    cluster_num++;
   
    
    /* go on to cluster  genes */

    for( ; fc_vec_iter != fc_vec.end(); fc_vec_iter++){
        vector<string> closest_genes;
        set<int>  closest__clusters_;
        for( int i = 0; i < clustered_genes.size(); i++){
            // the rule that a gene belongs to a cluster :
            // nearest one of its closest gene
            if( simi_mat( simi_mat_order[ fc_vec_iter->first] , simi_mat_order[clustered_genes[i]] ) >= simi_thresh ){//not less than
                closest_genes.push_back( clustered_genes[i] );
                closest__clusters_.insert( _gene_of_cluster_[ clustered_genes[i] ] );
            }
        }

		//is this gene a seed gene?
		if( seed_genes.find(fc_vec_iter->first) != seed_genes.end() )
			clustered_seeds++;

        if( closest_genes.empty() ){
            clustered_genes.push_back(  fc_vec_iter->first );
            /* the gene generate a new cluster */
            set<string> temp;
            temp.insert( fc_vec_iter->first);
            _clusters_.insert( map<int,set<string> >::value_type( make_pair(cluster_num, temp)) );
            _gene_of_cluster_.insert( map<string,int>::value_type( make_pair(fc_vec_iter->first, cluster_num) ) );

			if( _clusters_[cluster_num].size() > second_largest_module_size  ){
				second_largest_module = cluster_num;
				second_largest_module_size = 1;
				if( seed_genes.find( fc_vec_iter->first) != seed_genes.end() )
					second_largest_module_seeds = 1;
			}
			//record the current statistics; key: clustered genes; value : (1) clustered seeds (2) largest module size (3) largest module seeds (4) 2nd largest module (5) 2nd largest module seeds
			vector<int> curr_stat;
			curr_stat.push_back( clustered_seeds);
			curr_stat.push_back( largest_module_size);
			curr_stat.push_back( largest_module_seeds);
			curr_stat.push_back( second_largest_module_size );
			curr_stat.push_back( second_largest_module_seeds );
			int clustered_gene_num = clustered_genes.size();
			_clustering_process_.insert( map<int,vector<int> >::value_type( make_pair( clustered_gene_num, curr_stat) ) );

            cluster_num++;
        }
        else if( closest_genes.size() == 1 || closest__clusters_.size() == 1 ){
            clustered_genes.push_back( fc_vec_iter->first );
            /* only one genes are closest */
            ( _clusters_[ _gene_of_cluster_[ closest_genes[0] ] ] ).insert( fc_vec_iter->first );
            _gene_of_cluster_.insert( map<string,int>::value_type(make_pair(fc_vec_iter->first, _gene_of_cluster_[ closest_genes[0] ]) ) );

			//determine largest module and second largest module, perhaps mistakes
			// 1st situation: larger than the largest module
			int temp_size = _clusters_[ _gene_of_cluster_[ closest_genes[0] ] ].size();
			if(  temp_size > largest_module_size ){
				largest_module = _gene_of_cluster_[ closest_genes[0] ];
				largest_module_size = _clusters_[ _gene_of_cluster_[ closest_genes[0] ] ].size();
				vector<string> temp_vec;
				set_intersection( _clusters_[largest_module].begin(), _clusters_[largest_module].end(), seed_genes.begin(), seed_genes.end(), back_inserter(temp_vec) );
				largest_module_seeds = temp_vec.size();

				//second_largest_module_size = 0;
				int new_second_largest_module = 0;
				int new_second_largest_module_size = 0;
				for(map<int ,set<string> >::iterator clusters_iter = _clusters_.begin(); clusters_iter != _clusters_.end(); clusters_iter++ ){
					if(  (_clusters_[clusters_iter->first].size() < largest_module_size) && ( _clusters_[clusters_iter->first].size() > new_second_largest_module_size ) ){
						new_second_largest_module_size = _clusters_[clusters_iter->first].size();
						new_second_largest_module = clusters_iter->first;
					}
				}
				second_largest_module_size = new_second_largest_module_size;
				second_largest_module = new_second_largest_module;
				vector<string> temp_vec_2;
				set_intersection( _clusters_[second_largest_module].begin(), _clusters_[second_largest_module].end(), seed_genes.begin(), seed_genes.end(), back_inserter(temp_vec_2) );
				second_largest_module_seeds = temp_vec_2.size();
			}
			// 2nd situation: smaller than the largest but larger than the 2nd largest
			else if( ( temp_size < largest_module_size) && ( temp_size > second_largest_module_size ) ){
				if( second_largest_module != _gene_of_cluster_[ closest_genes[0] ] ){
					second_largest_module = _gene_of_cluster_[ closest_genes[0] ];
					second_largest_module_size = _clusters_[ _gene_of_cluster_[ closest_genes[0] ] ].size();
					vector<string> temp_vec;
					set_intersection( _clusters_[ _gene_of_cluster_[ closest_genes[0] ] ].begin(), _clusters_[ _gene_of_cluster_[ closest_genes[0] ] ].end(), seed_genes.begin(), seed_genes.end(), back_inserter(temp_vec) );
					second_largest_module_seeds = temp_vec.size();
				}
				else if ( second_largest_module == _gene_of_cluster_[closest_genes[0]] ){
					second_largest_module_size++;
					if( seed_genes.find(fc_vec_iter->first) != seed_genes.end() )
						second_largest_module_seeds++;
				}
			}
			else if( temp_size == largest_module_size ){
				second_largest_module = _gene_of_cluster_[closest_genes[0]];
				second_largest_module_size = temp_size;
				vector<string> temp_vec;
				set_intersection( _clusters_[ _gene_of_cluster_[closest_genes[0]]].begin(), _clusters_[ _gene_of_cluster_[closest_genes[0]]].end(), seed_genes.begin(), seed_genes.end(), back_inserter(temp_vec) );
				second_largest_module_seeds = temp_vec.size();
			}
			//record the current statistics; key: clustered genes; value : (1) clustered seeds (2) largest module size (3) largest module seeds (4) 2nd largest module (5) 2nd largest module seeds
			vector<int> curr_stat;
			curr_stat.push_back( clustered_seeds);
			curr_stat.push_back( largest_module_size);
			curr_stat.push_back( largest_module_seeds);
			curr_stat.push_back( second_largest_module_size );
			curr_stat.push_back( second_largest_module_seeds );
			int clustered_gene_num = clustered_genes.size();
			_clustering_process_.insert( map<int,vector<int> >::value_type( make_pair( clustered_gene_num, curr_stat) ) );

        }
        else{
            clustered_genes.push_back( fc_vec_iter->first );
            /* this gene is closed to several _clusters_ */
            //first merge these _clusters_
            //int module_merge = _gene_of_cluster_[ closest_genes[0] ];
            int module_merge = cluster_num;
            cluster_num++;
            set<int>::iterator closest__clusters__iter = closest__clusters_.begin();
            //start to merge
            for( ; closest__clusters__iter != closest__clusters_.end(); closest__clusters__iter++){
				if( *closest__clusters__iter == module_merge )
					continue;
                ( _clusters_[ module_merge ] ).insert( ( _clusters_[ *closest__clusters__iter ]).begin() , ( _clusters_[ *closest__clusters__iter ]).end() );
                for( set<string>::iterator temp_iter = ( _clusters_[ *closest__clusters__iter ]).begin(); temp_iter != ( _clusters_[ *closest__clusters__iter ]).end(); temp_iter++)
                    _gene_of_cluster_[ *temp_iter ] = module_merge;
                _clusters_.erase( *closest__clusters__iter );//delete the merge the j
            }
            //then cluster the gene into the merged module
            _gene_of_cluster_.insert( map<string,int>::value_type(make_pair( fc_vec_iter->first, module_merge )) );
            _clusters_[ module_merge ].insert( fc_vec_iter -> first );

			//determine largest
			if( _clusters_.find( largest_module) != _clusters_.end() ){
				// the largest module exists, but still the largest?
				if(_clusters_[module_merge].size() > largest_module_size ){
					largest_module = module_merge;
					largest_module_size = _clusters_[module_merge].size();
					vector<string> temp_vec;
					set_intersection( _clusters_[module_merge].begin(), _clusters_[module_merge].end(), seed_genes.begin(), seed_genes.end(), back_inserter(temp_vec) );
					largest_module_seeds = temp_vec.size();
				}
				else if( _clusters_[module_merge].size() == _clusters_[largest_module].size() ){
					second_largest_module = module_merge;
					second_largest_module_size = _clusters_[module_merge].size();
					vector<string> temp_vec;
					set_intersection( _clusters_[module_merge].begin(), _clusters_[module_merge].end(), seed_genes.begin(), seed_genes.end(), back_inserter(temp_vec) );
					second_largest_module_seeds = temp_vec.size();
				}
			}
			else{
				//still the largest, but module number is changed, because merge other modules
				largest_module = module_merge;
				largest_module_size = _clusters_[module_merge].size();
				vector<string> temp_vec;
				set_intersection( _clusters_[module_merge].begin(), _clusters_[module_merge].end(), seed_genes.begin(), seed_genes.end(), back_inserter(temp_vec) );
				largest_module_seeds = temp_vec.size();
			}

			//2nd largest module
			//second_largest_module_size = 0;
			if( largest_module_size != second_largest_module_size){
				int new_second_largest_module = 0;
				int new_second_largest_module_size = 0;
				for(map<int,set<string> >::iterator cluster_iter = _clusters_.begin(); cluster_iter != _clusters_.end(); cluster_iter++){
					if( ( _clusters_[cluster_iter->first].size() > new_second_largest_module_size ) && ( _clusters_[cluster_iter->first].size() < largest_module_size ) ){
						new_second_largest_module_size = _clusters_[cluster_iter->first].size();
						new_second_largest_module = cluster_iter->first;
					}
				}
				second_largest_module_size = new_second_largest_module_size;
				second_largest_module = new_second_largest_module;
				vector<string> temp_vec;
				set_intersection( _clusters_[second_largest_module].begin(), _clusters_[second_largest_module].end(), seed_genes.begin(), seed_genes.end(), back_inserter(temp_vec) );
				second_largest_module_seeds = temp_vec.size();
			}
			

			//record the current statistics; key: clustered genes; value : (1) clustered seeds (2) largest module size (3) largest module seeds (4) 2nd largest module (5) 2nd largest module seeds
			vector<int> curr_stat;
			curr_stat.push_back( clustered_seeds);
			curr_stat.push_back( largest_module_size);
			curr_stat.push_back( largest_module_seeds);
			curr_stat.push_back( second_largest_module_size );
			curr_stat.push_back( second_largest_module_seeds );
			int clustered_gene_num = clustered_genes.size();
			_clustering_process_.insert( map<int,vector<int> >::value_type( make_pair( clustered_gene_num, curr_stat) ) );

        }
    }



    //ofil_quantity << "merge times" << '\t' << "clustered gene counts\t" << "clustered differential gene counts\t"
    //	<< "largest module size\t" << "differential gene counts in the largest module\n";
    omp_set_lock(&lock);
    for( map<int, vector<int> >::iterator temp_iter = _clustering_process_.begin(); temp_iter != _clustering_process_.end(); temp_iter++){
        ofs << simi_thresh << '\t';
        ofs << temp_iter->first << '\t';
        for( int i = 0; i < temp_iter->second.size(); i++)
            ofs << (temp_iter->second)[i] << '\t';
        ofs << endl;
    }
    omp_unset_lock(&lock);

}


void Diff_modules::hma_merge_stop( MatrixXd &simi_mat, Gene_list &gene_list, Parameter &params, PPI &network){
	/*******************
	build time: 2013/08/16
	modification time: 2013/09/15; 2013/9/29
	function: the main function of HMA
	(1) use the two parameters to generate clustering results
    (2) the second part of our application clustex2
	*******************/
	gene_of_cluster.clear();
	clusters.clear();
	/* load the order of simi_mat, same with fold_change */
	set<string> seed_genes = gene_list.get_seeds();
	map<string, double> fold_change = gene_list.get_fold_change();
	int N = fold_change.size();
	map<string, int> simi_mat_order;
	int order = 0;
	for( map<string, double>::iterator fc_iter = fold_change.begin(); fc_iter != fold_change.end(); fc_iter++ ){
		simi_mat_order.insert( map<string,int>::value_type( make_pair( fc_iter->first, order) ) );
		order++;
	}

	//find the top x% similarity threshold
	map<string, set<string> > ppi = network.get_ppi();
	set<double> simi_set;
	for( map<string, set<string> >::iterator ppi_iter = ppi.begin(); ppi_iter != ppi.end(); ppi_iter++){
		//cout << ppi_iter->first ;
		set<string>::iterator set_iter = ppi_iter->second.begin();
		for( ; set_iter != ppi_iter->second.end(); set_iter++){
			//cout <<'\t' << *set_iter;
			simi_set.insert( simi_mat( simi_mat_order[ppi_iter->first] , simi_mat_order[*set_iter] ) );
		}
		//cout << endl;
	}
	//here is the place that can be improved
	vector<double> simi_vec;
	for( set<double>::iterator set_iter = simi_set.begin(); set_iter != simi_set.end(); set_iter++)
		simi_vec.push_back( *set_iter );
	sort( simi_vec.begin(), simi_vec.end() );//larger neighborhood in the back of the vector, namely smaller neighborhood

	double cut_off = params.simi_thresh;//do not know true or not
	double simi_cutoff = simi_vec[ integer( floor( (simi_vec.size() -1)* (1-cut_off) ) ) ];//neighborhood = 1 - similarity threshold
	
	//debug info
	//cout << "The cutoff is " << cut_off << endl;
	//

	/* use HMA method to hierarchically clustering genes */
	vector<string_double> fc_vec;
	for( map<string, double>::iterator fc_iter = fold_change.begin(); fc_iter != fold_change.end(); fc_iter++)
		fc_vec.push_back( make_pair( fc_iter->first, fc_iter->second ) );
	sort( fc_vec.begin(), fc_vec.end(), cmp );
	
	
	vector<string> clustered_genes;
	int largest_module_name;
	double largest_module_size;
	
	// make the first module
	vector<string_double>::iterator fc_vec_iter = fc_vec.begin();
	set<string> temp_1;
	temp_1.insert( fc_vec_iter->first );
	int cluster_num = 1;
	clusters.insert( map<int,set<string> >::value_type(make_pair(cluster_num, temp_1)) );
	gene_of_cluster.insert( map<string,int>::value_type(make_pair(fc_vec_iter->first, cluster_num)) );

	largest_module_name = cluster_num;
	largest_module_size = 1.0;

	clustered_genes.push_back( fc_vec_iter->first );
	fc_vec_iter++;
	cluster_num++;

	
	/* go on to cluster  genes */
	

	for( ; fc_vec_iter != fc_vec.end(); fc_vec_iter++){
		//determine the closest genes and corresponding modules
		vector<string> closest_genes;
		closest_genes.clear();
		set<int>  closest_clusters;
		closest_clusters.clear();
		for( int i = 0; i < clustered_genes.size(); i++){
			if( simi_mat( simi_mat_order[ fc_vec_iter->first ] , simi_mat_order[clustered_genes[i]] ) >= simi_cutoff ){//not less than
				closest_genes.push_back( clustered_genes[i] );
				closest_clusters.insert( gene_of_cluster[ clustered_genes[i] ] );
			}
		}
		if( closest_genes.empty() ){
			clustered_genes.push_back(  fc_vec_iter->first );
			/* the gene generate a new cluster */
			set<string> temp;
			temp.insert( fc_vec_iter->first);
			clusters.insert( map<int,set<string> >::value_type( cluster_num, temp ) );
			gene_of_cluster.insert( map<string,int>::value_type( make_pair(fc_vec_iter->first, cluster_num) ) );
			cluster_num++;
		}
		else if( closest_genes.size() == 1 || closest_clusters.size() == 1 ){
			// only one genes are closest
			clustered_genes.push_back( fc_vec_iter->first );
			(clusters[ gene_of_cluster[ closest_genes[0] ] ] ).insert( fc_vec_iter->first );
			gene_of_cluster.insert( map<string,int>::value_type(fc_vec_iter->first, gene_of_cluster[ closest_genes[0] ]) );
			//determine whether this module is the largest module
			if( clusters[gene_of_cluster[ closest_genes[0] ]].size() > largest_module_size ){
				largest_module_name = gene_of_cluster[ closest_genes[0] ];
				largest_module_size = clusters[ gene_of_cluster[ closest_genes[0] ] ].size();
			}
		}
		else{
			clustered_genes.push_back( fc_vec_iter->first );

			int module_merge = cluster_num;
			//int module_merge = gene_of_cluster[ closest_genes[0] ];
			cluster_num++;
			set<int>::iterator closest_clusters_iter = closest_clusters.begin();
			//here we have to know that after merging, the largest module's size is suddenly enlarged and more than the params.merge_gene_cutoff, after merge, the largest module is larger than 310, then break
			double largest_module_num_after = 0;
			for( ; closest_clusters_iter != closest_clusters.end(); closest_clusters_iter++){
				largest_module_num_after += clusters[*closest_clusters_iter].size();
			}
			//largest_module_num_after += clusters[module_merge].size();
			if( largest_module_num_after >= params.merge_gene_cutoff + 10 ){
				break;
			}
			

			closest_clusters_iter = closest_clusters.begin();
			for( ; closest_clusters_iter != closest_clusters.end(); closest_clusters_iter++){
				(clusters[ module_merge ] ).insert( (clusters[ *closest_clusters_iter ]).begin() , (clusters[ *closest_clusters_iter ]).end() );
				for( set<string>::iterator temp_iter = (clusters[ *closest_clusters_iter ]).begin(); temp_iter != (clusters[ *closest_clusters_iter ]).end(); temp_iter++)
					gene_of_cluster[ *temp_iter ] = module_merge;
				clusters.erase( *closest_clusters_iter );//delete the clustered module
			}
			
			//then cluster the gene into the merged module
			gene_of_cluster.insert( map<string,int>::value_type(make_pair( fc_vec_iter->first, module_merge )) );
			clusters[ module_merge ].insert( fc_vec_iter -> first );

			if( clusters[module_merge].size() > largest_module_size ){
				largest_module_size = clusters[module_merge].size();
				largest_module_name = module_merge;
			}
		}
		
		if( largest_module_size  >= params.merge_gene_cutoff )
			break;
	}
}

void Diff_modules::hma_merge_fix_neighborhood_process(MatrixXd &simi_mat, Gene_list &gene_list, Parameter &params, PPI &network){
	/*******************
	build time: 2013/10/06
	modification time: 2013/10/30;2013/11/15
	function: HMA, with fix neighborhood as 0.1, to cluster all the genes to the end.
	(1) use the two parameters to generate clustering results: neighborhood as 0.1, stop gene no use
	(2) the output file is used to draw the clustering process
	(3) record when cluster ONE gene: clustered genes; clustered seed genes;  largest module size; largest module seed genes; 2nd largest module; seed genes in 2nd largest module
	*******************/
	cout << "Start to show clustering process (-S step)......";
	gene_of_cluster.clear();
	clusters.clear();
	clustering_process.clear();
	/* load the order of simi_mat, same with fold_change */
	set<string> seed_genes = gene_list.get_seeds();
	map<string, double> fold_change = gene_list.get_fold_change();
	int N = fold_change.size();
	map<string, int> simi_mat_order;
	int order = 0;
	for( map<string, double>::iterator fc_iter = fold_change.begin(); fc_iter != fold_change.end(); fc_iter++ ){
		simi_mat_order.insert( map<string,int>::value_type( make_pair( fc_iter->first, order) ) );
		order++;
	}

	//find the top x% similarity threshold
	map<string, set<string> > ppi = network.get_ppi();
	set<double> simi_set;
	for( map<string, set<string> >::iterator ppi_iter = ppi.begin(); ppi_iter != ppi.end(); ppi_iter++){
		//cout << ppi_iter->first ;
		set<string>::iterator set_iter = ppi_iter->second.begin();
		for( ; set_iter != ppi_iter->second.end(); set_iter++){
			//cout <<'\t' << *set_iter;
			simi_set.insert( simi_mat( simi_mat_order[ppi_iter->first] , simi_mat_order[*set_iter] ) );
		}
		//cout << endl;
	}
	//here is the place that can be improved
	vector<double> simi_vec;
	for( set<double>::iterator set_iter = simi_set.begin(); set_iter != simi_set.end(); set_iter++)
		simi_vec.push_back( *set_iter );
	sort( simi_vec.begin(), simi_vec.end() );//larger neighborhood in the back of the vector, namely smaller neighborhood

    double cut_off = 0.1;//the default neighborhood value,, NOTICE now changed to 0.07, need to back to 0.1
	double simi_cutoff = simi_vec[ integer( floor( (simi_vec.size() -1)* (1-cut_off) ) ) ];//neighborhood = 1 - similarity threshold
	//for debug
	//cout << "The current 0.1 neighborhood is " << simi_cutoff << " when -S is on" <<endl;
	//

	vector<string_double> fc_vec;
	for( map<string, double>::iterator fc_iter = fold_change.begin(); fc_iter != fold_change.end(); fc_iter++)
		fc_vec.push_back( make_pair( fc_iter->first, fc_iter->second ) );
	sort( fc_vec.begin(), fc_vec.end(), cmp );


	

	vector<string> clustered_genes;
	// make the first module
	vector<string_double>::iterator fc_vec_iter = fc_vec.begin();
	set<string> temp_1;
	temp_1.insert( fc_vec_iter->first );
	int cluster_num = 1;
	clusters.insert( map<int,set<string> >::value_type(make_pair(cluster_num, temp_1)) );
	gene_of_cluster.insert( map<string,int>::value_type(make_pair(fc_vec_iter->first, cluster_num)) );
	clustered_genes.push_back( fc_vec_iter->first );

	//make the largest module and 2nd largest module
	int largest_module = cluster_num;
	int largest_module_size = 1;
	int largest_module_seeds = 0;
	int second_largest_module = 0;
	int second_largest_module_size = 0;
	int second_largest_module_seeds = 0;
	//clustered seed genes;
	int clustered_seeds = 0;
	if( seed_genes.find(fc_vec_iter->first) != seed_genes.end() )
	{
		largest_module_seeds++;
		clustered_seeds++;
	}
	//

	fc_vec_iter++;
	cluster_num++;

	/* go on to cluster  genes */
	for( ; fc_vec_iter != fc_vec.end(); fc_vec_iter++){
		vector<string> closest_genes;
		set<int>  closest_clusters;
		for( int i = 0; i < clustered_genes.size(); i++){
			// the rule that a gene belongs to a cluster :
			// nearest one of its closest gene
			if( simi_mat( simi_mat_order[ fc_vec_iter->first] , simi_mat_order[clustered_genes[i]] ) >= simi_cutoff ){//not less than
				closest_genes.push_back( clustered_genes[i] );
				closest_clusters.insert( gene_of_cluster[ clustered_genes[i] ] );
			}
		}
		//is this gene a seed gene?
		if( seed_genes.find(fc_vec_iter->first) != seed_genes.end() )
			clustered_seeds++;
		
		clustered_genes.push_back(  fc_vec_iter->first );

		if( closest_genes.empty() ){
			//could not be larger than the largest module
			//clustered_genes.push_back(  fc_vec_iter->first );
			/* the gene generate a new cluster */
			set<string> temp;
			temp.insert( fc_vec_iter->first);
			clusters.insert( map<int,set<string> >::value_type( make_pair(cluster_num, temp)) );
			gene_of_cluster.insert( map<string,int>::value_type( make_pair(fc_vec_iter->first, cluster_num) ) );
			if( clusters[cluster_num].size() > second_largest_module_size  ){
				second_largest_module = cluster_num;
				second_largest_module_size = 1;
				if( seed_genes.find( fc_vec_iter->first) != seed_genes.end() )
					second_largest_module_seeds = 1;
			}
			//record the current statistics; key: clustered genes; value : (1) clustered seeds (2) largest module size (3) largest module seeds (4) 2nd largest module (5) 2nd largest module seeds
			vector<int> curr_stat;
			curr_stat.push_back( clustered_seeds);
			curr_stat.push_back( largest_module_size);
			curr_stat.push_back( largest_module_seeds);
			curr_stat.push_back( second_largest_module_size );
			curr_stat.push_back( second_largest_module_seeds );
			int clustered_gene_num = clustered_genes.size();
			clustering_process.insert( map<int,vector<int> >::value_type( make_pair( clustered_gene_num, curr_stat) ) );
			
			cluster_num++;
		}
		else if( closest_genes.size() == 1 || closest_clusters.size() == 1 ){
			//clustered_genes.push_back( fc_vec_iter->first );
			/* only one genes are closest */ 
			( clusters[ gene_of_cluster[ closest_genes[0] ] ] ).insert( fc_vec_iter->first );
			gene_of_cluster.insert( map<string,int>::value_type(make_pair(fc_vec_iter->first, gene_of_cluster[ closest_genes[0] ]) ) );
			
			//determine largest module and second largest module
			//notice: whether the second largest module is as large as the largest module
			int temp_size = clusters[ gene_of_cluster[ closest_genes[0] ] ].size();
			//cout << "temp_size: " << temp_size << endl;
			if( temp_size > largest_module_size ){
				largest_module = gene_of_cluster[ closest_genes[0] ];
				largest_module_size = clusters[ gene_of_cluster[ closest_genes[0] ] ].size();
				vector<string> temp_vec;
				set_intersection( clusters[largest_module].begin(), clusters[largest_module].end(), seed_genes.begin(), seed_genes.end(), back_inserter(temp_vec) );
				largest_module_seeds = temp_vec.size();

				//second_largest_module_size = 0;
				int new_second_largest_module_size = 0;
				int new_second_largest_module = 0;
				for(map<int ,set<string> >::iterator clusters_iter = clusters.begin(); clusters_iter != clusters.end(); clusters_iter++ ){
					if(  (clusters[clusters_iter->first].size() < largest_module_size) && ( clusters[clusters_iter->first].size() > new_second_largest_module_size ) ){
						new_second_largest_module_size = clusters[clusters_iter->first].size();
						new_second_largest_module = clusters_iter->first;
					}
				}
				second_largest_module = new_second_largest_module;
				second_largest_module_size = new_second_largest_module_size;
				vector<string> temp_vec_2;
				set_intersection( clusters[second_largest_module].begin(), clusters[second_largest_module].end(), seed_genes.begin(), seed_genes.end(), back_inserter(temp_vec_2) );
				second_largest_module_seeds = temp_vec_2.size();

			}
			else if( (temp_size < largest_module_size) && (temp_size > second_largest_module_size ) ){
				if( second_largest_module != gene_of_cluster[ closest_genes[0] ] ){
					second_largest_module = gene_of_cluster[ closest_genes[0] ];
					second_largest_module_size = clusters[ gene_of_cluster[ closest_genes[0] ] ].size();
					vector<string> temp_vec;
					set_intersection( clusters[ gene_of_cluster[ closest_genes[0] ] ].begin(), clusters[ gene_of_cluster[ closest_genes[0] ] ].end(), seed_genes.begin(), seed_genes.end(), back_inserter(temp_vec) );
					second_largest_module_seeds = temp_vec.size();
				}
				else if ( second_largest_module == gene_of_cluster[closest_genes[0]] ){
					second_largest_module_size++;
					if( seed_genes.find(fc_vec_iter->first) != seed_genes.end() )
						second_largest_module_seeds++;
				}
			}
			else if( temp_size == largest_module_size  ){
				second_largest_module = gene_of_cluster[ closest_genes[0] ];
				second_largest_module_size = clusters[ gene_of_cluster[ closest_genes[0] ] ].size();
				vector<string> temp_vec;
				set_intersection( clusters[ gene_of_cluster[ closest_genes[0] ] ].begin(), clusters[ gene_of_cluster[ closest_genes[0] ] ].end(), seed_genes.begin(), seed_genes.end(), back_inserter(temp_vec) );
				second_largest_module_seeds = temp_vec.size();
			}
			//record the current statistics; key: clustered genes; value : (1) clustered seeds (2) largest module size (3) largest module seeds (4) 2nd largest module (5) 2nd largest module seeds
			vector<int> curr_stat;
			curr_stat.push_back( clustered_seeds);
			curr_stat.push_back( largest_module_size);
			curr_stat.push_back( largest_module_seeds);
			curr_stat.push_back( second_largest_module_size );
			curr_stat.push_back( second_largest_module_seeds );
			int clustered_gene_num = clustered_genes.size();
			clustering_process.insert( map<int,vector<int> >::value_type( make_pair( clustered_gene_num, curr_stat) ) );
		}
		else{
			//clustered_genes.push_back( fc_vec_iter->first );
			/* this gene is closed to several clusters */
			//first merge these clusters
			//int module_merge = gene_of_cluster[ closest_genes[0] ];
			int module_merge = cluster_num;
			cluster_num++;
			set<int>::iterator closest_clusters_iter = closest_clusters.begin();
			//start to merge
			for( ; closest_clusters_iter != closest_clusters.end(); closest_clusters_iter++){
				//merge
				( clusters[ module_merge ] ).insert( ( clusters[ *closest_clusters_iter ]).begin() , ( clusters[ *closest_clusters_iter ]).end() );
				for( set<string>::iterator temp_iter = ( clusters[ *closest_clusters_iter ]).begin(); temp_iter != ( clusters[ *closest_clusters_iter ]).end(); temp_iter++)
					gene_of_cluster[ *temp_iter ] = module_merge;
				clusters.erase( *closest_clusters_iter );//delete the merged module
			}
			//then cluster the gene into the merged module
			gene_of_cluster.insert( map<string,int>::value_type(make_pair( fc_vec_iter->first, module_merge )) );
			clusters[ module_merge ].insert( fc_vec_iter -> first );

			
			
			//determine largest
			if( clusters.find( largest_module) != clusters.end() ){// the largest module exists, but still the largest?
				if( clusters[module_merge].size() > clusters[largest_module].size() ){//former: clusters[module_merge].size() > largest_module_size
					largest_module = module_merge;
					largest_module_size = clusters[module_merge].size();
					vector<string> temp_vec;
					set_intersection( clusters[module_merge].begin(), clusters[module_merge].end(), seed_genes.begin(), seed_genes.end(), back_inserter(temp_vec) );
					largest_module_seeds = temp_vec.size();
				}
				else if( clusters[module_merge].size() == clusters[largest_module].size() ){
					second_largest_module == module_merge;
					second_largest_module_size = clusters[module_merge].size();
					vector<string> temp_vec;
					set_intersection( clusters[module_merge].begin(), clusters[module_merge].end(), seed_genes.begin(), seed_genes.end(), back_inserter(temp_vec) );
					second_largest_module_seeds = temp_vec.size();
				}// the condition: clusters[module_merge].size() == clusters[largest_module].size() : detemine second largest
			}
			else{//still the largest, but module number is changed, because merge other modules
				largest_module = module_merge;
				largest_module_size = clusters[module_merge].size();
				vector<string> temp_vec;
				set_intersection( clusters[module_merge].begin(), clusters[module_merge].end(), seed_genes.begin(), seed_genes.end(), back_inserter(temp_vec) );
				largest_module_seeds = temp_vec.size();
			}

			//2nd largest module
			//second_largest_module_size = 0
			if( largest_module_size != second_largest_module_size ){
				int new_second_largest_module_size = 0;
				int new_second_largest_module = 0;
				for(map<int,set<string> >::iterator cluster_iter = clusters.begin(); cluster_iter != clusters.end(); cluster_iter++){
					if( ( clusters[cluster_iter->first].size() > new_second_largest_module_size ) && ( clusters[cluster_iter->first].size() < largest_module_size ) ){
						new_second_largest_module_size = clusters[cluster_iter->first].size();
						new_second_largest_module = cluster_iter->first;
					}
				}
				second_largest_module_size = new_second_largest_module_size;
				second_largest_module = new_second_largest_module;
				vector<string> temp_vec;
				set_intersection( clusters[second_largest_module].begin(), clusters[second_largest_module].end(), seed_genes.begin(), seed_genes.end(), back_inserter(temp_vec) );
				second_largest_module_seeds = temp_vec.size();
			}

			//record the current statistics; key: clustered genes; value : (1) clustered seeds (2) largest module size (3) largest module seeds (4) 2nd largest module (5) 2nd largest module seeds
			vector<int> curr_stat;
			curr_stat.push_back( clustered_seeds);
			curr_stat.push_back( largest_module_size);
			curr_stat.push_back( largest_module_seeds);
			curr_stat.push_back( second_largest_module_size );
			curr_stat.push_back( second_largest_module_seeds );
			int clustered_gene_num = clustered_genes.size();
			clustering_process.insert( map<int,vector<int> >::value_type( make_pair( clustered_gene_num, curr_stat) ) );
		}
	}


	string str = params.job_id + "_fix_neighborhood_0.1.txt";
	ofstream ofs( str.c_str() );
	ofs << "clustered genes\t" << "clustered seed genes\t" 
		   << "largest module size\t" << "seed genes in largest\t" << "2nd largest module size\t" << "seed genes in 2nd largest\n";
	for( map<int, vector<int> >::iterator temp_iter = clustering_process.begin(); temp_iter != clustering_process.end(); temp_iter++){
		ofs << temp_iter->first << '\t';
		for( int i = 0; i < temp_iter->second.size(); i++)
			ofs << (temp_iter->second)[i] << '\t';
		ofs << endl;
	}
	cout << "done!\n";
}

void Diff_modules::hma_merge_fix_neighborhood(MatrixXd &simi_mat, Gene_list &gene_list, Parameter &params, PPI &network){
	/*******************
	build time: 2013/10/06
	modification time: 2013/10/30
	function: the main function of HMA, with fixed neighborhood 0.1
	(1) use the two parameters to generate clustering results, neighborhood as 0.1 and stop genes with user input
    (2) the second part of our application clustex2
	(3) record the largest module and second module
	*******************/
	gene_of_cluster.clear();
	clusters.clear();
	/* load the order of simi_mat, same with fold_change */
	set<string> seed_genes = gene_list.get_seeds();
	map<string, double> fold_change = gene_list.get_fold_change();
	int N = fold_change.size();
	map<string, int> simi_mat_order;
	int order = 0;
	for( map<string, double>::iterator fc_iter = fold_change.begin(); fc_iter != fold_change.end(); fc_iter++ ){
		simi_mat_order.insert( map<string,int>::value_type( make_pair( fc_iter->first, order) ) );
		order++;
	}

	//find the top x% similarity threshold
	map<string, set<string> > ppi = network.get_ppi();
	set<double> simi_set;
	for( map<string, set<string> >::iterator ppi_iter = ppi.begin(); ppi_iter != ppi.end(); ppi_iter++){
		set<string>::iterator set_iter = ppi_iter->second.begin();
		for( ; set_iter != ppi_iter->second.end(); set_iter++){
			simi_set.insert( simi_mat( simi_mat_order[ppi_iter->first] , simi_mat_order[*set_iter] ) );
		}
	}
	//here is the place that can be improved
	vector<double> simi_vec;
	for( set<double>::iterator set_iter = simi_set.begin(); set_iter != simi_set.end(); set_iter++)
		simi_vec.push_back( *set_iter );
	sort( simi_vec.begin(), simi_vec.end() );//larger neighborhood in the back of the vector, namely smaller neighborhood

	double cut_off = 0.1;//do not know true or not
	double simi_cutoff = simi_vec[ integer( floor( (simi_vec.size() -1)* (1-cut_off) ) ) ];//neighborhood = 1 - similarity threshold
	//for debug
	//cout << "The current 0.1 neighborhood is " << simi_cutoff << " when -C is on" <<endl;
	//


	/* use HMA method to hierarchically clustering genes */
	vector<string_double> fc_vec;
	for( map<string, double>::iterator fc_iter = fold_change.begin(); fc_iter != fold_change.end(); fc_iter++)
		fc_vec.push_back( make_pair( fc_iter->first, fc_iter->second ) );
	sort( fc_vec.begin(), fc_vec.end(), cmp );
	
	
	vector<string> clustered_genes;
	int largest_module_name;
	double largest_module_size;
	
	// make the first module
	vector<string_double>::iterator fc_vec_iter = fc_vec.begin();
	set<string> temp_1;
	temp_1.insert( fc_vec_iter->first );
	int cluster_num = 1;
	clusters.insert( map<int,set<string> >::value_type(make_pair(cluster_num, temp_1)) );
	gene_of_cluster.insert( map<string,int>::value_type(make_pair(fc_vec_iter->first, cluster_num)) );

	largest_module_name = cluster_num;
	largest_module_size = 1.0;

	clustered_genes.push_back( fc_vec_iter->first );
	fc_vec_iter++;
	cluster_num++;

	
	/* go on to cluster  genes */

	for( ; fc_vec_iter != fc_vec.end(); fc_vec_iter++){
		//determine the closest genes and corresponding modules
		vector<string> closest_genes;
		closest_genes.clear();
		set<int>  closest_clusters;
		closest_clusters.clear();
		for( int i = 0; i < clustered_genes.size(); i++){
			if( simi_mat( simi_mat_order[ fc_vec_iter->first ] , simi_mat_order[clustered_genes[i]] ) >= simi_cutoff ){//not less than
				closest_genes.push_back( clustered_genes[i] );
				closest_clusters.insert( gene_of_cluster[ clustered_genes[i] ] );
			}
		}
		if( closest_genes.empty() ){
			clustered_genes.push_back(  fc_vec_iter->first );
			/* the gene generate a new cluster */
			set<string> temp;
			temp.insert( fc_vec_iter->first);
			clusters.insert( map<int,set<string> >::value_type( cluster_num, temp ) );
			gene_of_cluster.insert( map<string,int>::value_type( make_pair(fc_vec_iter->first, cluster_num) ) );
			cluster_num++;
		}
		else if( closest_genes.size() == 1 || closest_clusters.size() == 1 ){
			// only one genes are closest
			clustered_genes.push_back( fc_vec_iter->first );
			(clusters[ gene_of_cluster[ closest_genes[0] ] ] ).insert( fc_vec_iter->first );
			gene_of_cluster.insert( map<string,int>::value_type(fc_vec_iter->first, gene_of_cluster[ closest_genes[0] ]) );
			//determine whether this module is the largest module
			if( clusters[gene_of_cluster[ closest_genes[0] ]].size() > largest_module_size ){
				largest_module_name = gene_of_cluster[ closest_genes[0] ];
				largest_module_size = clusters[ gene_of_cluster[ closest_genes[0] ] ].size();
			}
		}
		else{
			clustered_genes.push_back( fc_vec_iter->first );

			int module_merge = cluster_num;
			//int module_merge = gene_of_cluster[ closest_genes[0] ];
			cluster_num++;
			set<int>::iterator closest_clusters_iter = closest_clusters.begin();
			//here we have to know that after merging, the largest module's size is suddenly enlarged and more than the params.merge_gene_cutoff, after merge, the largest module is larger than 310, then break
			double largest_module_num_after = 0;
			for( ; closest_clusters_iter != closest_clusters.end(); closest_clusters_iter++){
				largest_module_num_after += clusters[*closest_clusters_iter].size();
			}
			//largest_module_num_after += clusters[module_merge].size();
			if( largest_module_num_after >= params.merge_gene_cutoff + 10 ){
				break;
			}
			//start to merge
			//for debug
			//ofil << "Merging module " << module_merge << "\tINFO: " << clusters[module_merge].size() << "\t";
			//

			closest_clusters_iter = closest_clusters.begin();
			for( ; closest_clusters_iter != closest_clusters.end(); closest_clusters_iter++){
				(clusters[ module_merge ] ).insert( (clusters[ *closest_clusters_iter ]).begin() , (clusters[ *closest_clusters_iter ]).end() );
				for( set<string>::iterator temp_iter = (clusters[ *closest_clusters_iter ]).begin(); temp_iter != (clusters[ *closest_clusters_iter ]).end(); temp_iter++)
					gene_of_cluster[ *temp_iter ] = module_merge;
				clusters.erase( *closest_clusters_iter );//delete the clustered module
			}
			
			//then cluster the gene into the merged module
			gene_of_cluster.insert( map<string,int>::value_type(make_pair( fc_vec_iter->first, module_merge )) );
			clusters[ module_merge ].insert( fc_vec_iter -> first );

			if( clusters[module_merge].size() > largest_module_size ){
				largest_module_size = clusters[module_merge].size();
				largest_module_name = module_merge;
			}
		}
		
		if( largest_module_size  >= params.merge_gene_cutoff )
			break;
	}
}


void Diff_modules::print_clusters(Gene_list &gene_list, PPI &network, Parameter &params){
	/***************
	build time: 2013/08/16
	function: print out the clustering results, including
	(1) the differential modules, ignore modules with only 1 gene; format: gene 1 TAB gene 2 TAB differetial module number(0 represents edges between two module); visualized by Cytoscape
	(2) clusters; format: differential module number TAB gene 1 .... n
	(3) gene in clusters; format: gene a TAB module number TAB seed gene or (1 == seed gene 0 != seed gene)
	******make sure that the output is consensus
	now some modules with only one genes is also output, which should not out output anymore
	***************/
	//(0) remove all the single gene module
	map<int,set<string> >::iterator cluster_iter_remove = clusters.begin();
	int cluster_num = 0;
	while( cluster_iter_remove != clusters.end() ){
		if( cluster_iter_remove->second.size() == 0 ){
			cerr << "The clustering results is WRONG! Module " << cluster_iter_remove->first << " size is 0!\n";
			exit(-1);
		}
		if( cluster_iter_remove->second.size() == 1){// this place could be set as a input parameter, default is 1 
			set<string>::iterator temp_iter = cluster_iter_remove->second.begin();
			gene_of_cluster.erase( *temp_iter );
			map<int,set<string> >::iterator cluster_iter_temp = cluster_iter_remove;
			cluster_iter_remove++;
			clusters.erase( cluster_iter_temp );
		}
		else{
			cluster_num++;
			cluster_iter_remove++;
		}
	}
	cout << "There are "  << cluster_num << " gene modules in total\n";
	//(1) the differential modules
	map<string, set<string> > ppi = network.get_ppi();
	set<string> seed_genes = gene_list.get_seeds();
	//set<int> single_gene_clusters;
	map<int , set<string> >::iterator clusters_iter = clusters.begin();

	string str_1 = "../../temp/" + params.job_id + "_modules_"; // testing-onfah
	stringstream sstrm; 
	string stop_gene;
	string neighborhood;
	sstrm << params.simi_thresh;
	sstrm >> neighborhood;
	sstrm.clear();
	sstrm << params.merge_gene_cutoff;
	sstrm >> stop_gene;
	string xiahuaxian = "_";
	string str;
	str = str_1 + stop_gene +xiahuaxian+ neighborhood+ ".txt" ;
	ofstream ofil_3( str.c_str() );
	if( !ofil_3 ){
		cerr << "Cannot open the file to output the gene modules!\n";
		exit(-1);
	}
	ofil_3 << "GeneA\tGeneB\tModule_num\n";
	for( clusters_iter = clusters.begin(); clusters_iter != clusters.end(); clusters_iter++){
		//if( clusters_iter->second.size() == 1){
		//	single_gene_clusters.insert( clusters_iter->first );
		//	continue;
		//}
		for( set<string>::iterator node_iter_1 = clusters_iter->second.begin(); node_iter_1 != clusters_iter->second.end(); node_iter_1++ ){
			set<string>::iterator node_iter_2 = node_iter_1;
			node_iter_2++;
			if( node_iter_2 == clusters_iter->second.end() )
				break;
			for( ; node_iter_2 != clusters_iter->second.end(); node_iter_2++)
				if( ( ppi[*node_iter_1].find( *node_iter_2) != ppi[*node_iter_1].end() ) && ( ppi[*node_iter_2].find( *node_iter_1) != ppi[*node_iter_2].end() ) )
					ofil_3 <<  *node_iter_1 << '\t' << *node_iter_2 << '\t' << clusters_iter->first << "\n";
		}
	}
	/*output the edges between modules
	for( clusters_iter = clusters.begin(); clusters_iter != clusters.end(); clusters_iter++){
		map<int, set<string> >::iterator clusters_iter_2 = clusters_iter;
		clusters_iter_2++;
		if( clusters_iter_2 == clusters.end() )
			break;
		else if ( clusters_iter_2->second.size() == 1 )
			continue;
		else{
			for( set<string>::iterator node_iter_1 = clusters_iter->second.begin(); node_iter_1 != clusters_iter->second.end(); node_iter_1++){
				for( set<string>::iterator node_iter_2 = clusters_iter_2->second.begin(); node_iter_2 != clusters_iter_2->second.end(); node_iter_2++){
					if( ppi[*node_iter_2].find( *node_iter_1) != ppi[*node_iter_2].end() )
						ofil_3 << *node_iter_1 << "\t" << *node_iter_2 << "0\n";
				}
			}
		}
	}
	*/
	ofil_3.close();

	//(2) the clusters
	str_1.clear();
	str_1 = "../../temp/" + params.job_id + "_clusters_"; // testing-onfah
	str.clear();
	str = str_1 + stop_gene +xiahuaxian+ neighborhood+ ".txt" ;
	ofstream ofil_1( str.c_str() );
	if( !ofil_1 ){
		cerr << "Cannot open the file to output modules!\n";
		exit(-1);
	}
	
	ofil_1 << "The parameters used to generate the following results are : \n"
		        << "largest module size: " << stop_gene << endl
				<< "neighborhood: " << neighborhood << endl
				<< "restart_prob: " << params.restart_prob << endl
				<< "diffusion kernel coefficient: " << params.diff_beta << endl
				<< "Modules: " << endl;
	for( clusters_iter = clusters.begin(); clusters_iter != clusters.end(); clusters_iter++){
		set<string>::iterator temp_iter = ( clusters_iter->second ).begin();
		ofil_1 << "module "  << clusters_iter->first  << "\t";
		int seed_num = 0;
		for( set<string>::iterator itt = (clusters_iter->second).begin(); itt != (clusters_iter->second).end(); itt++){
			if( seed_genes.find(*itt) != seed_genes.end() ){
				seed_num++;
			}
		}
		ofil_1 << "seed_gene_fraction " << ((double) seed_num )/ ( (double) (clusters[clusters_iter->first]).size() ) << "\t" ;
		ofil_1 << "number of genes: " << clusters[clusters_iter->first].size() << "\t";
		ofil_1 << "genes:";
		for( ; temp_iter != ( clusters_iter->second).end(); temp_iter++ ){
			ofil_1 << "\t" << *temp_iter;
		}
		ofil_1 << endl;
		
	}

	//(3) gene in clusters 
	str.clear();
	str_1.clear();
	str_1 = "../../temp/" + params.job_id +  "_genes_"; // testing-onfah
	str = str_1 + stop_gene +xiahuaxian+ neighborhood + ".txt" ;
	ofstream ofil_2( str.c_str() );
	if( !ofil_2 ){
		cerr << "Cannot open the file to output genes with corresponding modules!\n";
		exit(-1);
	}
	// ofil_2 << "Gene\tModule\tseed_gene\n"; // testing-onfah
	map<string, int >::iterator gene_iter = gene_of_cluster.begin();
	for( ; gene_iter != gene_of_cluster.end(); gene_iter++){
		ofil_2 << gene_iter->first << "\n"; // testing-onfah
		/*
		ofil_2 << gene_iter->first
			<< '\t'
			<< gene_iter->second
			<<'\t';
		//whether this gene is a seed gene
		if( seed_genes.find( gene_iter->first) != seed_genes.end() )
			ofil_2 << "1\n";
		else
			ofil_2 << "0\n";
		*/
	}
	ofil_2.close();
	ofil_1.close();
}

void Diff_modules::clustering_results_check(){
	/******************
	build time: 2013/09/26
	function:  check
	(1) consensus between gene_of_cluster, clusters
	(2) delete the single gene module
	*******************/
	//(1) consensus between gene_of_cluster, clusters
	ofstream ofil("clustering_consensus_error.txt");
	for(map<string, int>::iterator gene_of_cluster_iter = gene_of_cluster.begin(); gene_of_cluster_iter != gene_of_cluster.end(); gene_of_cluster_iter++){
		set<string>::iterator temp_iter = clusters[gene_of_cluster_iter->second].find( gene_of_cluster_iter->first );
		if( temp_iter == clusters[gene_of_cluster_iter->second].end() ){
			ofil << "module: " <<  gene_of_cluster_iter->second << "\t" << "Gene: " << gene_of_cluster_iter->first << "\t";
			ofil << "Info: " << "module size\t"  << clusters[gene_of_cluster_iter->second].size() << "\n";
		}
	}
	//(2)
	
}

Parameter::Parameter(){
	//parameter
	diff_beta = 0.01;
	restart_prob = 0.5;
	seed_gene_thresh = 0;
	simi_thresh = 0.1;
	simi_thresh_input = false;
	//merge_time_cutoff = 0;
	merge_gene_cutoff = 0;
	parameter_level = 1;
	//file manipulation
	network_file_input = false;
	gene_file_input = false;
	diff_kern_file_input = false;
	density_file_input = false;
	//step choice
	diff_kern_step = false;
	gene_density_step = false;
	parameter_calc_step = false;
	show_clustering_process = false;
	clustering_step = false;
	//function
	// out_middle_results = true;
	out_middle_results = false; // testing-onfah
	//job id
	time_t t;
	time(&t);
	char *temp;
	temp = asctime( localtime(&t) );
	job_id = temp;
	//job_id[job_id.size() -1 ] ='\0';
	job_id_input = false;
	replace_all( job_id, "\n","_");
	replace_all( job_id, " ", "_");
}



Parameter parameter_manipulate(int argc, char *argv[]){
	/**************************
	build time: 2013/08/17
	modification time: 2013/08/23
	function: manipulate all the parameters
	requirement:
	(1) four parts: four steps
	(2) principle: 
	       1) if read any of the two middle results: gene density file and diffusion kernel file, than do not need to calculate them
		   2) suggest the users to run in two ways :
		          a) first run the firs three parts, then the last part
				  b) if changed the gene list but not the network, the first principle works(diffusion kernel file exists)
	***************************/
	Parameter params;

	size_t i;
	string op;
	stringstream sstrm;
	if( argc == 1){
		help();
		exit(-1);
	}
	for( i =1; i < argc; i++){
		op = argv[i];
		if( op == "-h" || op == "--help"){
			help();
			exit(0);
		}
		//step choice
		else if( op == "-D" || op == "--diffusion_kernel_calc"){
			params.diff_kern_step = true;
		}
		else if( op == "-G" || op == "--gene_importance_calc"){
			params.gene_density_step = true;
		}
        else if( op == "-P" || op == "--parameter_greedy_search"){
			params.parameter_calc_step = true;
			i++;
			sstrm.clear();
			if( i == argc ){
				printf("Lack value of %s\n" ,argv[i-1]);
				help();
				exit(-1);
			}
			sstrm << argv[i];
			sstrm >> params.parameter_level;
			if( params.parameter_level != 1 && params.parameter_level != 2 && params.parameter_level !=3 && params.parameter_level !=4 ){
				printf("Invalid value of %s", argv[i-1]);
				help();
				exit(-1);
			}
		}
		//else if( op == "-S" || op == "--show_clustering_process" ){ 
		//	params.show_clustering_process = true;
		//}
		else if( op == "-C" || op == "--clustering_calc" ){
			params.clustering_step = true;
		}
		//job id
		else if( op == "-j" || op == "--job"){
			i++;
			if( i== argc ){
				printf("Lack value of %s\n",argv[i-1]);
				help();
				exit(-1);
			}
			params.job_id_input = true;
			params.job_id = argv[i];
		}
		//threshold choice
		else if( op == "-r" || op == "--restart_prob"){
			sstrm.clear();
			i++;
			if( i == argc ){
				printf("Lack value of %s\n", argv[i-1]);
				help();
				exit(-1);
			}
			sstrm << argv[i];
			sstrm >> params.restart_prob;
			if( !sstrm || params.restart_prob <= 0 || params.restart_prob > 1.0 ){
				printf("Invalid restart probability value %s \n", argv[i]);
				help();
				exit(-1);
			}
		}
		else if( op == "-w" || op == "--weight_threshold"){
			sstrm.clear();
			i++;
			if(i == argc ){
				printf("Lack value of %s", argv[i-1]);
				help();
				exit(-1);
			}
			sstrm << argv[i];
			sstrm >> params.seed_gene_thresh;
			if( !sstrm || params.seed_gene_thresh <= 0){
				printf("Invalid seed gene threshold value %s \n", argv[i]); // in the random walk should talk about it again
				help();
				exit(-1);
			}
		}
		else if( op == "-d" || op == "--diffusion_beta"){
			sstrm.clear();
			i++;
			if( i == argc ){
				printf("Lack value of %s\n", argv[i-1]);
				help();
				exit(-1);
			}
			sstrm <<  argv[i];
			sstrm >> params.diff_beta;
			if( !sstrm || params.diff_beta <= 0 ){ //diffusion kernel beta should be larger than 0 right?
				printf("Invalid diffusion kernel coefficient %s\n", argv[i]);
				help();
				exit(-1);
			}
		}
		else if( op == "-n" || op == "--neighborhood"){
			sstrm.clear();
			i++;
			if( i == argc ){
				printf("Lack value of %s\n", argv[i-1] );
				help();
				exit(-1);
			} 
			sstrm << argv[i];
			sstrm >> params.simi_thresh;
			params.simi_thresh_input = true;
			//params.simi_input = true;
			if( !sstrm || params.simi_thresh <= 0 || params.simi_thresh >= 1 ){
				printf("Invalid neighborhood %s", argv[i]);
				help();
				exit(-1);
			}
		}
		else if( op == "-s" || op == "--size"){
			sstrm.clear();
			i++;
			if( i == argc ){
				printf( "Lack value of %s\n", argv[i-1] );
				help();
				exit(-1);
			}
			sstrm << argv[i];
			sstrm >> params.merge_gene_cutoff;
			//params.merge_cutoff_input = true;
			if( !sstrm || params.merge_gene_cutoff < 0 ){
				printf("Invalid stop time %s", argv[i]);
				help();
				exit(-1);
			}
		}
		else if( op == "--output_middle_results" ){
			params.out_middle_results = true;
		}
		else if( op == "--gene_list" ){
			i++;
			if( i == argc ){
				printf("Lack value of %s\n", argv[i-1]);
				help();
				exit(-1);
			}
			params.gene_file = argv[i];
			params.gene_file_input = true;
		}
		else if( op == "--network"){
			i++;
			if( i == argc ){
				printf("Lack value of %s\n", argv[i-1]);
				help();
				exit(-1);
			}
			params.network_file = argv[i];
			params.network_file_input =true;
		}
		else if( op == "--diffusion_kernel"){
			params.diff_kern_file_input = true;
			i++;
			sstrm.clear();
			if( i == argc ){
				printf( "Lack value of %s", argv[i-1]);
				help();
				exit(-1);
			}
			sstrm << argv[i];
			sstrm >> params.diff_kern_file;
		}
        else if( op == "--gene_importance" ){
			params.density_file_input = true;
			i++;
			sstrm.clear();
			if( i == argc ){
				printf("Invalid value of %s",argv[i-1]);
				help();
				exit(-1);
			}
			sstrm << argv[i];
			sstrm >> params.density_file;
        }
        else{
            cout << "Invalid option: " << argv[i] << "!\n";
            help();
            exit(-1);
        }
	}
	//check out the parameters
	if( params.network_file_input == false || params.gene_file_input == false ){
		cout <<  "Lack of input file: gene list file and gene network file!\n";
		help();
		exit(-1);
	}
	if( params.density_file_input == true )
		params.gene_density_step = false;
	if( params.diff_kern_file_input == true )
		params.diff_kern_step = false;
	/*
	if( params.job_id_input == false ){
		cerr << "Lack of job id!\n";
		help();
		exit(-1);
	}
	*/
	
	if( params.clustering_step == true ){
		if( params.diff_kern_step == false && params.diff_kern_file_input == false){
			cerr << "Diffusion kernel file does not load successfully, please check your command line!\n";
			help();
			exit(-1);
		}
		if( params.gene_density_step == false && params.density_file_input == false ){
			cerr << "Gene importance file does not load successfully, please check your command line!\n";
			help();
			exit(-1);
		}
	}
	
	if( params.density_file_input == false && params.gene_density_step == false  ){
		cerr << "There is neither gene importance calculation nor gene importance input, please check your command line!\n";
		help();
		exit(-1);
	}
	if( params.diff_kern_file_input == false && params.diff_kern_step ==false ){
		cerr << "There is neither diffusion kernel calculation nor diffusion kernel input, please check your command line!\n";
		help();
		exit(-1);
	}

	return params;
}



void help(){
	cout << SOFTWARE << " " << VERSION <<endl << endl;
	cout << "Usage:" <<endl;
    cout << "-h, --help                                                                                                    print this help info" << endl<<endl;
	cout << "<steps choice>"<<endl;
    cout << "-D, --diffusion_kernel_calc                                                                                   calculate the diffusion kernel of the network"<< endl;
    cout << "-G, --gene_importance _calc                                                                                   calculate genes importance" << endl;
    cout << "-P, --parameter_greedy_search                          [int]                                                  suggests inputs for clustering in three levels namely 1, 2, 3 and 4" << endl;
    //cout << "-S, --show_clustering_process                                                                               show clustering process with specific neighborhood 0.1" << endl;
    cout << "-C, --clustering_calc                                                                                         cluster genes and output modules" << endl << endl;
    cout << "<adjustable parameters>"<<endl;
    cout << "-r, --restart_prob                                     [double]                                               set restart probability for random walk,from 0 to 1,default 0.5"<< endl;
    cout << "-d, --diffusion_beta                                   [double]                                               diffusion kernel parameter, default is 0.01" << endl;
    cout << "-n, --neighborhood                                     [double]                                               set neighborhood when -C is on, default is 0.1" << endl;
    cout << "<threshold setting>" << endl;
    cout << "-w, --weight_threshold                                 [double]                                               set a threshold to determine seed genes when given initial gene weights" << endl;
    cout << "-s, --size                                             [int]                                                  clustering stops by controlling  the size of  largest module, when -C is on" << endl << endl;
    cout << "<basic Input/Output>"<<endl;
    cout << "-j, --job                                              [string]                                               job name" << endl;
    cout << "--diffusion_kernel                                     [string]                                               diffusion kernel file, one of ClustEx 2.0's middle results, input when -P or -C is on" << endl;
    cout << "--gene_importance                                      [string]                                               gene importance score file, one of ClustEx 2.0 middle results, input along with --diffusion_kernel" << endl;
    cout << "--gene_list                                            [string]                                               input the candidate gene list file" << endl;
    cout << "--network                                              [string]                                               input the gene network file"<< endl << endl;

}

