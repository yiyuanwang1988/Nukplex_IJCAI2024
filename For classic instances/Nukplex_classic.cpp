#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <libgen.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <algorithm>
#include <sstream>
#include <list>
#include <utility>
#include "utils.hpp"
using namespace std;
#define LARGE_INT 2147483647
/*Parameters*/
int* time_stamp;
char param_graph_file_name[1024]="/home/chen/benchmarks/splex/2nd_dimacs/brock400_4.clq";
int param_s = 2;
int param_best = 9999;
int param_max_seconds = 120;
int param_cycle_iter = 4000;
unsigned int param_seed;
int total_add = 0;
int sai_add = 0;
int count_add = 0;
int total_swap = 0;
int sai_swap = 0;
int count_swap = 0;
/*original graph, the structure is kept the same,
 * but the id of vertices are renumbered, the original
 * id can be retrievaled by org_vid*/
int aaa;
int org_vnum;
int org_enum;
int org_fmt;
int* org_v_edge_cnt;
int** org_v_adj_vertex;
int *org_vid;	/*The real id of each original vertex*/
int count_per = 0;
int size_per = 0;
int** no_org_v_adj_vertex;		//no adjacent matrix data
/*reduced graph*/
int red_vnum;
int red_enum;
int *red_orgid; /*the corresponding id in original graph*/
int* red_v_edge_cnt;
int** red_v_adj_vertex;
int red_min_deg ;
int *freq;
int *freq1;
int *momentum;
int restart_pass;
int *threshold;
int *deposit;
int cur_iter;
int *cur_c_deg;
int *cur_c_consat; /*cur_c_consat[v] Then number of saturated neighboors of v*/
int cur_c_satu_num; //the number of saturate vertices
int *is_in_c;
int *plex_missing;		//plex_missing[v]: the number of non-adjacent vertices in cur_splex and v
int *critical_missing;	//critical_missing[v]: the number of non-adjacent critical vertices in cur_splex and v
int *unadj_satu_with_M2;
int *M2_type;
RandAccessList *vec_M1;
RandAccessList *vec_M2;
RandAccessList *cur_splex;
RandAccessList *cur_cand;
clock_t start_time;
int flag_remove = 1;
//fast construct
bool* can_add;
bool* is_satu;
int* remove_flag;
/*final result*/
int best_size;
int* best_plex;
clock_t best_time;
int best_found_iter;
int total_iter;
clock_t total_time;
int total_start_pass;
int* local_opt_score;
int local_best_size = 0;
int cycle_iter;
RandAccessList* ral_init(int capacity) {
	RandAccessList *ral = new RandAccessList;
	ral->vlist = new int[capacity];
	ral->vpos = new int[capacity];
	ral->vnum = 0;
	ral->capacity = capacity;
	for (int i = 0; i < capacity; i++) {
		ral->vpos[i] = capacity;
	}
	return ral;
}

void ral_add(RandAccessList *ral, int vid) {
	assert(ral->vpos[vid] >= ral->vnum || ral->vlist[ral->vpos[vid]] != vid);
	ral->vlist[ral->vnum] = vid;
	ral->vpos[vid] = ral->vnum;
	ral->vnum++;
}
void ral_delete(RandAccessList *ral, int vid) {
	// assert(ral->vpos[vid] < ral->vnum);
	int last_id = ral->vlist[ral->vnum - 1];
	int id_pos = ral->vpos[vid];
	ral->vlist[id_pos] = last_id;
	ral->vpos[last_id] = id_pos;
	ral->vnum--;
	//	ral->vpos[vid] = ral->vnum; /*It is not obligatory*/
}

void ral_clear(RandAccessList *ral) {
	ral->vnum = 0;
}
void ral_release(RandAccessList *ral) {
	delete[] ral->vlist;
	delete[] ral->vpos;
	delete ral;
}

int cmpfunc(const void * a, const void * b)
{
	return (*(int*)a - *(int*)b);
}

void ral_showList(RandAccessList *ral, FILE *f) {
	fprintf(f, "Total %d: ", ral->vnum);
	int *tmp_lst = new int[ral->capacity];
	memcpy(tmp_lst, ral->vlist, ral->vnum * sizeof(int));
	qsort(tmp_lst, ral->vnum, sizeof(int), cmpfunc);
	for (int i = 0; i < ral->vnum; i++) {
		fprintf(f, "%d ", tmp_lst[i]);
	}
	fprintf(f, "\n");
}


inline bool ver_exist(RandAccessList *ral, int ver, int flagg) {
	bool flag = (ral->vpos[ver] < ral->vnum &&  ral->vlist[ral->vpos[ver]] == ver);
	if (flag == false)
		cout << param_graph_file_name << "ver exist " << param_seed << " " << flagg << " " << restart_pass << " " << cycle_iter << endl;
	return flag;
	//return 1;
}

inline bool ver_exist2(RandAccessList *ral, int ver) {
	bool flag = (ral->vpos[ver] < ral->vnum &&  ral->vlist[ral->vpos[ver]] == ver);
	return flag;
	//return 1;
}
/**
 * load instances from  2nd Dimacs competetion
 *
 */
int load_clq_instance(char* filename){
	ifstream infile(filename);
	char line[1024];
	char tmps1[1024];
	char tmps2[1024];
	/*graph adjacent matrix*/
	int **gmat;
	if (!infile.is_open()){
		fprintf(stderr,"Can not find file %s\n", filename);
		return 0;
	}

	infile.getline(line,  1024);
	while (line[0] != 'p')	infile.getline(line,1024);
	sscanf(line, "%s %s %d %d", tmps1, tmps2, &org_vnum, &org_enum);

	gmat = new int*[org_vnum];
	for (int i = 0; i < org_vnum; i++){
		gmat[i] = new int[org_vnum];
		memset(gmat[i], 0, sizeof(int) * org_vnum);
	}

	int ecnt = 0;
	org_v_edge_cnt = new int[org_vnum];
	while (infile.getline(line, 1024)){
		int v1,v2;
		if (strlen(line) == 0)
			continue;
		if (line[0] != 'e')
			fprintf(stderr, "ERROR in line %d\n", ecnt+1);
		sscanf(line, "%s %d %d", tmps1, &v1, &v2);
		v1--,v2--;
		gmat[v1][v2] = 1;
		gmat[v2][v1] = 1;
		org_v_edge_cnt[v1]++;
		org_v_edge_cnt[v2]++;
		ecnt++;
	}
	assert(org_enum == ecnt);
	org_fmt = 0;
	org_vid = new int[org_vnum];
	org_v_adj_vertex = new int*[org_vnum];
	no_org_v_adj_vertex = new int*[org_vnum];	//no adjacent vertices matrix data
	for (int i = 0; i < org_vnum; i++){
		int adj_cnt = 0;
		int no_adj_cnt=0;
		org_v_adj_vertex[i] = new int[org_v_edge_cnt[i]];
		no_org_v_adj_vertex[i]= new int[org_vnum-org_v_edge_cnt[i]-1];  //no adjacent vertices data of i
		org_vid[i] = i+1;
		for (int j = 0; j < org_vnum; j++){
			if (gmat[i][j] == 1)
				org_v_adj_vertex[i][adj_cnt++] = j;
			else if (gmat[i][j] == 0 && (i != j)){
				no_org_v_adj_vertex[i][no_adj_cnt++] = j;		//no adjacent vertices data added to matrix
			}
		}
		assert(org_v_edge_cnt[i] == adj_cnt);
		delete[] gmat[i];
	}
	delete[] gmat;
	return 1;
}

/*load instances from  Stanford Large Network Dataset Collection
 * URL: http://snap.stanford.edu/data/*/
int load_snap_instance(char *filename){
	ifstream infile(filename);
	char line[1024];
	vector<pair<int,int> > *pvec_edges = new vector<pair<int, int> >();
	const int CONST_MAX_VE_NUM = 9999999;
	if (!infile.is_open()){
		fprintf(stderr,"Can not find file %s\n", filename);
		return 0;
	}
	int max_id = 0;
	int from, to;
	//int max_id = 0;
	while (infile.getline(line, 1024)){
		char *p = line;
		while (*p ==' ' && *p !='\0') p++;
		if (*p != '#'){
			stringstream ss(line);
			ss >> from >> to;
			//printf("%d %d \n", from, to);
			pvec_edges->push_back(make_pair(from, to));
			if (from > max_id)
				max_id = from;
			else if (to > max_id)
				max_id = to;
		}
	}
	int *newid = new int[max_id+1];
	org_vid = new int[CONST_MAX_VE_NUM];
		//init the newid map
	for (int i = 0; i < max_id+1; i++){
		newid[i] = -1;
	}

	int v_num = 0;
	//count edges,resign ids from 1...v_num
	for (int i = 0; i < (int)pvec_edges->size(); i++){
		from = (pvec_edges->at(i)).first;
		if (newid[from] == -1){
			newid[from] = v_num;
			org_vid[v_num] = from;
			v_num++;
		}
		(pvec_edges->at(i)).first = newid[from];
		to = (pvec_edges->at(i)).second;
		if (newid[to] == -1) {
			newid[to] = v_num;
			org_vid[v_num] = to;
			v_num++;
		}
		(pvec_edges->at(i)).second = newid[to];
	}
	org_vnum = v_num;
	int *estimate_edge_cnt = new int[org_vnum];
	memset(estimate_edge_cnt, 0, sizeof(int) * org_vnum);
	for (int i = 0; i < (int)pvec_edges->size(); i++){
		from = (pvec_edges->at(i)).first;
		to = (pvec_edges->at(i)).second;
		estimate_edge_cnt[from]++;
		estimate_edge_cnt[to]++;
	}
	org_enum = 0;
	org_v_edge_cnt = new int[org_vnum];
	memset(org_v_edge_cnt, 0, sizeof(int) * org_vnum);
	org_v_adj_vertex = new int*[org_vnum];
	for (int i = 0; i < (int)pvec_edges->size(); i++){
		from = (pvec_edges->at(i)).first;
		to = (pvec_edges->at(i)).second;
		if (from == to) continue; //self-edges are abandoned
		if (org_v_edge_cnt[from] == 0){
			org_v_adj_vertex[from] = new int[estimate_edge_cnt[from]];
		}
		int exist1 = 0;
		for (int j = 0; j < org_v_edge_cnt[from]; j++){
			if (org_v_adj_vertex[from][j] == to){
				exist1 = 1;
				break;
			}
		}
		if (!exist1){
			org_v_adj_vertex[from][org_v_edge_cnt[from]] = to;
			org_v_edge_cnt[from]++;
			org_enum++;
		}

		if (org_v_edge_cnt[to] == 0){
			org_v_adj_vertex[to] = new int[estimate_edge_cnt[to]];
		}
		int exist2 = 0;
		for (int j = 0; j < org_v_edge_cnt[to]; j++){
			if (org_v_adj_vertex[to][j] == from){
				exist2 = 1;
				break;
			}
		}
		if (!exist2){
			org_v_adj_vertex[to][org_v_edge_cnt[to]] = from;
			org_v_edge_cnt[to]++;
			org_enum++;
		}
	}
	assert(org_enum %2 == 0);
	org_enum = org_enum/2;
	printf("load SNAP graph %s with vertices %d, edges %d (directed edges %d)\n",
			basename(param_graph_file_name), org_vnum, org_enum, (int)pvec_edges->size());
	//print_org_graph();
	//build v_
	delete[] newid;
	delete[] estimate_edge_cnt;
	delete pvec_edges;
	return 1;
}

/*load instances of metis format from
 * instances are download from
 * http://www.cc.gatech.edu/dimacs10/downloads.shtml
 */
int load_metis_instance(char* filename){
	ifstream infile(filename);
	string line;

	if (!infile.is_open()){
		fprintf(stderr,"Can not find file %s\n", filename);
		exit(0);
	}
	org_fmt = 0;
	//ignore comment
	getline(infile, line);
	while(line.length()==0 || line[0] == '%')
		getline(infile, line);
	stringstream l1_ss(line);
	l1_ss >> org_vnum >> org_enum >> org_fmt;
	cout << line << endl;
	if (org_fmt == 100){
		cerr << "self-loops and/or multiple edges need to be considered";
		return 0;
	}
	/*allocate graph memory*/
	org_v_edge_cnt = new int[org_vnum];
	org_v_adj_vertex = new int*[org_vnum];
	org_vid = new int[org_vnum];

	int v_no = 0;
	int *vlst = new int[org_vnum];
	int n_adj = 0;
	//ignore the char
	//read the end of first line
	//infile.getline(line, LINE_LEN);
	while (getline(infile, line) && v_no < org_vnum){
		if (infile.fail()){
			fprintf(stderr, "Error in read file, vno %d\n",v_no);
			exit(-1);
		}
		stringstream ss(line);
		int ve_adj;

		n_adj = 0;
 		while (ss >> ve_adj){
			vlst[n_adj++] = ve_adj-1;
		}
		org_v_edge_cnt[v_no] = n_adj;
		if (n_adj == 0)
			org_v_adj_vertex[v_no] = NULL;
		else{
			org_v_adj_vertex[v_no] = new int[n_adj];
			memcpy(org_v_adj_vertex[v_no], vlst, sizeof(int) * n_adj);
		}
		org_vid[v_no] = v_no+1;
		v_no++;
	}
	assert(v_no == org_vnum);
	printf("load metis graph %s with vertices %d, edges %d \n",
			basename(param_graph_file_name), org_vnum, org_enum);
	//debug
	delete[] vlst;
	infile.close();
	return 1;
}
/*recursively remove all the vertices with degree less or equal than reduce_deg,
 * reconstructing the reduced graph  */
void reduce_graph(int reduce_deg){
	// cout << "the reduce_deg is " << reduce_deg << "\t";
	int critical_num = 0;
	
	int rm_cnt = 0;
	int *rmflag = new int[red_vnum];
	int *remain_edge_count = new int[red_vnum];
	/*Avoid unnecessary reduction*/
	if (red_min_deg > reduce_deg){
		return;
	}
	cout << "reduce" << param_graph_file_name << param_seed << endl;
	memcpy(remain_edge_count, red_v_edge_cnt, sizeof(int) * red_vnum);
	memset(rmflag, 0, sizeof(int) * red_vnum);
	queue<int> rm_que;

	for (int idx = 0; idx < red_vnum; idx++){
		if (remain_edge_count[idx] <= reduce_deg){
			if (remain_edge_count[idx] == reduce_deg) { critical_num++; }//count
			rmflag[idx] = 1;
			rm_que.push(idx);
			rm_cnt++;
		}
	}
	while(!rm_que.empty()){
		int idx = rm_que.front();
		rm_que.pop();
		for (int i = 0; i < red_v_edge_cnt[idx]; i++){
			int adjv = red_v_adj_vertex[idx][i];
			remain_edge_count[adjv]--;
			if (!rmflag[adjv] && remain_edge_count[adjv] <= reduce_deg){
				if (remain_edge_count[adjv] == reduce_deg) { critical_num++; }//count
				rmflag[adjv] = 1;
				rm_que.push(adjv);
				rm_cnt++;
			}
		}
	}
	/*rebuild the reduced graph*/
	int rest_vnum = red_vnum - rm_cnt;
	int *new_id = new int[red_vnum];
	/*Mapp the vertex id to the original id(in the initial graph) */
	int *org_id = new int[rest_vnum];
	int count = 0; //count the rest vertices
	int *last_freq = new int[rest_vnum];
	int *last_momentum = new int[rest_vnum];
	double *last_potential = new double[rest_vnum];
	//resign new id to the rest of the vertices
	for (int idx = 0; idx < red_vnum; idx++){
		if (!rmflag[idx]){
			new_id[idx] = count;
			org_id[count] = red_orgid[idx];
			last_freq[count] = freq[idx];
			last_momentum[count] = momentum[idx];
			count++;
		}
	}
	/**/
	int min_deg = LARGE_INT;
	int n_edges = 0;
	int **new_adj_tbl = new int*[rest_vnum];
	int *new_edge_count = new int[rest_vnum];
	for (int idx_prev = 0; idx_prev < red_vnum; idx_prev++){
		if (rmflag[idx_prev]){
			delete[] red_v_adj_vertex[idx_prev];
			continue;
		}
		int idx_new = new_id[idx_prev];
		new_adj_tbl[idx_new] = new int[remain_edge_count[idx_prev]];
		new_edge_count[idx_new] = remain_edge_count[idx_prev];
		int cnt = 0;
		for (int i = 0; i < red_v_edge_cnt[idx_prev]; i++){
			int vi_adj = red_v_adj_vertex[idx_prev][i];
			if (!rmflag[vi_adj]){
				new_adj_tbl[idx_new][cnt++] = new_id[vi_adj];
				n_edges++;
			}
		}
		if (new_edge_count[idx_new] < min_deg)
			min_deg = new_edge_count[idx_new];
		assert(cnt == remain_edge_count[idx_prev]);
		delete[] red_v_adj_vertex[idx_prev];
	}
	assert(n_edges % 2 == 0);
	/*assign to the new graph*/
	delete[] red_v_adj_vertex;
	red_v_adj_vertex = new_adj_tbl;
	delete[] red_v_edge_cnt;
	red_v_edge_cnt = new_edge_count;

	red_vnum = rest_vnum;
	red_enum = n_edges/2;

	/*reset the minimum deg*/
	red_min_deg = min_deg;

	/*reset the original id map*/
	delete[] red_orgid;
	red_orgid = org_id;

	/*reset the last join record*/
	delete[] freq;
	freq = last_freq;

	delete[] momentum;
	momentum = last_momentum;

	delete[] new_id;
	delete[] rmflag;
	delete[] remain_edge_count;

}

/*reinitial the data for a new start of the algorithm*/
void init_search(){
	/*copy the graph to reduce graph*/
	red_vnum = org_vnum;
	red_enum = org_enum;
	red_orgid = new int[red_vnum];
	red_v_edge_cnt = new int[red_vnum];
	red_v_adj_vertex = new int*[red_vnum];
	red_min_deg = LARGE_INT;
	memcpy(red_v_edge_cnt, org_v_edge_cnt, sizeof(int) * org_vnum);
	for (int v = 0; v < red_vnum; v++){
		red_orgid[v] = v;
		// potential[v] = 0;
		red_v_adj_vertex[v] = new int[red_v_edge_cnt[v]];
		memcpy(red_v_adj_vertex[v], org_v_adj_vertex[v], sizeof(int) * org_v_edge_cnt[v]);
		if (red_v_edge_cnt[v] < red_min_deg)
			red_min_deg = red_v_edge_cnt[v];
	}

	/*init search data*/
	cur_iter = 0;
	cur_c_deg = new int[red_vnum];
	is_in_c = new int[red_vnum];
	memset(cur_c_deg, 0, sizeof(int) * red_vnum);
	memset(is_in_c, 0, sizeof(int) * red_vnum);
	freq = new int[red_vnum];
	memset(freq, 0, sizeof(int) * red_vnum);
	freq1 = new int[red_vnum];
	memset(freq1, 0, sizeof(int) * red_vnum);
	momentum = new int [red_vnum];
	memset(momentum, 0, sizeof(int) * red_vnum);
	local_opt_score = new int[red_vnum];
	memset(local_opt_score, 0, sizeof(int) * red_vnum);
	time_stamp = new int[org_vnum];
	cur_splex = ral_init(red_vnum);
	cur_cand = ral_init(red_vnum);
	vec_M1 = ral_init(red_vnum);
	vec_M2 = ral_init(red_vnum);
	M2_type = new int[red_vnum];
	memset(M2_type, 0, sizeof(int) * org_vnum);

	threshold = new int[red_vnum];
	deposit = new int[red_vnum];

	can_add = new bool[red_vnum];
	is_satu = new bool[red_vnum];

	plex_missing=new int[org_vnum]; 		//initialize plex_missing and critical_missing matrix 
	critical_missing=new int[org_vnum];
	remove_flag = new int[red_vnum];
	memset(remove_flag, 0, sizeof(int) * red_vnum);
	memset(plex_missing, 0, sizeof(int) * org_vnum);
	memset(critical_missing, 0, sizeof(int) * org_vnum);
	unadj_satu_with_M2=new int[org_vnum]; 
	memset(unadj_satu_with_M2, -1, sizeof(int) * org_vnum);
	/*init best found data*/
	best_size = 0;
	/*We could reduce this size*/
	best_plex = new int[red_vnum];
	best_time = 0;
	best_found_iter = 0;
	total_start_pass = 0;

	start_time = clock();
    srand((unsigned int)param_seed);
	if (org_vnum == 400 || org_vnum == 4000) {
		flag_remove = 0;
	}
}

void restart_search(){
	local_best_size = 0;
	memset(time_stamp, 0, sizeof(int) * org_vnum);
	memset(plex_missing, 0, sizeof(int) * org_vnum);		//reset plex_missing and critical_missing when restart
	memset(critical_missing, 0, sizeof(int) * org_vnum);
	memset(unadj_satu_with_M2, -1, sizeof(int) * org_vnum);
	memset(M2_type, 0, sizeof(int) * org_vnum);
	memset(cur_c_deg, 0, sizeof(int) * red_vnum);
	memset(is_in_c, 0, sizeof(int) * red_vnum);
	memset(local_opt_score, 0, sizeof(int) * red_vnum);
	ral_clear(cur_splex);
	ral_clear(cur_cand);	
	for (int i = 0; i < red_vnum; ++i){
		threshold[i] = 1;
		deposit[i] = 1;
	}
}

void add_without_change_deposit(int v){	//new add function without change deposit
	aaa = 1;
	assert(!is_in_c[v]);
	is_in_c[v] = 1;

	ral_add(cur_splex, v);
	if (cur_c_deg[v] > 0){
		ral_delete(cur_cand, v);
	}

	for (int i = 0; i < red_v_edge_cnt[v]; i++){
		int adjv = red_v_adj_vertex[v][i];
		cur_c_deg[adjv]++;
		if (!is_in_c[adjv] && cur_c_deg[adjv] == 1){
			ral_add(cur_cand, adjv);
		}
	}
	ral_delete(vec_M1, v);
	for(int i = 0; i < org_vnum-org_v_edge_cnt[v]-1; i++){    //no-adjacent vertices times iteration 
		int no_adjv = no_org_v_adj_vertex[v][i];	//no-adjacent vertex[i] of v
		plex_missing[no_adjv]++;	//plex_missing++
		if(!is_in_c[no_adjv]){			//no-adjacent vertices out of cur_splex,if their plex_missing up to bound, will change to another set
			if(plex_missing[no_adjv] == param_s+1) {  //k->k+1
				if(critical_missing[no_adjv] <= 1){	//M2 to M3
					if (ver_exist(vec_M2, no_adjv, 2)){
						ral_delete(vec_M2, no_adjv);
					}
				}
			}
			else if(plex_missing[no_adjv] == param_s) {    // k-1->k
				if(critical_missing[no_adjv] == 0) {	//M1 to M2
					if (ver_exist(vec_M1, no_adjv, 1)){
						ral_delete(vec_M1, no_adjv);
						ral_add(vec_M2, no_adjv);
						M2_type[no_adjv] = 1;
					}
				}
			}
		}
		else if(is_in_c[no_adjv]){	//if no-adjacent vertex belong to cur_splex and plex_missing change to bound,their no-adjacent vertices will change critical_missing
			if(plex_missing[no_adjv] == param_s-1) { //plex_missing up to k-1,become a critical node,critical_missing of no-adjacent vertices of it will plus 1
				for(int j=0; j < org_vnum-org_v_edge_cnt[no_adjv]-1; j++) {	//no-adjacent vertices of critical node
					int second_no_adjv = no_org_v_adj_vertex[no_adjv][j];
					critical_missing[second_no_adjv]++;	//critical_missing++
					unadj_satu_with_M2[second_no_adjv]=no_adjv;
					if(critical_missing[second_no_adjv] == 2 && !is_in_c[second_no_adjv] && second_no_adjv != v){  // M2 to M3
						if(plex_missing[second_no_adjv] <= param_s){
							if (ver_exist(vec_M2, second_no_adjv, 2)) {
								ral_delete(vec_M2, second_no_adjv);
							}
						}
					}
					else if(critical_missing[second_no_adjv] == 1 && !is_in_c[second_no_adjv] && second_no_adjv != v){ //c_m(critical_missing)0->1
						if(plex_missing[second_no_adjv] == param_s){   //p_m(plex_missing)==k,but c_m(0->1),M2(second=1) to M2(second=0)
							if (ver_exist(vec_M2, second_no_adjv, 2) && M2_type[second_no_adjv] == 1) {
								M2_type[second_no_adjv] = 0;
							}
						}
						else if(plex_missing[second_no_adjv] <= param_s-1){    //M1 to M2
							if (ver_exist(vec_M1, second_no_adjv, 1)) {
								ral_delete(vec_M1, second_no_adjv);
								ral_add(vec_M2, second_no_adjv);
								M2_type[second_no_adjv] = 0;
							}
						}
					}
				}
			}
			else if(plex_missing[no_adjv] == param_s){	//k-1->k, the vertex will not be a critical node, the c_m of no-adjacent vertices of it will minus 1
				for(int j=0; j < org_vnum-org_v_edge_cnt[no_adjv]-1; j++) {	//no-adjacent vertices of this node
					int second_no_adjv = no_org_v_adj_vertex[no_adjv][j];
					critical_missing[second_no_adjv]--;	//c_m--
					if(unadj_satu_with_M2[second_no_adjv]== no_adjv){
						unadj_satu_with_M2[second_no_adjv]=-1;
					}
					if(critical_missing[second_no_adjv] == 1 && !is_in_c[second_no_adjv] && second_no_adjv != v){ //c_m 2->1
						if(plex_missing[second_no_adjv] <= param_s){    //M3 to M2
							ral_add(vec_M2, second_no_adjv);
							M2_type[second_no_adjv] = 0;
							if (unadj_satu_with_M2[second_no_adjv] == -1) {
								for (int k = 0; k < org_vnum - org_v_edge_cnt[second_no_adjv] - 1; k++) {
									if (is_in_c[no_org_v_adj_vertex[second_no_adjv][k]] && plex_missing[no_org_v_adj_vertex[second_no_adjv][k]] == param_s - 1) {
										unadj_satu_with_M2[second_no_adjv] = no_org_v_adj_vertex[second_no_adjv][k];
										break;
									}
								}
							}
						}
					}
					else if(critical_missing[second_no_adjv] == 0 && !is_in_c[second_no_adjv] && second_no_adjv != v){  //c_m 1->0
						if(plex_missing[second_no_adjv] <= param_s-1){	// M2 to M1
							if (ver_exist(vec_M2, second_no_adjv, 2)) {
								ral_delete(vec_M2, second_no_adjv);
								ral_add(vec_M1, second_no_adjv);
							}
							
						}
						else if(plex_missing[second_no_adjv] == param_s){ // M2 to M2
							if (ver_exist(vec_M2, second_no_adjv, 2) && M2_type[second_no_adjv] == 0) {
								M2_type[second_no_adjv] = 1;
							}
						}
					}
				}
			}

		}
	}
	if(plex_missing[v] == param_s-1){	//v is a critical node and now add to cur_splex
		for(int j=0; j < org_vnum-org_v_edge_cnt[v]-1; j++) {	//no-adjacent vertices of it
			int no_adjv = no_org_v_adj_vertex[v][j];
			critical_missing[no_adjv]++;	//c_m++
			unadj_satu_with_M2[no_adjv]=v;
			if(critical_missing[no_adjv] == 2 && !is_in_c[no_adjv] && no_adjv != v){ //c_m 1->2  
				if(plex_missing[no_adjv] <= param_s){// M2 to M3
					if (ver_exist(vec_M2, no_adjv, 2)) {
						ral_delete(vec_M2, no_adjv);
					}
				}
			}
			else if(critical_missing[no_adjv] == 1 && !is_in_c[no_adjv] && no_adjv != v){ //c_m 0->1
				if(plex_missing[no_adjv] == param_s){   //M2 to M2
					if (ver_exist(vec_M2, no_adjv, 2) && M2_type[no_adjv] == 1) {
						M2_type[no_adjv] = 0;
					}
				}
				else if(plex_missing[no_adjv] <= param_s-1){    //M1 to M2
					if (ver_exist(vec_M1, no_adjv, 1)) {
						ral_delete(vec_M1, no_adjv);
						ral_add(vec_M2, no_adjv);
						M2_type[no_adjv] = 0;
					}
				}
			}
		}
	}
}

int bms_thre_pro(int v) {
	int count = 0;
	for (int i = 0; i < 50; i++) {
		if (freq[v] > freq[rand() % red_vnum]) {
			count++;
		}
	}
	return count;//若count很大说明freq很小，则threshold不++
}

void add_cur_vertex(int v){		//old add function used in init solution create
	aaa = 2;
	assert(!is_in_c[v]);

	is_in_c[v] = 1;
	deposit[v] = 0;
	threshold[v] = (threshold[v] + 1) % 3;
	ral_add(cur_splex, v);
	if (cur_c_deg[v] > 0){
		ral_delete(cur_cand, v);
	}

	for (int i = 0; i < red_v_edge_cnt[v]; i++){
		int adjv = red_v_adj_vertex[v][i];
		++deposit[adjv];
		cur_c_deg[adjv]++;
		if (!is_in_c[adjv] && cur_c_deg[adjv] == 1){
			ral_add(cur_cand, adjv);
		}
	}
}

void add_cur_vertex_new(int v){   //new add function used in M1 and M3
	aaa = 3;
	assert(!is_in_c[v]);

	is_in_c[v] = 1;
	deposit[v] = 0;
	threshold[v] = threshold[v] + 1;
	if (threshold[v] >= 3) {
		int countttt = bms_thre_pro(v);
		if (countttt > 40) {
			threshold[v] = 1;
		}
		else {
			threshold[v] = 0;
		}
	}
	ral_add(cur_splex, v);
	if (cur_c_deg[v] > 0){
		ral_delete(cur_cand, v);
	}

	for (int i = 0; i < red_v_edge_cnt[v]; i++){
		int adjv = red_v_adj_vertex[v][i];
		++deposit[adjv];
		cur_c_deg[adjv]++;
		if (!is_in_c[adjv] && cur_c_deg[adjv] == 1){
			ral_add(cur_cand, adjv);
		}
	}
	for(int i = 0; i < org_vnum-org_v_edge_cnt[v]-1; i++){   //no-adjacent vertices times iteration  
		int no_adjv = no_org_v_adj_vertex[v][i];	//no-adjacent vertex[i] of v
		plex_missing[no_adjv]++;	//plex_missing++
		if(!is_in_c[no_adjv]){	//no-adjacent vertices out of cur_splex,if their plex_missing up to bound, will change to another set
			if(plex_missing[no_adjv] == param_s+1) {  //k->k+1
				if(critical_missing[no_adjv] <= 1){	//M2 to M3
					if (ver_exist(vec_M2, no_adjv, 2)) {
						ral_delete(vec_M2, no_adjv);
					}
				}
			}
			else if(plex_missing[no_adjv] == param_s) {    // k-1->k
				if(critical_missing[no_adjv] == 0) {	//M1 to M2
					if (ver_exist(vec_M1, no_adjv, 1)) {
						ral_delete(vec_M1, no_adjv);
						ral_add(vec_M2, no_adjv);
						M2_type[no_adjv] = 1;
					}
				}
			}
		}
		else if(is_in_c[no_adjv]){	//if no-adjacent vertex belong to cur_splex and plex_missing change to bound,their no-adjacent vertices will change critical_missing
			if(plex_missing[no_adjv] == param_s-1) {	//plex_missing up to k-1,become a critical node,critical_missing of no-adjacent vertices of it will plus 1
				for(int j=0; j < org_vnum-org_v_edge_cnt[no_adjv]-1; j++) {	//no-adjacent vertices of critical node
					int second_no_adj = no_org_v_adj_vertex[no_adjv][j];
					critical_missing[second_no_adj]++;	//critical_missing++
					unadj_satu_with_M2[second_no_adj]=no_adjv;
					if(critical_missing[second_no_adj] == 2 && !is_in_c[second_no_adj] && second_no_adj != v){  // M2 to M3
						if(plex_missing[second_no_adj] <= param_s){
							if (ver_exist(vec_M2, second_no_adj, 2)) {
								ral_delete(vec_M2, second_no_adj);
							}
						}
					}
					else if(critical_missing[second_no_adj] == 1 && !is_in_c[second_no_adj] && second_no_adj != v){	//c_m(critical_missing)0->1
						if(plex_missing[second_no_adj] == param_s){   //p_m(plex_missing)==k,but c_m(0->1),M2(second=1) to M2(second=0)
							if (ver_exist(vec_M2, second_no_adj, 2) && M2_type[second_no_adj] == 1) {
								M2_type[second_no_adj] = 0;
							}
						}
						else if(plex_missing[second_no_adj] <= param_s-1){    //M1 to M2
							if (ver_exist(vec_M1, second_no_adj, 1)) {
								ral_delete(vec_M1, second_no_adj);
								ral_add(vec_M2, second_no_adj);
								M2_type[second_no_adj] = 0;
							}
						}
					}
				}
			}
			else if(plex_missing[no_adjv] == param_s){	///是否为必须要删除的点？？？移除顶点是否考虑 plex_missing == param_s	//k-1->k, the vertex will not be a critical node, the c_m of no-adjacent vertices of it will minus 1
				for(int j=0; j < org_vnum-org_v_edge_cnt[no_adjv]-1; j++) {		//no-adjacent vertices of this node
					int second_no_adj = no_org_v_adj_vertex[no_adjv][j];
					critical_missing[second_no_adj]--;		//c_m-
					if(unadj_satu_with_M2[second_no_adj]==no_adjv){
						unadj_satu_with_M2[second_no_adj]=-1;
					}
					if(critical_missing[second_no_adj] == 1 && !is_in_c[second_no_adj] && second_no_adj != v){		//c_m 2->1
						if(plex_missing[second_no_adj] <= param_s){    //M3 to M2
							ral_add(vec_M2, second_no_adj);
							M2_type[second_no_adj] = 0;
							int curid = second_no_adj;
							if (unadj_satu_with_M2[curid] == -1) {
								for (int k = 0; k < org_vnum - org_v_edge_cnt[curid] - 1; k++) {
									if (is_in_c[no_org_v_adj_vertex[curid][k]] && plex_missing[no_org_v_adj_vertex[curid][k]] == param_s - 1) {
										unadj_satu_with_M2[curid] = no_org_v_adj_vertex[curid][k];
										break;
									}
								}
							}
						}
					}
					else if(critical_missing[second_no_adj] == 0 && !is_in_c[second_no_adj] && second_no_adj != v){  // M2 to M1
						if(plex_missing[second_no_adj] == param_s){ // M2 to M2
							if (ver_exist(vec_M2, second_no_adj, 2) && M2_type[second_no_adj] == 0) {
								M2_type[second_no_adj] = 1;
							}
							
						}
						else if(plex_missing[second_no_adj] <= param_s-1){
							if (ver_exist(vec_M2, second_no_adj, 2)) {
								ral_delete(vec_M2, second_no_adj);
								ral_add(vec_M1, second_no_adj);
							}
						}
					}
				}
			}
		}
	}
	if(plex_missing[v] == param_s-1){ 	//v is a critical node and now add to cur_splex
		for(int j=0; j < org_vnum-org_v_edge_cnt[v]-1; j++) {	//no-adjacent vertices of it
			int second_no_adj = no_org_v_adj_vertex[v][j];
			critical_missing[second_no_adj]++;	//c_m++
			unadj_satu_with_M2[second_no_adj]=v;
			if(critical_missing[second_no_adj] == 2 && !is_in_c[second_no_adj] && second_no_adj != v){  // M2 to M3
				if(plex_missing[second_no_adj] <= param_s){
					if (ver_exist(vec_M2, second_no_adj, 2)) {
						ral_delete(vec_M2, second_no_adj);
					}
				}
			}
			else if(critical_missing[second_no_adj] == 1 && !is_in_c[second_no_adj] && second_no_adj != v){
				if(plex_missing[second_no_adj] == param_s){   //M2 to M2
					if (ver_exist(vec_M2, second_no_adj, 2) && M2_type[second_no_adj] == 1) {
						M2_type[second_no_adj] = 0;
					}
				}
				else if(plex_missing[second_no_adj] <= param_s-1){    //M1 to M2
					if (ver_exist(vec_M1, second_no_adj, 1)) {
						ral_delete(vec_M1, second_no_adj);
						ral_add(vec_M2, second_no_adj);
						M2_type[second_no_adj] = 0;
					}
				}
			}
		}
	}
}

void remove_cur_vertex_new(int v){   //new remove function
	aaa = 4;
	assert(is_in_c[v]);

	int tmpflag=0;
	if(plex_missing[v] == param_s-1) tmpflag=1;		//judge the vertex will be removed whether is a critical,if true, tmpflag=1

	is_in_c[v] = 0;
	deposit[v] = 0;
	ral_delete(cur_splex, v);
	if (cur_c_deg[v] > 0){
		ral_add(cur_cand, v);
	}

	for (int i = 0; i < red_v_edge_cnt[v]; i++){
		int adjv = red_v_adj_vertex[v][i];
		cur_c_deg[adjv]--;
		if (!is_in_c[adjv] && cur_c_deg[adjv] == 0){
			ral_delete(cur_cand, adjv);
		}
	}
	for(int i = 0; i < org_vnum-org_v_edge_cnt[v]-1; i++){     	//no-adjacent vertices times iteration 
		int no_adjv = no_org_v_adj_vertex[v][i];	//no-adjacent vertex[i] of v
		plex_missing[no_adjv]--;		//plex_missing--
		if(!is_in_c[no_adjv]){	//no-adjacent vertices out of cur_splex,if their plex_missing down to bound, will change to another set
			if(plex_missing[no_adjv] == param_s) {    //M3 to M2
				if(critical_missing[no_adjv] == 1){
					ral_add(vec_M2, no_adjv);
					M2_type[no_adjv] = 0;
					if (unadj_satu_with_M2[no_adjv] == -1) {
						for (int k = 0; k < org_vnum - org_v_edge_cnt[no_adjv] - 1; k++) {
							if (is_in_c[no_org_v_adj_vertex[no_adjv][k]] && plex_missing[no_org_v_adj_vertex[no_adjv][k]] == param_s - 1) {
								unadj_satu_with_M2[no_adjv] = no_org_v_adj_vertex[no_adjv][k];
								break;
							}
						}
					}
				}
				else if(critical_missing[no_adjv] == 0) {
					ral_add(vec_M2, no_adjv);
					M2_type[no_adjv] = 1;
				}
			}
			else if(plex_missing[no_adjv] == param_s-1) {  //M2 to M1
				if(critical_missing[no_adjv] == 0){
					if (ver_exist(vec_M2, no_adjv, 2)) {
						ral_delete(vec_M2, no_adjv);
						ral_add(vec_M1, no_adjv);
					}
				}
			}
		}

		else if(is_in_c[no_adjv]){	//if no-adjacent vertex belong to cur_splex and plex_missing change to bound(k->k-1),their no-adjacent vertices will change critical_missing
			if(plex_missing[no_adjv] == param_s-2) {	//plex_missing down to k-2,become a not critical node,critical_missing of no-adjacent vertices of it will minus 1
				for(int j=0; j < org_vnum-org_v_edge_cnt[no_adjv]-1; j++) {		//no-adjacent vertices of critical node
					int second_no_adjv = no_org_v_adj_vertex[no_adjv][j];
					critical_missing[second_no_adjv]--;		//critical_missing--
					if(unadj_satu_with_M2[second_no_adjv]==no_adjv){
						unadj_satu_with_M2[second_no_adjv]=-1;
					}
					if(critical_missing[second_no_adjv] == 1 && !is_in_c[second_no_adjv] && second_no_adjv != v){
						if(plex_missing[second_no_adjv] <= param_s){    //M3 to M2
							ral_add(vec_M2, second_no_adjv);
							M2_type[second_no_adjv] = 0;
							if (unadj_satu_with_M2[second_no_adjv] == -1) {
								for (int k = 0; k < org_vnum - org_v_edge_cnt[second_no_adjv] - 1; k++) {
									if (is_in_c[no_org_v_adj_vertex[second_no_adjv][k]] && plex_missing[no_org_v_adj_vertex[second_no_adjv][k]] == param_s - 1) {
										unadj_satu_with_M2[second_no_adjv] = no_org_v_adj_vertex[second_no_adjv][k];
										break;
									}
								}
							}
						}
					}
					else if(critical_missing[second_no_adjv] == 0 && !is_in_c[second_no_adjv] && second_no_adjv != v){  // M2 to M1
						if(plex_missing[second_no_adjv] == param_s){ // M2 to M2
							if (ver_exist(vec_M2, second_no_adjv, 2) && M2_type[second_no_adjv] == 0) {
								M2_type[second_no_adjv] = 1;
							}
						}
						else if(plex_missing[second_no_adjv] <= param_s-1){
							if (ver_exist(vec_M2, second_no_adjv, 2)) {
								ral_delete(vec_M2, second_no_adjv);
								ral_add(vec_M1, second_no_adjv);
							}
						}
					}
				}
			}
			else if(plex_missing[no_adjv] == param_s-1) {	//plex_missing down to k-1(k->k-1),become a critical node,critical_missing of no-adjacent vertices of it will plus 1
				for(int j=0; j < org_vnum-org_v_edge_cnt[no_adjv]-1; j++) {	//no-adjacent vertices of critical node
					int second_no_adjv = no_org_v_adj_vertex[no_adjv][j];
					critical_missing[second_no_adjv]++;	//c_m++
					unadj_satu_with_M2[second_no_adjv]=no_adjv;
					if(critical_missing[second_no_adjv] == 2 && !is_in_c[second_no_adjv] && second_no_adjv != v){  // M2 to M3
						if(plex_missing[second_no_adjv] <= param_s){
							if (ver_exist(vec_M2, second_no_adjv, 2)) {
								ral_delete(vec_M2, second_no_adjv);
							}
						}
					}
					else if(critical_missing[second_no_adjv] == 1 && !is_in_c[second_no_adjv] && second_no_adjv != v){
						if(plex_missing[second_no_adjv] == param_s){   //M2 to M2
							if (ver_exist(vec_M2, second_no_adjv, 2) && M2_type[second_no_adjv] == 1) {
								M2_type[second_no_adjv] = 0;
							}
							
						}
						else if(plex_missing[second_no_adjv] <= param_s-1){    //M1 to M2
							if (ver_exist(vec_M1, second_no_adjv, 1)) {
								ral_delete(vec_M1, second_no_adjv);
								ral_add(vec_M2, second_no_adjv);
								M2_type[second_no_adjv] = 0;
							}
						}
					}
				}
			}
		}
	}
	if(tmpflag == 1) {	//if the node removed is a critical node before,c_m of no-adjacent vertices of it will minus 1
		for(int j=0; j < org_vnum-org_v_edge_cnt[v]-1; j++) {	//no-adjacent vertices of critical node
			int second_no_adjv = no_org_v_adj_vertex[v][j];
			critical_missing[second_no_adjv]--;	//c_m--
			if(unadj_satu_with_M2[second_no_adjv] == v){
				unadj_satu_with_M2[second_no_adjv]=-1;
			}
			if(critical_missing[second_no_adjv] == 1 && !is_in_c[second_no_adjv] && second_no_adjv != v){
				if(plex_missing[second_no_adjv] <= param_s){    //M3 to M2
					ral_add(vec_M2, second_no_adjv);
					M2_type[second_no_adjv] = 0;
					if (unadj_satu_with_M2[second_no_adjv] == -1) {
						for (int k = 0; k < org_vnum - org_v_edge_cnt[second_no_adjv] - 1; k++) {
							if (is_in_c[no_org_v_adj_vertex[second_no_adjv][k]] && plex_missing[no_org_v_adj_vertex[second_no_adjv][k]] == param_s - 1) {
								unadj_satu_with_M2[second_no_adjv] = no_org_v_adj_vertex[second_no_adjv][k];
								break;
							}
						}
					}
				}
			}
			else if(critical_missing[second_no_adjv] == 0 && !is_in_c[second_no_adjv] && second_no_adjv != v){  // M2 to M1
				if(plex_missing[second_no_adjv] <= param_s-1){
					if (ver_exist(vec_M2, second_no_adjv, 2)) {
						ral_delete(vec_M2, second_no_adjv);
						ral_add(vec_M1, second_no_adjv);
					}
					
				}
				else if(plex_missing[second_no_adjv] == param_s){ // M2 to M2
					if (ver_exist(vec_M2, second_no_adjv, 2) && M2_type[second_no_adjv] == 0) {
						M2_type[second_no_adjv] = 1;
					}
					
				}
			}
		}
	}	
	if(critical_missing[v]==0 && plex_missing[v]<=param_s-1) 	//add the removed vertex to the M(1/2/3) set
		ral_add(vec_M1, v);
	else if(critical_missing[v]==0 && plex_missing[v]==param_s) {
		ral_add(vec_M2, v);
		M2_type[v] = 1;
	}
	else if(critical_missing[v] == 1 && plex_missing[v] <= param_s) {
		ral_add(vec_M2, v);
		M2_type[v] = 0;
		if(unadj_satu_with_M2[v] == -1){
			for(int k = 0; k < org_vnum-org_v_edge_cnt[v]-1; k++){
				if(is_in_c[no_org_v_adj_vertex[v][k]] && plex_missing[no_org_v_adj_vertex[v][k]] == param_s-1){
					unadj_satu_with_M2[v]=no_org_v_adj_vertex[v][k];
					break;
				}
			}
		}

	}
}

#define Is_Saturated(v) (cur_c_deg[v] == cur_splex->vnum - param_s)
#define Is_Overflow(v) (cur_c_deg[v] < cur_splex->vnum - param_s)

int get_saturate_size(){
	int size = 0;
	for (int i = 0; i < cur_splex->vnum; i++){
		int vin = cur_splex->vlist[i];
		if (Is_Saturated(vin)){
			size++;
		}	
	}
	return size;
}

void record_best(){
	best_size = cur_splex->vnum;
	for (int i = 0; i < cur_splex->vnum; i++){
		int v = cur_splex->vlist[i];
		best_plex[i] = red_orgid[v];
	}
	best_time = clock();
	best_found_iter = cur_iter;
}

void record_local_best() {
	for (int i = 0; i < red_vnum; i++) {
		local_opt_score[i] = cur_c_deg[i];
	}
}

void fast_init_solution(){
	vector<int> vec_M1;
	bool* adj_flag = new bool[red_vnum];
	bool* satu_adj_flag = new bool[red_vnum];
	memset(can_add, 0, sizeof(bool) * red_vnum);
	memset(is_satu, 0, sizeof(bool) * red_vnum);
	int leastfreq = LARGE_INT;
	for (int i = 0; i < 100; i++){
		int randv = rand() % red_vnum;
		if (freq[randv] < leastfreq){
			vec_M1.clear();
			vec_M1.push_back(randv);
			leastfreq = freq[randv];
		}else if (freq[randv] == leastfreq){
			vec_M1.push_back(randv);
		}
	}
	assert (!vec_M1.empty());
	int vstart = vec_M1[rand()%vec_M1.size()];
	add_cur_vertex(vstart);   //use old add function here
	for (int i = 0; i < red_v_edge_cnt[vstart]; ++i){
		can_add[red_v_adj_vertex[vstart][i]] = true;
	}
	while(1){
		vec_M1.clear();
		leastfreq = LARGE_INT;
		for (int i = 0; i < red_vnum; ++i){
			if (can_add[i]){
				if (freq[i] < leastfreq){
					vec_M1.clear();
					vec_M1.push_back(i);
					leastfreq = freq[i];
				}else if(freq[i] == leastfreq){
					vec_M1.push_back(i);
				}
			}
		}
		if (vec_M1.empty()){
			break;
		}
		assert (!vec_M1.empty());
		int vadd = vec_M1[rand()%vec_M1.size()];
		add_cur_vertex(vadd);
		can_add[vadd] = false;
		memset(adj_flag, 0, sizeof(bool) * red_vnum);
		for (int i = 0; i < red_v_edge_cnt[vadd]; ++i){
			adj_flag[red_v_adj_vertex[vadd][i]] = true;
		}
		for (int i = 0; i < red_vnum; ++i){
			if (!adj_flag[i]){
				if (is_in_c[i] && !is_satu[i] && Is_Saturated(i)){
					// i is just saturated
					is_satu[i] = true;
					memset(satu_adj_flag, 0, sizeof(bool) * red_vnum);
					for (int j = 0; j < red_v_edge_cnt[i]; ++j){
						satu_adj_flag[red_v_adj_vertex[i][j]] = true;
					}
					for (int j = 0; j < red_vnum; ++j){
						if (can_add[j] && !satu_adj_flag[j]){
							can_add[j] = false;
						}
					}
				}else if (!is_in_c[i] && can_add[i] && cur_c_deg[i] <= cur_splex->vnum - param_s){
					can_add[i] = false;
				}
			}
		}

	}
	while( cur_splex->vnum < param_s){
		int vrand = rand() % red_vnum;
		while (is_in_c[vrand]) vrand = rand() % red_vnum;
		add_cur_vertex(vrand);
	}
	if (cur_splex->vnum > best_size){
		record_best();
	}
}

/*find the only unadjacent saturated vertex of v*/
int find_unadj_satu(int v){
	assert(!is_in_c[v]);
	int *mark = new int[cur_splex->vnum];
	int v_rt = -1;
	memset(mark, 0, sizeof(int) * cur_splex->vnum);
	for (int i = 0; i < red_v_edge_cnt[v]; i++){//邻居在解中
		int vcur = red_v_adj_vertex[v][i];
		if (is_in_c[vcur]){
			mark[cur_splex->vpos[vcur]] = 1;
		}
	}
	for (int i = 0; i < cur_splex->vnum; i++){
		int vin = cur_splex->vlist[i];
		if (Is_Saturated(vin) && mark[i] == 0){
			v_rt = vin;
			break;
		}
	}
	delete[] mark;
	return v_rt;
}

/*get a random vertex of C\N_C(v)*/
int random_unadj_with_exception(int v, int vexception){
	assert(!is_in_c[v]);
	int *mark = new int[cur_splex->vnum];
	vector<int> vec_unadj;
	memset(mark, 0, sizeof(int) * cur_splex->vnum);
	for (int i = 0; i < red_v_edge_cnt[v]; i++){
		int vcur = red_v_adj_vertex[v][i];
		if (is_in_c[vcur]){
			mark[cur_splex->vpos[vcur]] = 1;
		}
	}
	for (int i = 0; i < cur_splex->vnum; i++){
		if (!mark[i] && cur_splex->vlist[i] != vexception)
			vec_unadj.push_back(cur_splex->vlist[i]);
	}
	delete[] mark;
	if (vec_unadj.size() == 0)
		return -1;
	else
		return vec_unadj[rand() % vec_unadj.size()];
}

int get_most_momentum(vector<int>& perturb_set){
	int min_momentum_2 = -999999999;
	vector<int> cand;
	vector<int> cand2;
	cand.push_back(perturb_set[0]);
	int min_momentum = momentum[perturb_set[0]];
	for (int i = 1, size = perturb_set.size(); i < size; ++i){
		int score = momentum[perturb_set[i]];
		if (score < min_momentum && min_momentum_2 == -999999999) {
			min_momentum_2 = score;
			cand2.push_back(perturb_set[i]);
		}
		else if (score > min_momentum) {
			min_momentum_2 = min_momentum;
			cand2.clear();
			for (int i = 0; i < cand.size(); i++) {
				cand2.push_back(cand[i]);
			}

			min_momentum = score;
			cand.clear();
			cand.push_back(perturb_set[i]);
		}
		else if (score < min_momentum && score > min_momentum_2) {
			min_momentum_2 = score;
			cand2.clear();
			cand2.push_back(perturb_set[i]);
		}
		else if (score == min_momentum) {
			cand.push_back(perturb_set[i]);
		}
		else if (score == min_momentum_2) {
			cand2.push_back(perturb_set[i]);
		}
	}
	if (cand2.size() > 0 && rand() % 10 < 3)
		return cand2[rand() % cand2.size()];
	return cand[rand() % cand.size()];
}

int get_most_momentum_add(vector<int>& perturb_set) {
	double min_momentum = -999999999;
	vector<int> cand;
	for (int i = 0, size = perturb_set.size(); i < size; ++i) {
		double curr_score = (double)cur_c_deg[perturb_set[i]] + 0.2 * ((double)cur_c_deg[perturb_set[i]] / red_vnum)*(red_v_edge_cnt[perturb_set[i]] - cur_c_deg[perturb_set[i]]);
		if (curr_score > min_momentum) {
			min_momentum = cur_c_deg[perturb_set[i]];
			cand.clear();
			cand.push_back(perturb_set[i]);
		}
		else if (curr_score == min_momentum) {
			cand.push_back(perturb_set[i]);
		}
	}
	return cand[rand() % cand.size()];
}

void check_missing() {
	for (int i = 0; i < red_vnum; i++) {
		int temp_p = 0;
		int temp_c = 0;
		for (int j = 0; j < org_vnum - org_v_edge_cnt[i] - 1; j++) {	//no-adjacent vertices of critical node
			int second_no_adjv = no_org_v_adj_vertex[i][j];
			if (is_in_c[second_no_adjv])
				temp_p++;
			if (is_in_c[second_no_adjv] && Is_Saturated(second_no_adjv))
				temp_c++;
		}
		if (temp_p != plex_missing[i])
			cout << "temp_p != plex_missing[i]" << endl;
		if (temp_c != critical_missing[i])
			cout << "temp_c != critical_missing[i]" << endl;
	}
}

void check_vector() {
	for (int i = 0; i < red_vnum; i++) {
		if (is_in_c[i])
			continue;
		if (critical_missing[i] == 0 && plex_missing[i] <= param_s - 1) {
			if (!ver_exist2(vec_M1, i)) {
				cout << i << " " << " !ver_exist2(vec_M1, i)" << endl;
			}
		}
		else if (critical_missing[i] == 0 && plex_missing[i] == param_s) {
			if (!ver_exist2(vec_M2, i)) {
				cout << i << " " << " !ver_exist2(vec_M2, i)1" << endl;
			}
		}
		else if (critical_missing[i] == 1 && plex_missing[i] <= param_s) {
			if (!ver_exist2(vec_M2, i)) {
				cout << i << " " << " !ver_exist2(vec_M2, i)2" << endl;
			}
		}
	}
	for (int j = 0; j < vec_M1->vnum; j++) {
		int ver = vec_M1->vlist[j];
		if (critical_missing[ver] == 0 && plex_missing[ver] <= param_s - 1) {

		}
		else {
			cout << "M1 ELSE" << endl;
		}
	}
	for (int j = 0; j < vec_M2->vnum; j++) {
		int ver = vec_M2->vlist[j];
		if ((critical_missing[ver] == 0 && plex_missing[ver] == param_s) || (critical_missing[ver] == 1 && plex_missing[ver] <= param_s)) {

		}
		else {
			cout << "M2 ELSE" << endl;
		}
	}
}

void push_vertex_tabu(int vpush) {
	add_cur_vertex_new(vpush);
	time_stamp[vpush] = cycle_iter;
	/*repair*/
	int idx = 0;
	while (idx < cur_splex->vnum) {
		int vin = cur_splex->vlist[idx];
		if (Is_Overflow(vin)) {
			remove_cur_vertex_new(vin);
			freq[vin]++;
			freq1[vin]++;
			time_stamp[vin] = cycle_iter;
		}
		else
			idx++;
	}
}

void tabu_based_search(){
	ral_clear(vec_M1);
	ral_clear(vec_M2);
	cycle_iter = 0;
	int cycle_best = cur_splex->vnum;
	int fixed = -1;
	int vpush;

	for(int i = 0; i < cur_splex->vnum; i++){	//update plex_missing and critical_missing through the initial solution 
		for(int j=0; j < org_vnum-org_v_edge_cnt[cur_splex->vlist[i]]-1; j++){
			int no_adjv1 = no_org_v_adj_vertex[cur_splex->vlist[i]][j];
			plex_missing[no_adjv1]++;
			if(plex_missing[no_adjv1] == param_s-1 && is_in_c[no_adjv1]){
				for(int k=0; k < org_vnum-org_v_edge_cnt[no_adjv1]-1; k++){
					critical_missing[no_org_v_adj_vertex[no_adjv1][k]]++;
					unadj_satu_with_M2[no_org_v_adj_vertex[no_adjv1][k]] = no_adjv1;
				}
			}
		}
	}
	for(int i = 0; i < cur_cand->vnum; i++){	//add cur_cand vertices to M(1/2/3) with the principle
		if(critical_missing[cur_cand->vlist[i]]==0 && plex_missing[cur_cand->vlist[i]]<=param_s-1) {
			ral_add(vec_M1, cur_cand->vlist[i]);
		}
		else if(critical_missing[cur_cand->vlist[i]]==0 && plex_missing[cur_cand->vlist[i]]==param_s) {
			ral_add(vec_M2, cur_cand->vlist[i]);
			M2_type[cur_cand->vlist[i]] = 1;
		}
		else if(critical_missing[cur_cand->vlist[i]] == 1 && plex_missing[cur_cand->vlist[i]] <= param_s) {
			ral_add(vec_M2, cur_cand->vlist[i]);
			M2_type[cur_cand->vlist[i]] = 0;
		}
	}
	while (1){
		int end = 0;
		int non_improve_iter = 0;
		while (!end){
			vector<int> vec_M1_tmp;
			for(int i = 0; i < vec_M1->vnum; ++i){ //choose (deposit>=threshold) vertices from vec_M1
				if(deposit[vec_M1->vlist[i]] >= threshold[vec_M1->vlist[i]] || cycle_iter - time_stamp[vec_M1->vlist[i]] > 4 || cur_splex->vnum + 1 > best_size) {
					vec_M1_tmp.push_back(vec_M1->vlist[i]);
				}
			}
			total_add += vec_M1->vnum;
			sai_add += vec_M1_tmp.size();
			count_add++;
			if (!vec_M1_tmp.empty()){				// add vertex from vec_M1_tmp
				int vadd;
				vadd = get_most_momentum_add(vec_M1_tmp);
				ral_delete(vec_M1, vadd);
				add_cur_vertex_new(vadd);   //new add function
				time_stamp[vadd] = cycle_iter;
				if (cur_splex->vnum > cycle_best){
					cycle_best = cur_splex->vnum;
				}
				freq[vadd]++;
				freq1[vadd]++;
				vec_M1_tmp.clear();    //clear vec_M1_tmp
			}
			else{
				vector<int> vec_M2_tmp;
				for(int i = 0; i < vec_M2->vnum; ++i){ 	//choose (deposit>=threshold) vertices from vec_M2_out
					if(deposit[vec_M2->vlist[i]] >= threshold[vec_M2->vlist[i]]) {
						vec_M2_tmp.push_back(vec_M2->vlist[i]);		
					}
				}
				total_swap += vec_M2->vnum;
				sai_swap += vec_M2_tmp.size();
				count_swap++;
				if (!vec_M2_tmp.empty()){
					int pvswp;
					pvswp = get_most_momentum(vec_M2_tmp);
					int vswp_in;
					if (M2_type[pvswp] == 0){ //type 1
						if(unadj_satu_with_M2[pvswp] != -1)
							vswp_in=unadj_satu_with_M2[pvswp];
						else vswp_in = find_unadj_satu(pvswp);//pvswp不与解中的某个饱和顶点相邻
					}else{ //type 2
						vswp_in = random_unadj_with_exception(pvswp, fixed);//在解中，不与pvswp相邻，且不为fixed的顶点
					}
					if (flag_remove == 0 && (vswp_in == fixed || vswp_in == -1)) {
						end = 1;
					}
					else{
						remove_cur_vertex_new(vswp_in);   //new remove function
						time_stamp[vswp_in] = cycle_iter;
						add_without_change_deposit(pvswp);//SwapSet交换之后仍能保证当前解为feasible.
						time_stamp[pvswp] = cycle_iter;
						freq[vswp_in]++;
						freq[pvswp]++;
						freq1[vswp_in]++;
						freq1[pvswp]++;
						vec_M2_tmp.clear();		//clear vec_M2_tmp set
					}
				}
				else{
					if (flag_remove == 0) {
						vector<int> cand;
						int remove_count = 0;
						for (int i = 0; i < cur_splex->vnum; i++) {
							int ver = cur_splex->vlist[i];
							if (Is_Saturated(ver)) {
								cand.push_back(ver);
							}
						}
						while (cand.size() != 0) {
							int vvv = cand[rand() % cand.size()];
							remove_cur_vertex_new(vvv);
							time_stamp[vvv] = cycle_iter;
							cand.clear();
							remove_count++;
							remove_flag[vvv] = 1;
							for (int i = 0; i < cur_splex->vnum; i++) {
								int ver = cur_splex->vlist[i];
								if (Is_Saturated(ver)) {
									cand.push_back(ver);
								}
							}
						}
						vector<int> ttttt;
						vector<int> ttttt2;//移除的顶点
						int count1 = 0;
						do {
							ttttt.clear();
							ttttt2.clear();
							for (int i = 0; i < vec_M1->vnum; ++i) { //choose (deposit>=threshold) vertices from vec_M1
								if (remove_flag[vec_M1->vlist[i]] == 0) {
									ttttt.push_back(vec_M1->vlist[i]);
								}
								else {
									ttttt2.push_back(vec_M1->vlist[i]);
								}
							}

							if (ttttt.size() != 0 && ttttt2.size() != 0) {
								if (rand() % 10 < 5) {
									int vvv = ttttt[rand() % ttttt.size()];
									ral_delete(vec_M1, vvv);
									add_cur_vertex_new(vvv);
									time_stamp[vvv] = cycle_iter;
								}
								else {
									int vvv = ttttt2[rand() % ttttt2.size()];
									ral_delete(vec_M1, vvv);
									add_cur_vertex_new(vvv);
									time_stamp[vvv] = cycle_iter;
								}
								count1++;
								remove_count--;
							}
						} while (ttttt.size() != 0 && ttttt2.size() != 0);
						memset(remove_flag, 0, sizeof(int) * red_vnum);
						count_per++;
						non_improve_iter = 0;
					}
					else {
						end = 1;
					}
				}
			}
			if (cur_splex->vnum > local_best_size) {
				local_best_size = cur_splex->vnum;
				record_local_best();
			}
			if (cur_splex->vnum > best_size){
				record_best();
				if (best_size == param_best) //reach optimum
					goto ts_stop;
				non_improve_iter = 0;
			}else{
				non_improve_iter++;
			}
			if (non_improve_iter > param_s * best_size) {
				if (flag_remove == 0) {
					vector<int> m4;
					vector<int> m3;
					int max_freq = 0;
					for (int i = 0; i < red_vnum; i++) {
						if (is_in_c[i] == 1 || ver_exist2(vec_M1, i) || ver_exist2(vec_M2, i)) {
							continue;
						}
						m4.push_back(i);
						if (critical_missing[i] == 0 && plex_missing[i] <= param_s - 1) {
							cout << "critical_missing[i] == 0 && plex_missing[i] <= param_s - 1 " << restart_pass << " " << cycle_iter << " " << param_graph_file_name << endl;
						}
						else if (critical_missing[i] == 0 && plex_missing[i] == param_s) {
							cout << "critical_missing[i] == 0 && plex_missing[i] == param_s " << restart_pass << " " << cycle_iter << " " << param_graph_file_name << endl;
						}
						else if (critical_missing[i] == 1 && plex_missing[i] <= param_s) {
							cout << "critical_missing[cur_cand->vlist[i]] == 1 && plex_missing[cur_cand->vlist[i]] <= param_s " << restart_pass << " " << cycle_iter << " " << param_graph_file_name << endl;
						}

						if (deposit[i] >= threshold[i] && freq1[i] > max_freq) {
							m3.clear();
							m3.push_back(i);
							max_freq = freq1[i];
						}
						else if (deposit[i] >= threshold[i] && freq1[i] == max_freq) {
							m3.push_back(i);
						}
					}
					int vpush1;
					if (m3.empty()) {
						vpush1 = m4[rand() % m4.size()];
					}
					else {
						vpush1 = m3[rand() % m3.size()];
					}
					freq1[vpush1] = 0;
					fixed = vpush1;
					push_vertex_tabu(vpush1);
				}
				else {
					vector<int> cand;
					int remove_count = 0;
					for (int i = 0; i < cur_splex->vnum; i++) {
						int ver = cur_splex->vlist[i];
						if (Is_Saturated(ver)) {
							cand.push_back(ver);
						}
					}
					while (cand.size() != 0) {
						int vvv = cand[rand() % cand.size()];
						remove_cur_vertex_new(vvv);
						time_stamp[vvv] = cycle_iter;
						cand.clear();
						remove_count++;
						remove_flag[vvv] = 1;
						for (int i = 0; i < cur_splex->vnum; i++) {
							int ver = cur_splex->vlist[i];
							if (Is_Saturated(ver)) {
								cand.push_back(ver);
							}
						}
					}
					vector<int> ttttt;
					vector<int> ttttt2;//移除的顶点
					int count1 = 0;
					do {
						ttttt.clear();
						ttttt2.clear();
						for (int i = 0; i < vec_M1->vnum; ++i) { //choose (deposit>=threshold) vertices from vec_M1
							if (remove_flag[vec_M1->vlist[i]] == 0) {
								ttttt.push_back(vec_M1->vlist[i]);
							}
							else {
								ttttt2.push_back(vec_M1->vlist[i]);
							}
						}

						if (ttttt.size() != 0 && ttttt2.size() != 0) {
							if (rand() % 10 < 5) {
								int vvv = ttttt[rand() % ttttt.size()];
								ral_delete(vec_M1, vvv);
								add_cur_vertex_new(vvv);
								time_stamp[vvv] = cycle_iter;
							}
							else {
								int vvv = ttttt2[rand() % ttttt2.size()];
								ral_delete(vec_M1, vvv);
								add_cur_vertex_new(vvv);
								time_stamp[vvv] = cycle_iter;
							}
							count1++;
							remove_count--;
						}
					} while (ttttt.size() != 0 && ttttt2.size() != 0);
					memset(remove_flag, 0, sizeof(int) * red_vnum);
					count_per++;
					non_improve_iter = 0;
				}
			}
			cur_iter++;
			cycle_iter++;

			if(cycle_iter > param_cycle_iter ||
					(float)(clock() - start_time) / CLOCKS_PER_SEC > param_max_seconds){
				goto ts_stop;
			}
		}
		if (cur_splex->vnum == 0) {
			break;
		}
		//vpush = cur_splex->vlist[rand() % cur_splex->vnum];//随机移除一个顶点。
		//remove_cur_vertex_new(vpush);
		//time_stamp[vpush] = cycle_iter;
		if (flag_remove == 0) {
			vector<int> m4;
			vector<int> m3;
			int max_freq = 0;
			for (int i = 0; i < red_vnum; i++) {
				if (is_in_c[i] == 1 || ver_exist2(vec_M1, i) || ver_exist2(vec_M2, i)) {
					continue;
				}
				m4.push_back(i);
				if (critical_missing[i] == 0 && plex_missing[i] <= param_s - 1) {
					cout << "critical_missing[i] == 0 && plex_missing[i] <= param_s - 1 " << restart_pass << " " << cycle_iter << " " << param_graph_file_name << endl;
				}
				else if (critical_missing[i] == 0 && plex_missing[i] == param_s) {
					cout << "critical_missing[i] == 0 && plex_missing[i] == param_s " << restart_pass << " " << cycle_iter << " " << param_graph_file_name << endl;
				}
				else if (critical_missing[i] == 1 && plex_missing[i] <= param_s) {
					cout << "critical_missing[cur_cand->vlist[i]] == 1 && plex_missing[cur_cand->vlist[i]] <= param_s " << restart_pass << " " << cycle_iter << " " << param_graph_file_name << endl;
				}

				if (deposit[i] >= threshold[i] && freq1[i] > max_freq) {
					m3.clear();
					m3.push_back(i);
					max_freq = freq1[i];
				}
				else if (deposit[i] >= threshold[i] && freq1[i] == max_freq) {
					m3.push_back(i);
				}
			}
			int vpush1;
			if (m3.empty()) {
				vpush1 = m4[rand() % m4.size()];
			}
			else {
				vpush1 = m3[rand() % m3.size()];
			}
			freq1[vpush1] = 0;
			fixed = vpush1;
			push_vertex_tabu(vpush1);
		}
		else {
			vpush = cur_splex->vlist[rand() % cur_splex->vnum];//随机移除一个顶点。
			remove_cur_vertex_new(vpush);
			time_stamp[vpush] = cycle_iter;
		}
	}
ts_stop:
	return;
}

/*TODO:Entrance of the whole search*/
void search_main(){
	restart_pass = 0;
	/*Initial data structure*/
	init_search();
	// start from a new solution
	fast_init_solution();
	for (int i = 0; i < red_vnum; i++) {
		momentum[i] += cur_c_deg[i];
	}
	reduce_graph(best_size - param_s);
	if (red_vnum == 0 ){
		goto end;
	}
	while (1){
		restart_search();
		fast_init_solution();
		tabu_based_search();
		if (best_size == param_best)
			goto end;
		if (best_size - param_s >= red_min_deg){
			reduce_graph(best_size - param_s);
		}
		if (red_vnum <= best_size )
			goto end;
		restart_pass++;
		for (int i = 0; i < red_vnum; i++) {
			momentum[i] += local_opt_score[i];
		}
		if ((float)(clock() - start_time) / CLOCKS_PER_SEC > param_max_seconds ){
			break;
		}
	}
end:
	total_start_pass = restart_pass;
	total_time = clock();
	total_iter = cur_iter;
}

void report_result(){
	cout << " Instance: " << param_graph_file_name << endl;
	cout << " k: " << param_s << endl;
	cout << " seed: " << param_seed << endl;
	cout << " best solution size: " << best_size << endl;
	cout << " best solution time: " << (float)((best_time - start_time) / CLOCKS_PER_SEC) << endl;
	cout << " Best solution: ";
	for (int i = 0; i < best_size; i++) {
		cout << best_plex[i] << " ";
	}
}

void showUsage(){
	fprintf(stderr, "splex -f <filename> -s <paramete s> -t <max seconds>  [-o optimum object] [-c seed]");
}

void read_params(int argc, char **argv){
	int hasFile = 0;
	int hasTimeLimit = 0;
	int hasS = 0;
	for (int i = 1; i < argc; i+=2){
		if (argv[i][0] != '-' || argv[i][2] != 0){
			showUsage();
			exit(0);
		}else if (argv[i][1] == 'f'){
			strncpy(param_graph_file_name, argv[i+1],1000);
			hasFile = 1;
		}else if(argv[i][1] == 's'){
			param_s = atoi(argv[i+1]);
			if (param_s >= 1)	hasS = 1;
		}else if (argv[i][1] == 'o'){
			param_best = atoi(argv[i+1]);
		}else if (argv[i][1] == 'c'){
			param_seed = atoi(argv[i+1]);
		}else if(argv[i][1] == 't'){
			param_max_seconds = atoi(argv[i+1]);
			hasTimeLimit = 1;
		}
	}
	/*check parameters*/
	if (!hasFile){
		fprintf(stderr,"No file name\n");
		showUsage();
		exit(1);
	}
	if (!hasTimeLimit){
		fprintf(stderr,"No time limit \n");
		showUsage();
		exit(1);
	}
	if (!hasS){
		fprintf(stderr,"No paramete s\n");
		showUsage();
		exit(1);
	}
}

int check_solution(){
	int *mark = new int[org_vnum];
	memset(mark, 0, sizeof(int) * org_vnum);
	for (int i = 0; i < best_size; i++){
		mark[best_plex[i]] = 1;
	}
	for (int i = 0; i < best_size; i++){
		int v = best_plex[i];
		int indeg = 0;
		for (int j = 0; j < org_v_edge_cnt[v]; j++){
			int vadj = org_v_adj_vertex[v][j];
			if (mark[vadj])
				indeg++;
		}
		if (indeg < best_size - param_s)
			return 0;
	}
	delete[] mark;
	return 1;
}

const char* file_suffix(char* filename){
    const char *dot = strrchr(filename, '.');
    if(!dot || dot == filename) return "";
    return dot + 1;
}

int main(int argc, char** argv) {
	int load = 0;

	read_params(argc, argv);
	const char* fileext = file_suffix(param_graph_file_name);
	//printf("%s\n",fileext);
	if (0 == strcmp(fileext, "graph")){
		load = load_metis_instance(param_graph_file_name);
	}else if(0 == strcmp(fileext, "txt")){
		load = load_snap_instance(param_graph_file_name);
	}else if(0 == strcmp(fileext, "clq")){
		load = load_clq_instance(param_graph_file_name);
	}
	if (load != 1){
		fprintf(stderr, "failed in loading graph %s\n",param_graph_file_name);
		exit(-1);
	}

	search_main();
	int chk = check_solution();
	if (chk == 0){
		fprintf(stderr,"ERROR! Final solution is infeasible\n");
		exit(-1);
	}
	report_result();
}
