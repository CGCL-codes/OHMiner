#pragma once
#include<vector>
#include "map"
#include "unordered_map"
#include "common.h"
#include "string"
#include "pattern.h"
#include "set"
#include "unordered_set"

class Pattern;
struct Stable{
public:
    std::vector<int> hedge_list;
    std::unordered_map<int, std::vector<int>> sedge;                      //sedge表示的是一系列的超边
    std::unordered_map<int, std::vector<int>> invert_index;              //倒排索引
};

class hypergraph{
public:
    //没有hyperedge 包含超边的信息啊
    hypergraph(int e_cnt, int v_cnt);
    int e_cnt;
    int v_cnt;
    std::vector<int>hedge_to_kid;
    std::vector<int>v_to_lable;
    std::unordered_map<int , Stable> HyperGraph_table;


    //hyperedge2node、node2hyperedge、hyperedge_adj是为了random_walk
    std::vector< std::vector<int> >hyperedge2node;        //the hyperedge of the query
    std::vector< std::vector<int> >node2hyperedge;        //表示的是node 和 hyperedge之间的关系
    std::vector< std::vector<int> > hyperedge_adj;        //表示的是超边之间的连接关系

    std::vector< std::vector<int> >csr;
    std::vector< std::vector<int> >deg;
    //std::vector< std::unordered_map<int, std::vector<std::pair<int, int>>> > hyperedgeAdj;
    std::vector< std::unordered_map<int, std::vector<int> >> hyperedgeAdj;

    std::unordered_map<int64_t, std::vector<int> > hyper_inter;


    // int max_size_DAG = 0;
    std::vector<std::vector<int>> DAG;
    std::vector<std::vector<std::vector<int>>> DAG_vector;
    static double valid_time;
    static double generate_time;
    static double intersect_time;
    static long long candidate_num;

    int64_t overlap_aware_pattern_matching(Pattern* p);
    int64_t overlap_aware_pattern_matching_parallel(Pattern* p);
    int64_t overlap_aware_pattern_matching_parallel_openmp(Pattern *p);

    void Enumerate_overlap(std::vector<int>&matching_order, std::vector<int>& candidate, int depth, Pattern* p, std::vector<int>& subtraction_set, int64_t& count);
    void Enumerate_overlap_parallel(int thread_id, std::vector<int>&matching_order, std::vector<int>& candidate, int depth, Pattern* p, std::vector<int>& subtraction_set,std::vector<std::vector<int>> &DAG_temp);
    void Enumerate_overlap_parallel_openmp(int thread_id, std::vector<int> &matching_order, std::vector<int> &candidate, int depth , Pattern *p, std::vector<int> &subtraction_set,std::vector<std::vector<int>>& DAG_temp);
    
    void get_no_incident(int query_edge, std::vector<int>& v_no_incident, std::vector<int>& matching_order, std::vector<int>& subtraction_set, Pattern *p);
    void get_incident(int u, std::vector<int>&v_no_incident, int hedge, std::vector<int>& v_incident, Pattern* p, std::vector<int>& subtraction_set);

    bool isValidEmbedding_overlap(std::vector<int>&subtrction_set, std::vector<int>&matching_order,int depth, Pattern* p);
    bool isValidEmbedding_overlap_old(std::vector<int>&subtrction_set, std::vector<int>&matching_order,int depth, Pattern* p);
    bool isValidEmbedding_overlap_parallel(std::vector<int>&subtrction_set, std::vector<int>&matching_order,int depth, Pattern* p,std::vector<std::vector<int>> &DAG_temp);
    bool isValidEmbedding_overlap_parallel_old(std::vector<int>&subtrction_set, std::vector<int>&matching_order,int depth, Pattern* p,std::vector<std::vector<int>> &DAG_temp);

    void generate_hyperedge_candidates_overlap(Pattern* p, std::vector<int>& matching_order, int query_edge, std::vector<int>& subtraction_set, std::vector<int>& candidate);
    void generate_hyperedge_candidates_unlable(Pattern* p, std::vector<int>& matching_order, int query_edge, std::vector<int>& subtraction_set, std::vector<int>& candidate);
    
    void random_walk(int num_edges_to_sample, int min_vertex_count, int max_vertex_count, const char* pattern_path, const char* node_lable_path);
    void random_walk_dense(int num_edges_to_sample, int min_vertex_count, int max_vertex_count, const char* pattern_path, const char* node_lable_path);
};


