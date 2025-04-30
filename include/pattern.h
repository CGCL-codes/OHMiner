#pragma once
#include "vector"
#include "string"
#include "hypergraph.h"
#include <set>
#include <unordered_set>

class hypergraph;


class Pattern{
public:
    int v_cnt = 0;
    int e_cnt = 0;
    int overlap_node = 0;
    std::vector<int>v_to_lable;       //存储的是顶点和lable之间的映射关系(其中lable是从  datagraph中获取的)
    std::vector<int> edge_to_lable;    //存储边的lable
    std::vector< std::vector<int> >hyperedge2node;       //the hyperedge of the query
    std::vector< std::vector<int> >node2hyperedge;        //表示的是node 和 hyperedge之间的关系
    std::vector< std::vector<int> > hyperedge_adj;     //表示的是超边之间的连接关系
    std::vector< std::vector<int> > hyperedge_matrix; //表示超边之间的连接关系
    std::vector<int> vertives_num;
    std::unordered_map<int, std::vector<int>> pattern_profile;

    std::vector<int>rank;
    std::vector<int>matching_order;
    //表示的是对应超边包含的顶点
    std::vector< std::vector<int> >overlap_graph;
    std::unordered_set<int> leaf_node;
    std::vector< std::vector<int> >root_leaf;
    //表示overlap_graph之间超边之间的连接关系
    std::vector< std::vector<int> >overlap_graph_adj;
    std::vector< std::vector<int> >overlap_graph_adj_in;
    std::vector< std::vector<int> >EP;

    int matrix[100][100];


    int intersection(std::set<int>v_list, int hedge);

    void ComputeMatchingOrder(hypergraph* g);
    int get_cardinality(hypergraph* g, int hyperedge);
    void generate_overlap_graph();
    void generate_overlap_graph_sub();
    void expand_overlap_graph(const std::vector<int>&HE);

    void expand_overlap_graph_merge(const std::vector<std::pair<int, std::vector<int>>>& HE, int depth);
    void generate_execution_plan();
    bool is_tuan(std::vector<int>&);
    void intersect(const std::vector<int>&set1, const std::vector<int>&set2, std::vector<int>&result);
    void print();

};