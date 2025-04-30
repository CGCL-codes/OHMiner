#include "hypergraph.h"
#include <fstream>
#include <iostream>
#include <random>
#include <chrono>
#include <algorithm>
#include "unordered_set"
#include "vector"
#include "common.h"

#include "parallel.h"
#include <thread>
#include "omp.h"
#define __ALIGN 16

inline size_t ROUND_UP(size_t bytes) {
	return (((bytes) + __ALIGN-1) & (~(__ALIGN - 1)));
}

double hypergraph::valid_time = 0;
double hypergraph::generate_time = 0;
double hypergraph::intersect_time = 0;
long long hypergraph::candidate_num = 0;


hypergraph::hypergraph(int e_cnt = 0, int v_cnt = 0) {
    this->e_cnt = e_cnt;
    this->v_cnt = v_cnt;
}



void hypergraph::get_incident(int u, std::vector<int>& v_no_incident, int hedge, std::vector<int> &v_incident, Pattern *p, std::vector<int>& subtraction_set) {
    int lable = hedge_to_kid[hedge];
    std::vector<int> re = unordered_sub_hash(HyperGraph_table[lable].sedge[hedge], v_no_incident);
    for(int v : re){
        int d_v = 0;
        int d_u = 0;
        if(this->v_to_lable[v] != p->v_to_lable[u])
            continue;
        else{
            for(int e : subtraction_set){
                if(has_data(HyperGraph_table[hedge_to_kid[e]].sedge[e], HyperGraph_table[hedge_to_kid[e]].sedge[e].size(), v))
                    d_v++;
            }
            //提前计算好
            for(int e : p->node2hyperedge[u]){
                if(has_data(p->matching_order, subtraction_set.size(), e))
                    d_u++;
            }
            if(d_u == d_v)
                v_incident.push_back(v);
        }
    }
    return;
}

bool hypergraph::isValidEmbedding_overlap(std::vector<int>&subtrction_set, std::vector<int>&matching_order, int depth, Pattern* p) {
    for(auto e : p->EP[matching_order[depth]]){
        if(!p->overlap_graph_adj_in[e].empty()){
            int e1 = p->overlap_graph_adj_in[e][0], e2 = p->overlap_graph_adj_in[e][1];
            if(e1 < p->e_cnt &&  e2 < p->e_cnt){
                intersect_simd512(hyperedge2node[subtrction_set[p->rank[e1]]], hyperedge2node[subtrction_set[p->rank[e2]]],DAG[e]);
            }
            else
               intersect_simd512(DAG[p->overlap_graph_adj_in[e][0]], DAG[p->overlap_graph_adj_in[e][1]], DAG[e]);
            for(int i = 2; i <  p->overlap_graph_adj_in[e].size(); ++i){
                std::vector<int>temp;
                intersect_simd512(DAG[e], DAG[p->overlap_graph_adj_in[e][i]], temp);
                DAG[e].swap(temp);
            }
            if(p->overlap_graph[e].size() != DAG[e].size())
                return false;
        }
    }
    return true;
}

bool hypergraph::isValidEmbedding_overlap_old(std::vector<int>&subtrction_set, std::vector<int>&matching_order, int depth, Pattern* p) {
    for(auto e : p->EP[matching_order[depth]]){
        if(!p->overlap_graph_adj_in[e].empty()){
            DAG[e] = intersect(DAG[p->overlap_graph_adj_in[e][0]], DAG[p->overlap_graph_adj_in[e][1]]);
            for(int i = 2; i <  p->overlap_graph_adj_in[e].size(); ++i)
                DAG[e] = intersect(DAG[e], DAG[p->overlap_graph_adj_in[e][i]]);
        }
        else
            DAG[e] = hyperedge2node[subtrction_set.back()];
        if(DAG[e].size() != p->overlap_graph[e].size())
            return false;
    }
    return true;
}

bool hypergraph::isValidEmbedding_overlap_parallel_old(std::vector<int>&subtrction_set, std::vector<int>&matching_order, int depth, Pattern* p,std::vector<std::vector<int>>& DAG_temp) {
    for(auto e : p->EP[matching_order[depth]]){
        if(!p->overlap_graph_adj_in[e].empty()){
            DAG_temp[e] = intersect(DAG_temp[p->overlap_graph_adj_in[e][0]], DAG_temp[p->overlap_graph_adj_in[e][1]]);
            for(int i = 2; i <  p->overlap_graph_adj_in[e].size(); ++i)
                DAG_temp[e] = intersect(DAG_temp[e], DAG_temp[p->overlap_graph_adj_in[e][i]]);
        }
        else
            DAG_temp[e] = hyperedge2node[subtrction_set.back()];
        if(DAG_temp[e].size() != p->overlap_graph[e].size())
            return false;
    }
    return true;
}

bool hypergraph::isValidEmbedding_overlap_parallel(std::vector<int>&subtrction_set, std::vector<int>&matching_order, int depth, Pattern* p,std::vector<std::vector<int>>& DAG_temp) {
    for(auto e : p->EP[matching_order[depth]]){
        if(!p->overlap_graph_adj_in[e].empty()){
            int e1 = p->overlap_graph_adj_in[e][0], e2 = p->overlap_graph_adj_in[e][1];
            if(e1 < p->e_cnt &&  e2 < p->e_cnt){
                intersect_simd512(hyperedge2node[subtrction_set[p->rank[e1]]], hyperedge2node[subtrction_set[p->rank[e2]]],DAG_temp[e]);
            }
            else
               intersect_simd512(DAG_temp[p->overlap_graph_adj_in[e][0]], DAG_temp[p->overlap_graph_adj_in[e][1]], DAG_temp[e]);
            for(int i = 2; i <  p->overlap_graph_adj_in[e].size(); ++i){
                std::vector<int>temp;
                intersect_simd512(DAG_temp[e], DAG_temp[p->overlap_graph_adj_in[e][i]], temp);
                DAG_temp[e].swap(temp);
            }
            if(p->overlap_graph[e].size() != DAG_temp[e].size())
                return false;
        }
    }
    return true;
}

//matching order代表的是   原始的匹配顺序  subtraction_set中是 和原始的匹配顺序对应的  datagraph中的超边  我们要获取的是
void hypergraph::get_no_incident(int query_edge, std::vector<int> &v_no_incident, std::vector<int>& matching_order, std::vector<int>& subtraction_set, Pattern* p) {
    for (int i = 0; i < subtraction_set.size(); ++i) {
        if (p->hyperedge_matrix[query_edge][matching_order[i]] == 0) {
            int hlable = hedge_to_kid[subtraction_set[i]];
            for (int& v: HyperGraph_table[hlable].sedge[subtraction_set[i]]) {
                v_no_incident.push_back(v);
            }
        }
    }
    std::sort(v_no_incident.begin(), v_no_incident.end());
}

void hypergraph::generate_hyperedge_candidates_overlap(Pattern *p, std::vector<int> &matching_order, int query_edge,
                                                       std::vector<int> &subtraction_set, std::vector<int> &candidate) {
    
    std::vector<int>v_noincdt;
    std::vector< std::set<int> >candidate_union;

    get_no_incident(query_edge, v_noincdt, matching_order, subtraction_set, p);

    for(int i = 0; i < subtraction_set.size(); ++i){      
        if( !p->hyperedge_matrix[query_edge][matching_order[i]] )
            continue;
        else{
            for(int& u : p->hyperedge2node[matching_order[i]]){
                if( !has_data(p->hyperedge2node[query_edge], p->hyperedge2node[query_edge].size(), u))
                    continue;
                std::vector<int>v_incident;
                get_incident(u, v_noincdt, subtraction_set[i], v_incident, p, subtraction_set);
                std::set<int> tmp;
                for(int i = 0; i < v_incident.size(); ++i){
                    for(int val : this->HyperGraph_table[p->edge_to_lable[query_edge]].invert_index[v_incident[i]])
                        tmp.insert(val);  
                }
                candidate_union.push_back(tmp);
            }
        }
    }
    intersection(candidate_union, candidate);
    std::vector<int>unnei;
    for(int i = 0; i < subtraction_set.size(); ++i){
        if(p->hyperedge_matrix[query_edge][matching_order[i]]){
        }
        else
            unnei.push_back(subtraction_set[i]);
    }
    int query_size = p->hyperedge2node[query_edge].size();

    std::vector<int>temp;
    for(int i = 0; i < unnei.size(); ++i){
        temp.clear();
        unordered_sub_hash(candidate, hyperedgeAdj[unnei[i]][query_size], temp);
        candidate.swap(temp);
    }
    return;
}

void hypergraph::generate_hyperedge_candidates_unlable(Pattern *p, std::vector<int> &matching_order, int query_edge,
                                                       std::vector<int> &subtraction_set, std::vector<int> &candidate) {
    std::vector<int>nei;
    std::vector<int>unnei;
    for(int i = 0; i < subtraction_set.size(); ++i){
        if(p->hyperedge_matrix[query_edge][matching_order[i]]){
            nei.push_back(subtraction_set[i]);
        }
        else
            unnei.push_back(subtraction_set[i]);
    }
    int query_size = p->hyperedge2node[query_edge].size();
    if(nei.size() > 1){
        intersect_smid512(hyperedgeAdj[nei[0]][query_size], hyperedgeAdj[nei[1]][query_size], candidate);
    }
    else
        candidate = hyperedgeAdj[nei[0]][query_size];
    for(int i = 2; i < nei.size(); ++i){
         int pos = binarySerach(deg[nei[i]], query_size);
         if(pos == -1)
             continue;
         int l = csr[nei[i]][pos];
         int r = csr[nei[i]][pos + 1];
        std::vector<int>temp;
        temp = intersect_simd512(candidate, hyperedge_adj[nei[i]], l, r);
        candidate.swap(temp);
    }

    for(int i = 0; i < unnei.size(); ++i){
         int pos = binarySerach(deg[unnei[i]], query_size);
         if(pos == -1)
             continue;
         int l = csr[unnei[i]][pos];
         int r = csr[unnei[i]][pos + 1];
         std::vector<int>temp;
         temp = unordered_sub_hash(candidate, hyperedge_adj[unnei[i]], l, r);
         candidate.swap(temp);
    }
    return;
}




void hypergraph::Enumerate_overlap(std::vector<int> &matching_order, std::vector<int> &candidate, int depth,
                           Pattern *p, std::vector<int> &subtraction_set, int64_t &count) {
    if (depth == p->e_cnt - 1){
        for (int i = 0; i < candidate.size(); ++i) {
            subtraction_set.push_back(candidate[i]);
            if (isValidEmbedding_overlap(subtraction_set, matching_order, depth, p)){
                count++;
            }
                
            subtraction_set.pop_back(); 
        }
    }
    else { 
        for (int i = 0; i < candidate.size(); ++i) { 
            subtraction_set.push_back(candidate[i]);
            if (!isValidEmbedding_overlap(subtraction_set, matching_order, depth, p)) {
                subtraction_set.pop_back();
            } else {
                std::vector<int> new_candidate;
                generate_hyperedge_candidates_unlable(p, matching_order, matching_order[depth + 1], subtraction_set, new_candidate);
                Enumerate_overlap(matching_order, new_candidate, depth + 1, p, subtraction_set, count);
                subtraction_set.pop_back();
            }
        }
    }
    return;
}

void hypergraph::Enumerate_overlap_parallel_openmp(int thread_id, std::vector<int> &matching_order, std::vector<int> &candidate, int depth,
                           Pattern *p, std::vector<int> &subtraction_set,std::vector<std::vector<int>>& DAG_temp) {
    if (depth == p->e_cnt - 1){
        for (int i = 0; i < candidate.size(); ++i) {
            subtraction_set.push_back(candidate[i]);
            if (isValidEmbedding_overlap_parallel_old(subtraction_set, matching_order, depth, p,DAG_temp)){
                counts_overlap[thread_id*8]++;
            }
            subtraction_set.pop_back(); 
        }
    }
    else { 
        for (int i = 0; i < candidate.size(); ++i) { 
            subtraction_set.push_back(candidate[i]);
            if (!isValidEmbedding_overlap_parallel_old(subtraction_set, matching_order, depth, p,DAG_temp)) {
                subtraction_set.pop_back();
            } else {
                std::vector<int> new_candidate;
                generate_hyperedge_candidates_unlable(p, matching_order, matching_order[depth + 1], subtraction_set, new_candidate);
                Enumerate_overlap_parallel_openmp(thread_id,matching_order, new_candidate, depth + 1, p, subtraction_set,DAG_temp);
                subtraction_set.pop_back();
            }
        }
    }
    return;
}

void hypergraph::Enumerate_overlap_parallel(int thread_id,std::vector<int> &matching_order, std::vector<int> &candidate, int depth,
                           Pattern *p, std::vector<int> &subtraction_set,std::vector<std::vector<int>>& DAG_temp) {  
    if (depth == p->e_cnt - 1){
        for (int i = 0; i < candidate.size(); ++i) {
            subtraction_set.push_back(candidate[i]);
            if (isValidEmbedding_overlap_parallel(subtraction_set, matching_order, depth, p,DAG_temp))
                counts_overlap[thread_id*8]++;
            subtraction_set.pop_back(); 
        }
    }
    else { 
        for (int i = 0; i < candidate.size(); ++i) {
            subtraction_set.push_back(candidate[i]);
            if (!isValidEmbedding_overlap_parallel(subtraction_set, matching_order, depth, p,DAG_temp)) {
                subtraction_set.pop_back();
            } else {
                std::vector<int> new_candidate;
                new_candidate.clear();
                task_overlap *temp = new task_overlap(depth + 1,new_candidate,subtraction_set,DAG_temp);
                generate_hyperedge_candidates_unlable(p, matching_order, matching_order[depth + 1], subtraction_set, temp->candidate); //TODO:new_candidate粒度
                task_queues_overlap[thread_id].push(temp);
                // std::vector<int> new_candidate;
                // generate_hyperedge_candidates_unlable(p, matching_order, matching_order[depth + 1], subtraction_set, new_candidate);
                // Enumerate_overlap(int thread_id,matching_order, new_candidate, depth + 1, p, subtraction_set);
                subtraction_set.pop_back();
            }
        }
    }
    return;
}

int64_t hypergraph::overlap_aware_pattern_matching(Pattern *p) {
    DAG.resize(ROUND_UP(p->overlap_node));
    int64_t count = 0;
    int key = p->edge_to_lable[p->matching_order[0]];
    std::vector<int>&candidate = this->HyperGraph_table[key].hedge_list;
    std::vector<int>subtraction_set;
    subtraction_set.clear();
    Enumerate_overlap(p->matching_order, candidate, 0, p, subtraction_set, count);
    std::cout<<"count is "<<count<<std::endl;
    return count;
}

int64_t hypergraph::overlap_aware_pattern_matching_parallel(Pattern *p) {
    DAG.resize(ROUND_UP(p->overlap_node));
    int64_t count = 0;
    int key = p->edge_to_lable[p->matching_order[0]];
    std::vector<int>&initial_candidate = this->HyperGraph_table[key].hedge_list;
    
    std::vector<std::thread> threads;
    std::vector<int> subtraction_set;
    subtraction_set.clear();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::vector<int> indices(NUM_THREADS);
    std::iota(indices.begin(), indices.end(), 0); //用0-NUM_THREADS-1填充数组

    int num_elements = initial_candidate.size();
    int base_size = num_elements / NUM_THREADS;
    int extra = num_elements % NUM_THREADS;
    int start_index = 0;

    std::vector<std::vector<int>> threads_candidate;

    for (int i = 0; i < NUM_THREADS; ++i) {
        int end_index = start_index + base_size + (i < extra ? 1 : 0);
        std::vector<int> thread_candidate(initial_candidate.begin() + start_index, initial_candidate.begin() + end_index);
        threads_candidate.push_back(std::move(thread_candidate));
        start_index = end_index;
    }

    for (int i = 0; i < NUM_THREADS; ++i) {
        task_queues_overlap[i].push(new task_overlap(0,threads_candidate[i],subtraction_set,DAG)); //初始数据集的task，如果不push到队列中前面的线程会修改vector和map
    }
    
    auto computation_thread_main = [&](int thread_id) {
        WorkStealingQueue<task_overlap*>& task_queue_now = task_queues_overlap[thread_id];
        std::vector<WorkStealingQueue<task_overlap*>>& task_queues_now = task_queues_overlap; 
        while (true) {
            while (!task_queue_now.empty()) {
                auto temp_task_optional = task_queues_now[thread_id].pop();
                if(temp_task_optional.has_value()) {
                    auto temp_task=temp_task_optional.value();
                    Enumerate_overlap_parallel(thread_id,p->matching_order, temp_task->candidate, temp_task->depth, p, temp_task->subtraction_set,temp_task->DAG);
                    delete temp_task;
                }
                else break;
            }
            if(isending_overlap()) break;
            else{
                std::shuffle(indices.begin(), indices.end(), gen);//打乱数组，相当于random()%NUM_THREADS，但不重复
                for (int i : indices) { 
                // for (int i=thread_id+1;i!=thread_id;i=(i+1)%NUM_THREADS){
                    if(!task_queues_now[i].empty()){
                        int steal_size =task_queues_now[i].size()/2; //steal选定队列一半数据
                        for(int j=0; j<steal_size;j++){  
                            auto temp =task_queues_now[i].steal();
                            if(temp.has_value()) task_queue_now.push(temp.value());  
                            else break;
                        }

                        break; //steal成功后跳出for (int i : indices)
                    }
                }
            }
        }
        // std::cout << "Thread " << thread_id << " finished." << std::endl;
    };

    // 创建并启动线程
    for (int i = 0; i < NUM_THREADS; ++i) {
        threads.push_back(std::thread(computation_thread_main, i));
    }

    // 等待所有线程完成
    for (auto& t : threads) {
        t.join();
    }
    
    for(int i=0;i<NUM_THREADS;i++) {
        // std::cout << "Thread "<<i<<" count is " << counts[i] << std::endl;
        count+=counts_overlap[i*8];
    }
    std::cout << "Total count is " << count << std::endl;
    return count;
}

int64_t hypergraph::overlap_aware_pattern_matching_parallel_openmp(Pattern *p) {
    int64_t count = 0;
    int key = p->edge_to_lable[p->matching_order[0]];
    std::vector<int>&initial_candidate = this->HyperGraph_table[key].hedge_list;
    
#pragma omp parallel for schedule(dynamic) 
    for(int i=0;i<initial_candidate.size();i++) {
        int tid = omp_get_thread_num();

        std::vector<int> subtraction_set;
        subtraction_set.clear();
        subtraction_set.push_back(initial_candidate[i]);

        std::vector<std::vector<int>> DAG_t;
        DAG_t.resize(ROUND_UP(p->overlap_node));
        

        if (!isValidEmbedding_overlap_parallel_old(subtraction_set, p->matching_order, 0, p,DAG_t)) {
            subtraction_set.pop_back();
        } else {
            std::vector<int> new_candidate;
            generate_hyperedge_candidates_unlable(p, p->matching_order, p->matching_order[1], subtraction_set, new_candidate);
            Enumerate_overlap_parallel_openmp(tid,p->matching_order, new_candidate, 1, p, subtraction_set,DAG_t);
            subtraction_set.pop_back();
        }
    }

    for(int i=0;i<NUM_THREADS;i++) {
        count+=counts_overlap[i*8];
    }
    std::cout << "Total count is " << count << std::endl;
    return count;
}


void hypergraph::random_walk(int num_edges_to_sample, int min_vertex_count, int max_vertex_count, const char* pattern_path, const char* node_lable_path) {
    srand((unsigned)time(nullptr));
    std::vector<int> sampled_edges;
    std::set<int> sampled_nodes;

    while(true) {
        sampled_edges.clear();
        sampled_nodes.clear();
        int start_edge;
        start_edge = rand() % e_cnt;
        while (hyperedge2node[start_edge].size() > max_vertex_count / 2)
            start_edge = rand() % e_cnt;
        sampled_edges.push_back(start_edge);
        std::vector<int> &start_nodes_in_edge = hyperedge2node[start_edge];
        for (int node: start_nodes_in_edge)
            sampled_nodes.insert(node);
        std::vector<int> available_edges;
        for (auto edge : hyperedge_adj[start_edge])
            available_edges.emplace_back(edge);
        int attempts = 0;
        while (sampled_edges.size() < num_edges_to_sample) {
            attempts++;
            if(!available_edges.empty()) {
                int new_edge = available_edges[rand() % available_edges.size()];
                if(std::find(sampled_edges.begin(), sampled_edges.end(), new_edge) != sampled_edges.end())
                    continue;
                std::vector<int> &nodes_in_edge = hyperedge2node[new_edge];
                std::set<int> tmp_nodes(sampled_nodes);
                for (int node: nodes_in_edge) {
                    tmp_nodes.insert(node);
                }
                if (tmp_nodes.size() <= max_vertex_count) {
                    sampled_edges.emplace_back(new_edge);
                    for (int node: nodes_in_edge)
                        sampled_nodes.insert(node);
                    available_edges.clear();
                    for(auto e : hyperedge_adj[new_edge])
                        available_edges.emplace_back(e);
                }
            }
            if(attempts > 100)
                break;
        }
        if(sampled_nodes.size() >= min_vertex_count && sampled_edges.size() == num_edges_to_sample)
            break;
    }

    printf("node num is %d, edge num is %d\n", sampled_nodes.size(), sampled_edges.size());

    printf("sample edges is\n");
    for (auto val : sampled_edges){
        printf("%d ", val);
    }
    printf("\n");
    std::ofstream vertexFile(node_lable_path);
    std::ofstream edgeFile(pattern_path);

    std::map<int,int>mp;
    int node_num = 0;

    for (int node : sampled_nodes){
        vertexFile << this->v_to_lable[node] << "\n";
        mp[node] = node_num++;
    }
    for (int edge : sampled_edges) {
        edgeFile << hedge_to_kid[edge] << "\n";
        const std::vector<int>& nodes_in_edge = HyperGraph_table[hedge_to_kid[edge]].sedge[edge];
        for (int i = 0; i < nodes_in_edge.size() - 1; ++i) {
            edgeFile << mp[nodes_in_edge[i]] << ",";
        }
        edgeFile << mp[nodes_in_edge.back()];
        edgeFile << "\n";
    }
    vertexFile.close();
    edgeFile.close();
}

void hypergraph::random_walk_dense(int num_edges_to_sample, int min_vertex_count, int max_vertex_count, const char* pattern_path, const char* node_lable_path){
    srand((unsigned)time(nullptr));
    std::vector<int> sampled_edges;
    std::set<int> sampled_nodes;

    while(true) {
        sampled_edges.clear();
        sampled_nodes.clear();
        int start_edge;
        start_edge = rand() % e_cnt;
        while (hyperedge2node[start_edge].size() > max_vertex_count / 2)
            start_edge = rand() % e_cnt;
        sampled_edges.push_back(start_edge);
        std::vector<int> &start_nodes_in_edge = hyperedge2node[start_edge];
        for (int node: start_nodes_in_edge)
            sampled_nodes.insert(node);
        std::vector<int> available_edges;
        for (auto edge : hyperedge_adj[start_edge])
            available_edges.emplace_back(edge);
        int attempts = 0;
        while (sampled_edges.size() < num_edges_to_sample) {
            attempts++;
            if(!available_edges.empty()) {
                int new_edge = available_edges[rand() % available_edges.size()];
                if(std::find(sampled_edges.begin(), sampled_edges.end(), new_edge) != sampled_edges.end())
                    continue;
                std::vector<int> &nodes_in_edge = hyperedge2node[new_edge];
                std::set<int> tmp_nodes(sampled_nodes);
                for (int node: nodes_in_edge) {
                    tmp_nodes.insert(node);
                }
                if (tmp_nodes.size() <= max_vertex_count) {
                    sampled_edges.emplace_back(new_edge);
                    for (int node: nodes_in_edge)
                        sampled_nodes.insert(node);
                    available_edges.clear();
                    available_edges = hyperedge_adj[new_edge];
                    for(int i = 0; i < sampled_edges.size() - 1; ++i){
                        available_edges = intersect(available_edges, hyperedge_adj[sampled_edges[i]]);
                    }
                }
            }
            if(attempts > 100)
                break;
        }
        if(sampled_nodes.size() >= min_vertex_count && sampled_edges.size() == num_edges_to_sample)
            break;
    }

    printf("node num is %d, edge num is %d\n", sampled_nodes.size(), sampled_edges.size());

    printf("sample edges is\n");
    for (auto val : sampled_edges){
        printf("%d ", val);
    }
    printf("\n");
    std::ofstream vertexFile(node_lable_path);
    std::ofstream edgeFile(pattern_path);

    std::map<int,int>mp;
    int node_num = 0;

    for (int node : sampled_nodes){
        vertexFile << this->v_to_lable[node] << "\n";
        mp[node] = node_num++;
    }
    for (int edge : sampled_edges) {
        edgeFile << hedge_to_kid[edge] << "\n";
        const std::vector<int>& nodes_in_edge = HyperGraph_table[hedge_to_kid[edge]].sedge[edge];
        for (int i = 0; i < nodes_in_edge.size() - 1; ++i) {
            edgeFile << mp[nodes_in_edge[i]] << ",";
        }
        edgeFile << mp[nodes_in_edge.back()];
        edgeFile << "\n";
    }
    vertexFile.close();
    edgeFile.close();
}

