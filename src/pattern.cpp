#include "pattern.h"
#include "hypergraph.h"
#include "unordered_set"
#include "queue"
#include "iostream"
#include <limits.h>
#include "algorithm"



void Pattern::intersect(const std::vector<int>&set1, const std::vector<int>&set2, std::vector<int>&result){
    auto size1 = set1.size();
    auto size2 = set2.size();
    int i = 0, j = 0;
    while(i < size1 && j < size2){
        if(set1[i] == set2[j]){
            result.emplace_back(set1[i]);
            i++;
            j++;
        }
        else if(set1[i] > set2[j])
            j++;
        else
            i++;
    }
    //这个是有用的
    if(set1.size() == set2.size() && result.size() == set1.size())
        result.clear();
    return;
}

int Pattern::get_cardinality(hypergraph *g, int hyperedge) {
    int lable = edge_to_lable[hyperedge];
    int count = 0;
    for (int i = 0; i < g->e_cnt; ++i){
        if(lable == g->hedge_to_kid[i])
            count++;
    }
    count = g->HyperGraph_table[lable].sedge.size();
    return count;
}

int Pattern::intersection(std::set<int> v_list, int hedge) {
    std::vector<int>node_list = hyperedge2node[hedge];
    std::set<int>::iterator i = v_list.begin();
    std::vector<int>::iterator j = node_list.begin();
    int count = 0;
    while(i != v_list.end() && j != node_list.end()){
        if(*i < *j){
            i++;
        }
        else if(*i > *j){
            j++;
        }
        else{
            i++;
            j++;
            count++;
        }
    }
    return count;
}


void Pattern::ComputeMatchingOrder(hypergraph *g) {
    bool used[this->e_cnt];
    std::fill(used, used + this->e_cnt, false);
    int min_cardinality = INT_MAX, e1 = -1;
    for (int i = 0; i < e_cnt; ++i) {
        int cardinality = get_cardinality(g, i);
        if (min_cardinality > cardinality) {
            min_cardinality = cardinality;
            e1 = i;
        }
    }
    matching_order.push_back(e1);
    used[e1] = true;
    std::set<int> v_list;
    while (matching_order.size() != e_cnt) {
        for (int i = 0; i < hyperedge2node[matching_order.back()].size(); ++i)
            v_list.insert(hyperedge2node[matching_order.back()][i]); 
        int min_cardinality = INT_MAX, e = -1;
        for (int i = 0; i < e_cnt; ++i) {
            if (!used[i]) {
                int cardinality = get_cardinality(g, i);
                double intersect_num = intersection(v_list, i);
                if( !intersect_num )
                    intersect_num = 1e-3;
                if (min_cardinality > cardinality / intersect_num) {
                    min_cardinality = cardinality / intersect_num;
                    e = i;
                }
            }
        }
        used[e] = true;
        matching_order.push_back(e);
    }


    this->rank.resize(e_cnt);

    for(int i = 0; i < e_cnt; ++i) rank[matching_order[i]] = i;

    for(auto e : matching_order){
        for(auto v : hyperedge2node[e])
            pattern_profile[v].push_back(rank[e]);
        vertives_num.push_back(pattern_profile.size());
    }

    printf("\n");
    for(int i = 0; i < e_cnt; ++i) printf("%d ", matching_order[i]);
    printf("\n");
    return;
}



void Pattern::generate_overlap_graph_sub() {
    std::vector<std::pair<int, std::vector<int>>>HE;
    overlap_graph.resize(e_cnt);
    for(auto he : matching_order){
        overlap_graph[he] = hyperedge2node[he];
        overlap_graph_adj.emplace_back(std::vector<int>());
        overlap_graph_adj_in.emplace_back(std::vector<int>());
        HE.emplace_back(std::make_pair(he, std::vector<int>()));
    }
    expand_overlap_graph_merge(HE, 0);

    // overlap_graph_adj_in.resize(overlap_graph_adj.size());
    // for(int i = 0; i < overlap_graph_adj.size(); ++i){
    //     for(int j = 0; j < overlap_graph_adj[i].size(); ++j){
    //         overlap_graph_adj_in[overlap_graph_adj[i][j]].emplace_back(i);
    //     }
    // }
    overlap_node = overlap_graph.size();
    return;
}

void Pattern::generate_overlap_graph() {
    std::vector<int>HE;
    overlap_graph.resize(e_cnt);
    for(auto he : matching_order){
        overlap_graph[he] = hyperedge2node[he];
        overlap_graph_adj.emplace_back(std::vector<int>());
        HE.emplace_back(he);
    }
    expand_overlap_graph(HE);

    overlap_graph_adj_in.resize(overlap_graph_adj.size());
    for(int i = 0; i < overlap_graph_adj.size(); ++i){
        for(int j = 0; j < overlap_graph_adj[i].size(); ++j){
            overlap_graph_adj_in[overlap_graph_adj[i][j]].emplace_back(i);
        }
    }
    overlap_node = overlap_graph.size();
    return;
}

void Pattern::expand_overlap_graph(const std::vector<int>& HE) {
    std::vector<int>new_HE;
    if(!HE.empty()){
        auto size = HE.size();
        std::map<std::vector<int>, std::vector<int>>mp;
        for(int i = 0; i < size; ++i)
            for(int j = i + 1; j < size; ++j){
                std::vector<int>overlap;
                intersect(overlap_graph[HE[i]], overlap_graph[HE[j]], overlap);
                if(!overlap.empty()){
                    overlap_graph_adj.emplace_back(std::vector<int>());
                    overlap_graph.emplace_back(overlap);
                    auto node_num = overlap_graph.size() - 1;
                    overlap_graph_adj[HE[i]].emplace_back(node_num);
                    overlap_graph_adj[HE[j]].emplace_back(node_num);
                    mp[overlap].emplace_back(node_num);
                }
            }
        for(auto kv : mp){
            if(kv.second.size() > 1){
                overlap_graph_adj.emplace_back(std::vector<int>());
                overlap_graph.emplace_back(kv.first);
                auto node_num = overlap_graph.size() - 1;
                for(auto v : kv.second){
                    overlap_graph_adj[v].push_back(node_num);
                }
                new_HE.emplace_back(node_num);
            }
            else{
                new_HE.emplace_back(kv.second.front());
            }
        }
        expand_overlap_graph(new_HE);
    }
    return;
}




void Pattern::expand_overlap_graph_merge(const std::vector<std::pair<int, std::vector<int>>>& HE, int depth) {
    //定义一个new_HE,存储下一层节点，其中每一个pair的first是 节点的编号， second是结点的parent, 如果是merge了节点，则是结点的父亲的父亲
    std::vector<std::pair<int, std::vector<int>>> new_HE;
    if(!HE.empty()){
        auto size = HE.size();
        std::map<std::vector<int>, std::vector<int>>mp;
        for(int i = 0; i < size; ++i)
            for(int j = i + 1; j < size; ++j){
                std::vector<int>overlap;
                intersect(overlap_graph[HE[i].first], overlap_graph[HE[j].first], overlap);
                if(!overlap.empty()){
                    //逻辑：先增加一个节点，然后在对应的overlap_graph_adj和overlap_graph_adj_in中新增一个节点，然后是一些添加的逻辑
                    overlap_graph.emplace_back(overlap);
                    overlap_graph_adj.emplace_back(std::vector<int>());
                    overlap_graph_adj_in.emplace_back(std::vector<int>());
                    auto node_num = overlap_graph.size() - 1;
                    overlap_graph_adj[HE[i].first].emplace_back(node_num);
                    overlap_graph_adj[HE[j].first].emplace_back(node_num);
                    overlap_graph_adj_in[node_num].push_back(HE[i].first);
                    overlap_graph_adj_in[node_num].push_back(HE[j].first);
                    matrix[HE[i].first][HE[j].first] = 1;
                    matrix[HE[j].first][HE[i].first] = 1;
                    mp[overlap].emplace_back(node_num);
                }
                else{
                    //获取祖先节点
                    std::vector<int>nodes;
                    for(auto v : HE[i].second){
                        nodes.emplace_back(v);
                    }
                    for(auto v : HE[j].second){
                        nodes.emplace_back(v);
                    }
                    std::sort(nodes.begin(), nodes.end());
                    auto last = std::unique(nodes.begin(), nodes.end());
                    nodes.erase(last, nodes.end());


                    if(depth != 0 && is_tuan(nodes)){
                        overlap_graph.emplace_back(overlap);
                        //还是需要成对出现
                        overlap_graph_adj.emplace_back(std::vector<int>());
                        overlap_graph_adj_in.emplace_back(std::vector<int>());
                        auto node_num = overlap_graph.size() - 1;
                        overlap_graph_adj[HE[i].first].emplace_back(node_num);
                        overlap_graph_adj[HE[j].first].emplace_back(node_num);
                        overlap_graph_adj_in[node_num].emplace_back(HE[i].first);
                        overlap_graph_adj_in[node_num].emplace_back(HE[j].first);
                    }
                }
            }
        for(auto kv : mp){
            std::pair<int, std::vector<int>>new_node;
            if(kv.second.size() > 1){
                overlap_graph.emplace_back(kv.first);
                overlap_graph_adj.emplace_back(std::vector<int>());
                overlap_graph_adj_in.emplace_back(std::vector<int>());
                auto node_num = overlap_graph.size() - 1;
                new_node.first = node_num;
                for(auto v : kv.second){
                    overlap_graph_adj[v].emplace_back(node_num);
                    overlap_graph_adj_in[node_num].emplace_back(v);
                    for(int parent : overlap_graph_adj_in[v]){
                        new_node.second.emplace_back(parent);
                    }
                } 
            }
            else{
                new_node.first = kv.second.front();
                for(int parent : overlap_graph_adj_in[new_node.first]){
                    new_node.second.emplace_back(parent);
                }
            }
            new_HE.emplace_back(new_node);
        }
        expand_overlap_graph_merge(new_HE, depth + 1);
    }
    return;
}

void Pattern::generate_execution_plan() {
    std::queue<int>q;
    bool visited[overlap_node];
    std::fill(visited, visited + overlap_node, false);
    EP.resize(e_cnt);
    root_leaf.resize(e_cnt);
    for(int i = 0; i < matching_order.size(); ++i){
        q.push(matching_order[i]);
        while(!q.empty()){
            int v = q.front();
            q.pop();
            EP[matching_order[i]].emplace_back(v);
            if(overlap_graph_adj[v].empty()){
                root_leaf[matching_order[i]].emplace_back(v);
            }
            visited[v] = true;
            for(auto out : overlap_graph_adj[v]){
                bool flag = true;
                for(auto in : overlap_graph_adj_in[out]){
                    if(!visited[in]){
                        flag = false;
                        break;
                    }
                }
                if(flag)
                    q.push(out);
            }
        }
    }
    return;
}

bool Pattern::is_tuan(std::vector<int>& nodes){
    for(int i = 0; i < nodes.size(); ++i){
        for(int j = i + 1; j < nodes.size(); ++j){
            if(matrix[nodes[i]][nodes[j]] == 0)
                return false;
        }
    }
    return true;
}

void Pattern::print() {
    for(int i = 0; i < overlap_node; ++i){
        std::cout << "node num is " << i << std::endl;
        for(auto v : overlap_graph_adj[i])
            std::cout << v << " ";
        std::cout << std::endl;
        for(auto v : overlap_graph[i])
            std::cout << v << " ";
        std::cout << std::endl;
        for(auto v : overlap_graph_adj_in[i])
            std::cout << v << " ";
        std::cout << std::endl;
    }

    std::cout << "excucate plan:" << std::endl;
    for(auto node : EP){
        for(auto v : node)
            std::cout << v << " ";
        std::cout << std::endl;
    }
}