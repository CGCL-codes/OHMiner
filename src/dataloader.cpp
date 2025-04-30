#include "dataloader.h"
#include "fstream"
#include <iostream>
#include <unordered_set>
#include "vector"
#include "algorithm"
#include <sstream>
#include "omp.h"
#include "common.h"

void Dataloader::load_DataGraph(hypergraph* &hg, const char *hyperedge_data, const char *node_lable_data) {
    std::ifstream nodeFile(node_lable_data);
    std::string line;
    std::unordered_map< std::string, std::vector<int> >label;
    int hyper_node = 0, vlabel_id = 0;
    while (getline(nodeFile, line))
    {
        label[line].push_back(hyper_node);
        hyper_node++;
    }
    hg->v_to_lable.clear();
    hg->v_to_lable.resize(hyper_node);
    for (auto tmp : label) {
        for(int v : tmp.second)
            hg->v_to_lable[v] = vlabel_id;
        vlabel_id++;
    }
    nodeFile.close();
    std::ifstream hyperFile(hyperedge_data);
    std::set<std::vector<int>> hedges;
    std::set<std::vector<int>> labels;
    std::map<std::vector<int>,int> key_to_kid;
    int origin_hyperedges = 0;           // hyperedge 的数量
    int num_hyperedges = 0;           // hyperedge 的数量
    int hlabel_id = 0;
    int v;
    while (getline(hyperFile, line)){
      origin_hyperedges++;
      for (int i = 0; i < line.size(); ++i) {
        if (line[i] == ',') {
          line[i] = ' ';
        }
      }

      std::stringstream sstream(line);
      std::vector<int> he;
      while (sstream >> v) {
        he.push_back(v-1);
      }
      //需要把有重复结点的超边给直接去掉
      std::sort(he.begin(), he.end());
      auto last = std::unique(he.begin(), he.end());
      he.erase(last, he.end());

      if (hedges.find(he) != hedges.end()) continue;

      hedges.insert(he);
      
      for (int i = 0; i < he.size(); ++i) {
        while(he[i] >= hg->node2hyperedge.size()){           //idx 是 vertex
            hg->node2hyperedge.push_back(std::vector<int>());
        }
        hg->node2hyperedge[he[i]].push_back(num_hyperedges);
      }
      std::vector<int> key;
      for (int i = 0; i < he.size(); ++i) {
        key.push_back(hg->v_to_lable[he[i]]);
      }
      std::sort(key.begin(), key.end());
      labels.insert(key);

      if (!key_to_kid.count(key)) {
        key_to_kid[key] = hlabel_id++;
      }
      
      hg->hedge_to_kid.push_back(key_to_kid[key]);
      hg->HyperGraph_table[key_to_kid[key]].sedge[num_hyperedges] = he;

      hg->HyperGraph_table[key_to_kid[key]].hedge_list.push_back(num_hyperedges);
      hg->hyperedge2node.push_back(he);
      num_hyperedges++;
    }
    printf("origin edge size = %d, edges size = %d, table size = %d\n", origin_hyperedges, hedges.size(), labels.size());
    hyperFile.close();

    for (auto& val : hg->HyperGraph_table) {
        for (auto& e : val.second.sedge){
            for (int j = 0; j < e.second.size(); ++j) {
                val.second.invert_index[e.second[j]].push_back(e.first);
            }
        }
    }

    hg->v_cnt = hyper_node;
    hg->e_cnt = num_hyperedges;

    


    hg->hyperedge_adj.resize(hg->e_cnt);
    hg->hyperedgeAdj.resize(hg->e_cnt);

    int num_threads = 128;

    //auto time1 = get_wall_time();
    std::vector< std::vector<long long> > upd_time(num_threads);
    for (int i = 0; i < num_threads; i++)
      upd_time[i].resize(hg->e_cnt, -1LL);

    #pragma omp parallel for schedule(dynamic) 
    for (int hyperedge_a = 0; hyperedge_a < hg->e_cnt; hyperedge_a++){
      int tid = omp_get_thread_num();
      int deg_a = 0;
      long long l_hyperedge_a = (long long)hyperedge_a;
      for (const int &node: hg->hyperedge2node[hyperedge_a]){
        for (const int &hyperedge_b: hg->node2hyperedge[node]){
          if (hyperedge_b == hyperedge_a) continue;
          if ((upd_time[tid][hyperedge_b] >> 31) ^ hyperedge_a){
            upd_time[tid][hyperedge_b] = (l_hyperedge_a << 31) + (long long)hg->hyperedge_adj[hyperedge_b].size();
            hg->hyperedge_adj[hyperedge_a].emplace_back(hyperedge_b);
            hg->hyperedgeAdj[hyperedge_a][hg->hyperedge2node[hyperedge_b].size()].push_back(hyperedge_b);
          }
          //hg->hyperedgeAdj[hyperedge_a][(int)(upd_time[tid][hyperedge_b] & 0x7FFFFFFFLL)];
        }
      }
    }

    #pragma omp parallel for schedule(dynamic) 
    for (int hyperedge_a = 0; hyperedge_a < hg->e_cnt; hyperedge_a++){
      auto & adjs = hg->hyperedgeAdj[hyperedge_a];
      for (auto &adj : adjs) {
        if (adj.second.size() > 1) {
          std::sort(adj.second.begin(), adj.second.end());
        }
      }
    }

    #pragma omp parallel for
    for(int hyperedge_a = 0; hyperedge_a < hg->e_cnt; hyperedge_a++){
      int64_t l_hyperedge_a = (int64_t)hyperedge_a;
      for(auto& val : hg->hyperedgeAdj[hyperedge_a]){
        for(auto v : val.second){
          int64_t key;
          if(v <= hyperedge_a)
            key = l_hyperedge_a << 32 | (int64_t)v;
          else
            key = (int64_t)v << 32 | l_hyperedge_a;
          if(hg->hyper_inter.find(key) != hg->hyper_inter.end())
            hg->hyper_inter[key] = intersect(hg->hyperedge2node[hyperedge_a], hg->hyperedge2node[v]);
        }
      }
    }

    hg->csr.resize(hg->e_cnt);
    hg->deg.resize(hg->e_cnt);

    #pragma omp parallel for
    for (int hyperedge_a = 0; hyperedge_a < hg->e_cnt; hyperedge_a++){
        std::sort(hg->hyperedge_adj[hyperedge_a].begin(), hg->hyperedge_adj[hyperedge_a].end(),[&](int& a, int& b){
          int sizea = hg->hyperedge2node[a].size(), sizeb = hg->hyperedge2node[b].size();
          return sizea != sizeb ? sizea < sizeb : a < b;
        });
        int last_degree = -1;
        for(int i = 0; i < hg->hyperedge_adj[hyperedge_a].size(); ++i){
          if(last_degree != hg->hyperedge2node[hg->hyperedge_adj[hyperedge_a][i]].size()){
            hg->deg[hyperedge_a].emplace_back(hg->hyperedge2node[hg->hyperedge_adj[hyperedge_a][i]].size());
            hg->csr[hyperedge_a].emplace_back(i);
          }
          last_degree = hg->hyperedge2node[hg->hyperedge_adj[hyperedge_a][i]].size();
        }
        hg->csr[hyperedge_a].emplace_back(hg->hyperedge_adj[hyperedge_a].size());
    }

    printf("stable_num:%d\nhyperedge_num:%d\nlable_num:%d\nnode_num:%d\n",hg->HyperGraph_table.size(),num_hyperedges,vlabel_id, hyper_node);
}

void Dataloader::load_PatternGraph(Pattern *&p, const char *pattern_data, const char *pattern_node_lable) {
    std::ifstream graphNodeFile(pattern_node_lable);
    std::string line;
    int hyper_node = 0;

    while (getline(graphNodeFile, line))
    {
        p->v_to_lable.push_back(std::stoi(line));        //相同的lable中的node
        hyper_node++;
    }

    p->v_cnt = hyper_node;
    graphNodeFile.close();

    std::ifstream graphFile(pattern_data);
    std::string lable;
    int v;
    p->e_cnt = 0;
    while (getline(graphFile, lable) && getline(graphFile, line)){
      for (int i = 0; i < line.size(); ++i) {
        if (line[i] == ',') {
          line[i] = ' ';
        }
      }

      std::stringstream sstream(line);
      std::vector<int> he;
      while (sstream >> v) {
        he.push_back(v);
      }
      
      // construct data structure
      for (int i = 0; i < he.size(); ++i) {
        while(he[i] >= p->node2hyperedge.size()){           //idx 是 vertex
            p->node2hyperedge.push_back(std::vector<int>());
        }
        p->node2hyperedge[he[i]].push_back(p->e_cnt);
      }
      p->edge_to_lable.push_back(std::stoi(lable));
      p->hyperedge2node.push_back(he);
      p->e_cnt++;
    }
    graphFile.close();

    p->hyperedge_adj.resize(p->e_cnt);
    std::vector<long long> upd_time(p->e_cnt, -1LL);

    //std::cout << "p->e_cnt number is "<< p->e_cnt << std::endl;
    p->hyperedge_matrix = std::vector< std::vector<int> >(p->e_cnt, std::vector<int>(p->e_cnt, 0));

    for (int hyperedge_a = 0; hyperedge_a < p->e_cnt; hyperedge_a++){
        long long l_hyperedge_a = (long long)hyperedge_a;
        for (const int &node: p->hyperedge2node[hyperedge_a]){
            for (const int &hyperedge_b: p->node2hyperedge[node]){
                if (hyperedge_b == hyperedge_a) continue;
                if ((upd_time[hyperedge_b] >> 31) ^ hyperedge_a){
                    upd_time[hyperedge_b] = (l_hyperedge_a << 31) + (long long)p->hyperedge_adj[hyperedge_b].size();
                    p->hyperedge_adj[hyperedge_b].push_back(hyperedge_a);
                    p->hyperedge_matrix[hyperedge_b][hyperedge_a] = 1;
                }
            }
        }
    }
}