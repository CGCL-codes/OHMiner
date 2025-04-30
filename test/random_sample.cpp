#include "hypergraph.h"
#include "dataloader.h"
#include "iostream"
#include "common.h"
int main(int argc, char* argv[]){
    //超边保存的路径
    const char* hypergraph_path = argv[1];
    //对应的结点的lable保存的路径
    const char* node_lable = argv[2];
    hypergraph* hg = new hypergraph(0,0);
    
    Dataloader dataloader;

    auto time1 = get_wall_time();
    dataloader.load_DataGraph(hg, hypergraph_path, node_lable);
    auto time2 = get_wall_time();

    std::cout << time2 - time1 << std::endl;


    const char* pattern_path = argv[3];
    const char* pattern_node_lable = argv[4];

    std::cout << pattern_path << std::endl;
    std::cout << pattern_node_lable << std::endl;

    int hyper_edge_num = atoi(argv[5]);
    int min_node_num = atoi(argv[6]);
    int max_node_num = atoi(argv[7]);
    hg->random_walk_dense(hyper_edge_num, min_node_num, max_node_num, pattern_path, pattern_node_lable);
    delete hg;
}