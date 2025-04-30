#include "hypergraph.h"
#include "dataloader.h"

int main(int argc, char* argv[]){
    //超边保存的路径
    const char* hypergraph_path = argv[1];
    //对应的结点的lable保存的路径
    const char* node_lable = argv[2];
    hypergraph* hg = new hypergraph(0,0);
    
    Dataloader dataloader;
    dataloader.load_DataGraph(hg, hypergraph_path, node_lable);
    const char* pattern_path = argv[3];
    const char* pattern_node_lable = argv[4];

    Pattern *pattern = new Pattern();
    dataloader.load_PatternGraph(pattern, pattern_path, pattern_node_lable);

    pattern->ComputeMatchingOrder(hg);

    pattern->generate_overlap_graph_sub();
    pattern->generate_execution_plan();

    auto time1 = get_wall_time();
    hg->overlap_aware_pattern_matching(pattern);
    auto time2 = get_wall_time();

    printf("\ntime is %.8f\n", time2 - time1);    
    fflush(stdout);
     
    delete pattern;
    delete hg;
}