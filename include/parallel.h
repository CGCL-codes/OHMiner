#include "wsq.h"
#include "vector"
#include "unordered_map"

const int NUM_THREADS = 128;
int threshold = 0;
struct task {
    int depth; 
    std::vector<int> candidate; //TODO:use array
    std::vector<int> subtraction_set;
    std::unordered_map<int, std::vector<int>> num_vertices;

    // 默认构造函数
    task() : depth(0) {}
    
    // 参数化构造函数
    task(int d, const std::vector<int>& c, const std::vector<int>& s, const std::unordered_map<int, std::vector<int>>& nv)
        : depth(d), candidate(c), subtraction_set(s), num_vertices(nv) {}

    // 拷贝构造函数
    task(const task& other)
        : depth(other.depth), candidate(other.candidate), subtraction_set(other.subtraction_set), num_vertices(other.num_vertices) {}

    // 拷贝赋值操作符
    task& operator=(const task& other) {
        if (this == &other) {
            return *this;
        }
        depth = other.depth;
        candidate = other.candidate;
        subtraction_set = other.subtraction_set;
        num_vertices = other.num_vertices;

        return *this;
    }

    // 移动构造函数
    task(task&& other) noexcept
        : depth(other.depth), candidate(std::move(other.candidate)), subtraction_set(std::move(other.subtraction_set)), num_vertices(std::move(other.num_vertices)) {}

    // 移动赋值操作符
    task& operator=(task&& other) noexcept {
        if (this == &other) {
            return *this;
        }
        depth = other.depth;
        candidate = std::move(other.candidate);
        subtraction_set = std::move(other.subtraction_set);
        num_vertices = std::move(other.num_vertices);

        return *this;
    }
    ~task(){
        candidate.clear();
        subtraction_set.clear();
        num_vertices.clear();
    }

};

struct task_overlap {
    int depth; 
    std::vector<int> candidate; //TODO:use array
    std::vector<int> subtraction_set;
    std::vector<std::vector<int>> DAG;

    // 默认构造函数
    task_overlap() : depth(0) {}
    
    // 参数化构造函数
    task_overlap(int d, const std::vector<int>& c, const std::vector<int>& s,const std::vector<std::vector<int>>& D)
        : depth(d), candidate(c), subtraction_set(s),DAG(D) {}

    // 拷贝构造函数
    task_overlap(const task_overlap& other)
        : depth(other.depth), candidate(other.candidate), subtraction_set(other.subtraction_set),DAG(other.DAG) {}

    // 拷贝赋值操作符
    task_overlap& operator=(const task_overlap& other) {
        if (this == &other) {
            return *this;
        }
        depth = other.depth;
        candidate = other.candidate;
        subtraction_set = other.subtraction_set;
        DAG = other.DAG;
        return *this;
    }

    // 移动构造函数
    task_overlap(task_overlap&& other) noexcept
        : depth(other.depth), candidate(std::move(other.candidate)), subtraction_set(std::move(other.subtraction_set)),DAG(std::move(other.DAG)) {}

    // 移动赋值操作符
    task_overlap& operator=(task_overlap&& other) noexcept {
        if (this == &other) {
            return *this;
        }
        depth = other.depth;
        candidate = std::move(other.candidate);
        subtraction_set = std::move(other.subtraction_set);
        DAG = std::move(other.DAG);
        return *this;
    }
    ~task_overlap(){
        candidate.clear();
        subtraction_set.clear();
        DAG.clear();
    }

};

std::vector<WorkStealingQueue<task*>> task_queues(NUM_THREADS); 
std::vector<int64_t> counts(NUM_THREADS*8,0);

bool isending(){
    for(int i=0;i<NUM_THREADS;i++)
        if(!task_queues[i].empty()) return false;
    return true;
}

std::vector<WorkStealingQueue<task_overlap*>> task_queues_overlap(NUM_THREADS); 
std::vector<int64_t> counts_overlap(NUM_THREADS*8,0);


bool isending_overlap(){
    for(int i=0;i<NUM_THREADS;i++)
        if(!task_queues_overlap[i].empty()) return false;
    return true;
}
