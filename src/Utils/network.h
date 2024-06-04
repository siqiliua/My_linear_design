
#include <map>
#include <unordered_map>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <memory>
#include <string>
#include <limits>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <array>

#include "utility_v.h"
#include "common.h"
#include "codon.h"

using namespace std;

// #define is_verbose

namespace LinearDesign {

template <typename IndexType,
          typename IndexWType = tuple<IndexType, double>,
          typename NodeType = tuple<IndexType, NumType, NucType>,
          typename NodeWType = tuple<NodeType, double>>
class Lattice {
public:
    unordered_map<IndexType, vector<NodeType>> nodes;
    unordered_map<NodeType, vector<NodeWType>, std::hash<NodeType>> left_edges;
    unordered_map<NodeType, vector<NodeWType>, std::hash<NodeType>> right_edges;
    
    Lattice(): nodes(), left_edges(), right_edges() {};

    void add_edge(NodeType n1, NodeType n2, double weight = 0.0f){
        right_edges[n1].push_back(make_pair(n2, weight));
        left_edges[n2].push_back(make_pair(n1, weight));
    }

    void add_node(NodeType n1){
        IndexType pos = get<0>(n1);
        nodes[pos].push_back(n1);

    }
};

template <typename IndexType,
          typename IndexWType = tuple<IndexType, double>,
          typename NodeType = tuple<IndexType, NumType, NucType>,
          typename NodeWType = pair<NodeType, double>>
class DFA {
public:
    unordered_map<IndexType, vector<NodeType>> nodes;
    unordered_map<NodeType, vector<NodeWType>, std::hash<NodeType>> left_edges;
    unordered_map<NodeType, vector<NodeWType>, std::hash<NodeType>> right_edges;
    unordered_map<NodeType, unordered_map<NodeType, double, std::hash<NodeType>>, std::hash<NodeType>> auxiliary_left_edges;
    unordered_map<NodeType, unordered_map<NodeType, double, std::hash<NodeType>>, std::hash<NodeType>> auxiliary_right_edges;
    unordered_map<NodeType, double, std::hash<NodeType>> node_rightedge_weights;

    DFA(): nodes(), left_edges(), right_edges(), auxiliary_left_edges(), auxiliary_right_edges() {};

    void add_edge(NodeType n1, NodeType n2, double weight = 0.0f){
        right_edges[n1].push_back(make_pair(n2, weight));
        left_edges[n2].push_back(make_pair(n1, weight));
        //auxiliary_right_edges[n1][n2].push_back(weight);
        //auxiliary_left_edges[n2][n1].push_back(weight);
        auxiliary_right_edges[n1][n2]= weight;
        auxiliary_left_edges[n2][n1]= weight;
        node_rightedge_weights[n1]=weight;
    }

    void add_node(NodeType n1){
        IndexType pos = get<0>(n1);
        nodes[pos].push_back(n1);
    }    
};


template <typename IndexType,
          typename NodeType = tuple<IndexType, NumType, NucType>,
          typename LatticeType = Lattice<IndexType>>
unordered_map<string, LatticeType> read_wheel(string const &filename) {
    unordered_map<string, LatticeType> aa_graphs;
    ifstream inFile;
    inFile.open(filename);
    if (!inFile) {
        printf("Unable to open coding_wheel file\n");
        exit(1);   // call system to stop
    }

    vector<string> stuff;
    vector<string> option_splited;
    string aa;
    IndexType i;
    IndexType j;

    for (string line; getline(inFile, line);) {
        stuff = split(line, '\t');

        aa = stuff[0];
        LatticeType graph = LatticeType();
        graph.add_node(make_tuple(0,0,0)); // always initialize with node (0,0)

        char last_first = 0;
        vector<string>::iterator iter = stuff.begin();
        ++iter; // position 0 is aa name
        i = 0;
        j = 0;
        while(iter != stuff.end()){
            string option = *iter;
            option_splited = split(option, ' ');
            char first = option_splited[0][0];
            char second = option_splited[1][0];
            string thirds = option_splited[2];
            NodeType n2 = make_tuple(2, i, GET_ACGU_NUC(second));
            graph.add_node(n2);
            NodeType n1;
            if (first != last_first) {
                n1 = make_tuple(1, i, GET_ACGU_NUC(first));
                graph.add_node(n1);
                graph.add_edge(make_tuple(0, 0, 0), n1);
            }
            else {
                n1 = make_tuple(1, i-1, GET_ACGU_NUC(first));
            }
            last_first = first;
            graph.add_edge(n1, n2);
            for (auto& third : thirds) {
                NodeType n3;
                n3 = make_tuple(3, j, GET_ACGU_NUC(third));
                graph.add_node(n3);
                graph.add_edge(n2, n3);
                graph.add_edge(n3, make_tuple(0, 0, 0));
                j++;
            }
            i++; iter++;
        }
        aa_graphs[aa] = graph;
#ifdef is_verbose
        printf("-----------------Lattice------------------------\n");
        for(IndexType pos = 0; pos <= 2; pos++){
            for(auto &node : graph.nodes[pos]){
                IndexType p = get<0>(node);
                IndexType num = get<1>(node);
                printf("node, (%d, %d)\n", p, num);
                for(auto &item : graph.right_edges[node]){
                    NodeType n2 = get<0>(item);
                    IndexType p2 = get<0>(n2); IndexType num2 = get<1>(n2);
                    IndexType nuc = get<1>(item);
                    double weight = get<2>(item);
                    printf("              (%d, %d) -(%d,%lf)-> (%d, %d)\n", p, num, nuc,  weight, p2, num2);
                }
                for(auto &item : graph.left_edges[node]){
                    NodeType n1 = get<0>(item);
                    IndexType p1 = get<0>(n1); IndexType num1 = get<1>(n1);
                    IndexType nuc = get<1>(item);
                    double weight = get<2>(item);
                    printf("  (%d, %d) <-(%d,%lf)- (%d, %d)\n", p1, num1, nuc, weight, p, num);
                }
            }
        }        
#endif
    }
    inFile.close();
    return aa_graphs;
}


template <typename IndexType,
          typename NodeType = tuple<IndexType, NumType, NucType>,
          typename LatticeType = Lattice<IndexType>,
          typename NucType = IndexType,
          typename NodeNucNodeType = std::tuple<NodeType, NucType, NodeType>>
unordered_map<string, LatticeType> read_wheel_with_weights(const std::string& filename,
        std::unordered_map<std::string, std::unordered_map<NodeType, double, std::hash<NodeType>>>& nodes_with_best_weight,
        std::unordered_map<std::string, std::unordered_map<NodeNucNodeType, double, std::hash<NodeNucNodeType>>>& edges_with_best_weight,
        const Codon& codon) {
    unordered_map<string, LatticeType> aa_graphs;
    ifstream inFile;
    inFile.open(filename);
    if (!inFile) 
        throw std::runtime_error("Unable to open coding_wheel file\n");

    vector<string> stuff;
    vector<string> option_splited;
    string aa;
    IndexType i;
    IndexType j;

    for (string line; getline(inFile, line);) {
        stuff = split(line, '\t');
        aa = stuff[0];
        LatticeType graph = LatticeType();
        graph.add_node(make_tuple(0,0,0)); // always initialize with node (0,0)

        char last_first = 0;
        vector<string>::iterator iter = stuff.begin();
        ++iter; // position 0 is aa name
        i = 0;
        j = 0;
        while(iter != stuff.end()){
            string option = *iter;
            option_splited = split(option, ' ');
            char first = option_splited[0][0];
            char second = option_splited[1][0];
            string thirds = option_splited[2];

            NodeType n2 = make_tuple(2, i, GET_ACGU_NUC(second));
            graph.add_node(n2);    

            NodeType n1;
            if (first != last_first) {
                n1 = make_tuple(1, i, GET_ACGU_NUC(first));
                graph.add_node(n1);
                double weight = 0.0f;
                graph.add_edge(make_tuple(0, 0, 0), n1, weight);
            }
            else {
                n1 = make_tuple(1, i-1, GET_ACGU_NUC(first));
            }
            
            last_first = first;

            double weight = 0.0f;
            if (nodes_with_best_weight[aa].count(make_tuple(0, 0, 0))) {
                weight = edges_with_best_weight[aa][make_tuple(n1, n2)] / nodes_with_best_weight[aa][make_tuple(0, 0, 0)];
            }
            graph.add_edge(n1, n2, weight);

            for (auto& third : thirds) {
                NodeType n3;
                n3 = make_tuple(3, j, GET_ACGU_NUC(third));
                graph.add_node(n3);
                double weight2 = 0.0f;
                if (nodes_with_best_weight[aa].count(n1)) {
                    weight2 = edges_with_best_weight[aa][make_tuple(n2, n3)] / nodes_with_best_weight[aa][n1];
                }
                graph.add_edge(n2, n3, weight2);

                std::string three_nums = std::string(1, first) + std::string(1, second) + std::string(1, third);
                double weight = 0.0f;
                if (nodes_with_best_weight[aa].count(n2)) {
                    weight = codon.get_weight(aa, three_nums) / nodes_with_best_weight[aa][n2];
                } else {
                    weight = codon.get_weight(aa, three_nums);
                }
                graph.add_edge(n3, make_tuple(0, 0, 0), weight);
                j++;
            }
            i++; iter++;
        }
        aa_graphs[aa] = graph;
    }

    inFile.close();
    return aa_graphs;
}


template <typename IndexType,
          typename NodeType = tuple<IndexType, NumType, NucType>,
          typename LatticeType = Lattice<IndexType>,
          typename NucType = IndexType,
          typename NodeNucNodeType = std::tuple<NodeType, NucType, NodeType>>
unordered_map<string, LatticeType> read_wheel_with_weights_log(const std::string& filename,
        std::unordered_map<std::string, std::unordered_map<NodeType, double, std::hash<NodeType>>>& nodes_with_best_weight,
        std::unordered_map<std::string, std::unordered_map<NodeNucNodeType, double, std::hash<NodeNucNodeType>>>& edges_with_best_weight,
        const Codon& codon, double lambda_) {
    unordered_map<string, LatticeType> aa_graphs;
    ifstream inFile;
    inFile.open(filename);
    if (!inFile) 
        throw std::runtime_error("Unable to open coding_wheel file\n");

    vector<string> stuff;
    vector<string> option_splited;
    string aa;
    IndexType i;
    IndexType j;

    for (string line; getline(inFile, line);) {
        stuff = split(line, '\t');
        aa = stuff[0];
        LatticeType graph = LatticeType();
        graph.add_node(make_tuple(0,0,0)); // always initialize with node (0,0)

        char last_first = 0;
        vector<string>::iterator iter = stuff.begin();
        ++iter; // position 0 is aa name
        i = 0;
        j = 0;
        while(iter != stuff.end()){
            string option = *iter;
            option_splited = split(option, ' ');
            char first = option_splited[0][0];
            char second = option_splited[1][0];
            string thirds = option_splited[2];
            NodeType n2 = make_tuple(2, i, GET_ACGU_NUC(second));
            graph.add_node(n2);
            NodeType n1;
            if (first != last_first) {
                n1 = make_tuple(1, i, GET_ACGU_NUC(first));
                graph.add_node(n1);
                //auto first_num = GET_ACGU_NUC(first);
                double weight = 1.0f;
                graph.add_edge(make_tuple(0, 0, 0), n1, lambda_ * log(weight));
            }
            else {
                n1 = make_tuple(1, i-1, GET_ACGU_NUC(first));
            }
            
            last_first = first;
            double weight = 1.0f;
            if (nodes_with_best_weight[aa].count(make_tuple(0, 0, 0))) {
                weight = lambda_ * log(edges_with_best_weight[aa][make_tuple(n1, n2)] / nodes_with_best_weight[aa][make_tuple(0, 0, 0)]);
            }
            graph.add_edge(n1, n2, weight);

            for (auto& third : thirds) {
                NodeType n3;
                n3 = make_tuple(3, j, GET_ACGU_NUC(third));
                graph.add_node(n3);
                double weight2 = 1.0f;
                if (nodes_with_best_weight[aa].count(n1)) {
                    weight2 = lambda_ *  log(edges_with_best_weight[aa][make_tuple(n2, n3)] / nodes_with_best_weight[aa][n1]);
                }
                graph.add_edge(n2, n3, weight2);

                std::string three_nums = std::string(1, first) + std::string(1, second) + std::string(1, third);
                double weight = 1.0f;
                if (nodes_with_best_weight[aa].count(n2)) {
                    weight = lambda_ *  log(codon.get_weight(aa, three_nums) / nodes_with_best_weight[aa][n2]);
                } else {
                    weight = lambda_ *  log(codon.get_weight(aa, three_nums));
                }
                graph.add_edge(n3, make_tuple(0, 0, 0), weight);
                j++;
            }
            i++; iter++;
        }
        aa_graphs[aa] = graph;
    }
    inFile.close();
    return aa_graphs;
}

template <typename IndexType,
          typename NucType = IndexType,
          typename NodeType = tuple<IndexType, NumType, NucType>,
          typename NodeToNodeType = std::tuple<NodeType, NodeType>,
          typename WeightType = double,
          typename LatticeType = Lattice<IndexType>>
void prepare_codon_unit_lattice(const std::string& wheel_path, const Codon& codon,
        std::unordered_map<string, LatticeType>& aa_graphs_with_ln_weights_ret,
        std::unordered_map<std::string, std::unordered_map<std::tuple<NodeType, NodeType>, std::tuple<double, NucType, NucType, NucType>, std::hash<std::tuple<NodeType, NodeType>>>>&
                best_path_in_one_codon_unit_ret,
        std::unordered_map<std::string, std::string>& aa_best_path_in_a_whole_codon_ret, double lambda_) {

    std::unordered_map<std::string, std::unordered_map<NodeType, WeightType, std::hash<NodeType>>> nodes_with_best_weight;
    std::unordered_map<std::string, std::unordered_map<NodeToNodeType, WeightType, std::hash<NodeToNodeType>>> edges_with_best_weight;

    std::unordered_map<std::string, LatticeType> aa_graphs_with_ln_weights;
    std::unordered_map<std::string, LatticeType> aa_graphs_with_weights = read_wheel_with_weights<IndexType>(wheel_path, nodes_with_best_weight, edges_with_best_weight, codon);
    
    for (auto& aa_aa_elem : aa_graphs_with_weights) {
        auto& aa = aa_aa_elem.first;
        auto& aa_elem = aa_aa_elem.second;
        for (auto& node_at_3 : aa_elem.nodes[3]) {
            for (auto& node_at_4_weight : aa_elem.right_edges[node_at_3]) {
                auto node_at_4 = std::get<0>(node_at_4_weight);
                auto weight = std::get<1>(node_at_4_weight);
                nodes_with_best_weight[aa][node_at_3] = weight;
                edges_with_best_weight[aa][make_tuple(node_at_3,node_at_4)] = weight;
            }
        }

        for (auto& node_at_2 : aa_elem.nodes[2]) {
            for (auto& node_at_3_nuc_weight : aa_elem.right_edges[node_at_2]) {
                auto node_at_3 = std::get<0>(node_at_3_nuc_weight);
                nodes_with_best_weight[aa][node_at_2] = max(nodes_with_best_weight[aa][node_at_2], nodes_with_best_weight[aa][node_at_3]);
                //edges_with_best_weight[aa][make_tuple(node_at_2, node_at_3)] = max(nodes_with_best_weight[aa][node_at_2], nodes_with_best_weight[aa][node_at_3]);
            }
            for (auto& node_at_3_nuc_weight : aa_elem.right_edges[node_at_2]) {
                auto node_at_3 = std::get<0>(node_at_3_nuc_weight);
                //nodes_with_best_weight[aa][node_at_2] = max(nodes_with_best_weight[aa][node_at_2], nodes_with_best_weight[aa][node_at_3]);
                edges_with_best_weight[aa][make_tuple(node_at_2, node_at_3)] = nodes_with_best_weight[aa][node_at_2];
            }
        }

        for (auto& node_at_1 : aa_elem.nodes[1]) {
            for (auto& node_at_2_nuc_weight : aa_elem.right_edges[node_at_1]) {
                auto node_at_2 = std::get<0>(node_at_2_nuc_weight);
                nodes_with_best_weight[aa][node_at_1] = max(nodes_with_best_weight[aa][node_at_1], nodes_with_best_weight[aa][node_at_2]);
                }
            for (auto& node_at_2_nuc_weight : aa_elem.right_edges[node_at_1]) {
                auto node_at_2 = std::get<0>(node_at_2_nuc_weight);
                edges_with_best_weight[aa][make_tuple(node_at_1,node_at_2)] = nodes_with_best_weight[aa][node_at_1];
            }
            
        }

        for (auto& node_at_0 : aa_elem.nodes[0]) {
            for (auto& node_at_1_nuc_weight : aa_elem.right_edges[node_at_0]) {
                auto node_at_1 = std::get<0>(node_at_1_nuc_weight);
                nodes_with_best_weight[aa][node_at_0] = max(nodes_with_best_weight[aa][node_at_0], nodes_with_best_weight[aa][node_at_1]);
                }
            for (auto& node_at_1_nuc_weight : aa_elem.right_edges[node_at_0]) {
                auto node_at_1 = std::get<0>(node_at_1_nuc_weight);
                edges_with_best_weight[aa][make_tuple(node_at_0,node_at_1)] = nodes_with_best_weight[aa][node_at_0] ;
            }
        }
    }

    aa_graphs_with_ln_weights = read_wheel_with_weights_log<IndexType>(wheel_path,  nodes_with_best_weight, edges_with_best_weight, codon, lambda_);

    // for(int k=0;k<4;k++){
    //     for(auto &node : aa_graphs_with_ln_weights["Thr"].nodes[k]){
    //         IndexType index1;
    //         NumType num1;
    //         NucType nuc1;
    //         std::tie(index1,num1,nuc1) = node;
    //         for (auto &node_weight : aa_graphs_with_ln_weights["Thr"].right_edges[node]){
    //             IndexType index_nr;
    //             NumType num_nr;
    //             NucType nuc_nr;
    //             NodeType node_right = std::get<0>(node_weight);
    //             double weight = std::get<1>(node_weight);
    //             std::tie(index_nr,num_nr,nuc_nr) = node_right;
    //             std::cout << index1 << "," << num1 << "," << GET_ACGU(nuc1) << "," << weight << std::endl;
    //             std::cout << index_nr << "," << num_nr << "," << GET_ACGU(nuc_nr)  << std::endl;               
    //         }
    //         std::cout <<  "........."  << std::endl;
    //     }
    // }

    std::unordered_map<std::string, 
                       std::unordered_map<std::tuple<NodeType, NodeType>, 
                       std::tuple<double, NucType, NucType,NucType>,
                       std::hash<std::tuple<NodeType, NodeType>>>>
                       best_path_in_one_codon_unit;

    for (auto& aa_graph : aa_graphs_with_ln_weights) {
        auto& aa = aa_graph.first;
        auto& graph = aa_graph.second;
        for (auto& node_0 : graph.nodes[0]) {
            for (auto& node_1_log_w : graph.right_edges[node_0]) {
                auto node_1 = std::get<0>(node_1_log_w);
                auto log_weight = std::get<1>(node_1_log_w);
                auto nuc = std::get<2>(node_0);
                //std::cout <<  nuc << std::endl;

                if (!best_path_in_one_codon_unit[aa].count(make_tuple(node_0,node_1)))
                    best_path_in_one_codon_unit[aa][make_tuple(node_0,node_1)] = make_tuple(util::value_min<double>(),k_void_nuc,k_void_nuc,k_void_nuc);

                double current_log_weight = std::get<0>(best_path_in_one_codon_unit[aa][make_tuple(node_0,node_1)]);
                if (current_log_weight < log_weight) {
                    best_path_in_one_codon_unit[aa][make_tuple(node_0,node_1)] = make_tuple(log_weight,nuc,k_void_nuc,k_void_nuc);
                }
            }
        }

        for (auto& node_1 : graph.nodes[1]) {
            for (auto& node_2_log_w : graph.right_edges[node_1]) {
                auto node_2 = std::get<0>(node_2_log_w);
                auto log_weight = std::get<1>(node_2_log_w);
                auto nuc = std::get<2>(node_1);
                //std::cout << nuc << std::endl;

                if (!best_path_in_one_codon_unit[aa].count(make_tuple(node_1,node_2)))
                    best_path_in_one_codon_unit[aa][make_tuple(node_1,node_2)] = make_tuple(util::value_min<double>(),k_void_nuc,k_void_nuc,k_void_nuc);

                double current_log_weight = std::get<0>(best_path_in_one_codon_unit[aa][make_tuple(node_1,node_2)]);
                if (current_log_weight < log_weight) {
                    best_path_in_one_codon_unit[aa][make_tuple(node_1,node_2)] = make_tuple(log_weight,nuc,k_void_nuc,k_void_nuc);
                }

                auto temp = best_path_in_one_codon_unit[aa][make_tuple(node_1,node_2)];
            }
        }

        for (auto& node_2 : graph.nodes[2]) {
            for (auto& node_3_log_w : graph.right_edges[node_2]) {
                auto node_3 = std::get<0>(node_3_log_w);
                auto log_weight = std::get<1>(node_3_log_w);
                auto nuc = std::get<2>(node_2);

                if (!best_path_in_one_codon_unit[aa].count(make_tuple(node_2,node_3)))
                    best_path_in_one_codon_unit[aa][make_tuple(node_2,node_3)] = make_tuple(util::value_min<double>(),k_void_nuc,k_void_nuc,k_void_nuc);

                double current_log_weight = std::get<0>(best_path_in_one_codon_unit[aa][make_tuple(node_2,node_3)]);
                if (current_log_weight < log_weight) {
                    best_path_in_one_codon_unit[aa][make_tuple(node_2,node_3)] = make_tuple(log_weight,nuc, k_void_nuc ,k_void_nuc);
                }
            }
        }

        for (auto& node_3 : graph.nodes[3]) {
            for (auto& node_4_log_w : graph.right_edges[node_3]) {
                auto node_4 = std::get<0>(node_4_log_w);
                auto log_weight = std::get<1>(node_4_log_w);
                auto nuc = std::get<2>(node_3);

                if (!best_path_in_one_codon_unit[aa].count(make_tuple(node_3,node_4)))
                    best_path_in_one_codon_unit[aa][make_tuple(node_3,node_4)] = make_tuple(util::value_min<double>(),k_void_nuc,k_void_nuc,k_void_nuc);

                double current_log_weight = std::get<0>(best_path_in_one_codon_unit[aa][make_tuple(node_3,node_4)]);
                if (current_log_weight < log_weight) {
                    best_path_in_one_codon_unit[aa][make_tuple(node_3,node_4)] = make_tuple(log_weight,nuc,k_void_nuc,k_void_nuc);
                }
            }
        }


        for (auto& node_0 : graph.nodes[0]) {
            for (auto& node_1_log_weight_0 : graph.right_edges[node_0]) {
                auto& node_1 = std::get<0>(node_1_log_weight_0);
                auto log_weight_0 = std::get<1>(node_1_log_weight_0);
                auto& nuc_0 = std::get<2>(node_0);
                for (auto& node_2_log_weight_1 : graph.right_edges[node_1]) {
                    auto& node_2 = std::get<0>(node_2_log_weight_1);
                    auto log_weight_1 = std::get<1>(node_2_log_weight_1);
                    auto& nuc_1 = std::get<2>(node_1);

                    if (!best_path_in_one_codon_unit[aa].count(make_tuple(node_0,node_2)))
                        best_path_in_one_codon_unit[aa][make_tuple(node_0,node_2)] = make_tuple(util::value_min<double>(),k_void_nuc,k_void_nuc,k_void_nuc);

                    if (std::get<0>(best_path_in_one_codon_unit[aa][make_tuple(node_0,node_2)]) < log_weight_0 + log_weight_1)
                        best_path_in_one_codon_unit[aa][make_tuple(node_0,node_2)] = make_tuple(log_weight_0 + log_weight_1, nuc_0, nuc_1,k_void_nuc);
                }
            }
        }

        for (auto& node_1 : graph.nodes[1]) {
            for (auto& node_2_log_weight_1 : graph.right_edges[node_1]) {
                auto& node_2 = std::get<0>(node_2_log_weight_1);
                auto log_weight_1 = std::get<1>(node_2_log_weight_1);
                auto& nuc_1 = std::get<2>(node_1);
                for (auto& node_3_log_weight_2 : graph.right_edges[node_2]) {
                    auto& node_3 = std::get<0>(node_3_log_weight_2);
                    auto log_weight_2 = std::get<1>(node_3_log_weight_2);
                    auto& nuc_2 = std::get<2>(node_2);

                    if (!best_path_in_one_codon_unit[aa].count(make_tuple(node_1,node_3)))
                        best_path_in_one_codon_unit[aa][make_tuple(node_1,node_3)] = make_tuple(util::value_min<double>(),k_void_nuc,k_void_nuc,k_void_nuc);

                    if (std::get<0>(best_path_in_one_codon_unit[aa][make_tuple(node_1,node_3)]) < log_weight_1 + log_weight_2)
                        best_path_in_one_codon_unit[aa][make_tuple(node_1,node_3)] = make_tuple(log_weight_1 + log_weight_2, nuc_1, nuc_2,k_void_nuc);
                }
            }
        }

        for (auto& node_2 : graph.nodes[2]) {
            for (auto& node_3_log_weight_2 : graph.right_edges[node_2]) {
                auto& node_3 = std::get<0>(node_3_log_weight_2);
                auto log_weight_2 = std::get<1>(node_3_log_weight_2);
                auto& nuc_2 = std::get<2>(node_2);
                for (auto& node_4_log_weight_3 : graph.right_edges[node_3]) {
                    auto& node_4 = std::get<0>(node_4_log_weight_3);
                    auto log_weight_3 = std::get<1>(node_4_log_weight_3);
                    auto& nuc_3 = std::get<2>(node_3);

                    if (!best_path_in_one_codon_unit[aa].count(make_tuple(node_2,node_4)))
                        best_path_in_one_codon_unit[aa][make_tuple(node_2,node_4)] = make_tuple(util::value_min<double>(),k_void_nuc,k_void_nuc,k_void_nuc);

                    if (std::get<0>(best_path_in_one_codon_unit[aa][make_tuple(node_2,node_4)]) < log_weight_2 + log_weight_3)
                        best_path_in_one_codon_unit[aa][make_tuple(node_2,node_4)] = make_tuple(log_weight_2 + log_weight_3, nuc_2, nuc_3,k_void_nuc);
                }
            }
        }

        for (auto& node_0 : graph.nodes[0]) {
            for (auto& node_1_log_weight_0 : graph.right_edges[node_0]) {
                auto& node_1 = std::get<0>(node_1_log_weight_0);
                auto log_weight_0 = std::get<1>(node_1_log_weight_0);
                auto& nuc_0 = std::get<2>(node_0);
                for (auto& node_2_log_weight_1 : graph.right_edges[node_1]) {
                    auto& node_2 = std::get<0>(node_2_log_weight_1);
                    auto log_weight_1 = std::get<1>(node_2_log_weight_1);
                    auto& nuc_1 = std::get<2>(node_1);
                    for (auto& node_3_log_weight_2 : graph.right_edges[node_2]) {
                        auto& node_3 = std::get<0>(node_3_log_weight_2);
                        auto log_weight_2 = std::get<1>(node_3_log_weight_2);
                        auto& nuc_2 = std::get<2>(node_2);

                        if (!best_path_in_one_codon_unit[aa].count(make_tuple(node_0,node_3)))
                            best_path_in_one_codon_unit[aa][make_tuple(node_0,node_3)] = make_tuple(util::value_min<double>(),k_void_nuc,k_void_nuc,k_void_nuc);

                        if (std::get<0>(best_path_in_one_codon_unit[aa][make_tuple(node_0,node_3)]) < log_weight_0 + log_weight_1 + log_weight_2)
                            best_path_in_one_codon_unit[aa][make_tuple(node_0,node_3)] = make_tuple(log_weight_0 + log_weight_1+log_weight_2, nuc_0, nuc_1, nuc_2);
                    }
                }
            }
        }
        
        for (auto& node_1 : graph.nodes[1]) {
            for (auto& node_2_log_weight_1 : graph.right_edges[node_1]) {
                auto& node_2 = std::get<0>(node_2_log_weight_1);
                auto log_weight_1 = std::get<1>(node_2_log_weight_1);
                auto& nuc_1 = std::get<2>(node_1);
                for (auto& node_3_log_weight_2 : graph.right_edges[node_2]) {
                    auto& node_3 = std::get<0>(node_3_log_weight_2);
                    auto log_weight_2 = std::get<1>(node_3_log_weight_2);
                    auto& nuc_2 = std::get<2>(node_2);
                    for (auto& node_4_log_weight_3 : graph.right_edges[node_3]) {
                        auto& node_4 = std::get<0>(node_4_log_weight_3);
                        auto log_weight_3 = std::get<1>(node_4_log_weight_3);
                        auto& nuc_3 = std::get<2>(node_3);

                        if (!best_path_in_one_codon_unit[aa].count(make_tuple(node_1,node_4)))
                            best_path_in_one_codon_unit[aa][make_tuple(node_1,node_4)] = make_tuple(util::value_min<double>(),k_void_nuc,k_void_nuc,k_void_nuc);

                        if (std::get<0>(best_path_in_one_codon_unit[aa][make_tuple(node_1,node_4)]) < log_weight_1 + log_weight_2 + log_weight_3)
                            best_path_in_one_codon_unit[aa][make_tuple(node_1,node_4)] = make_tuple(log_weight_1+log_weight_2+log_weight_3, nuc_1, nuc_2, nuc_3);
                    }
                }
            }
        }
    }

    // for(int k=0;k<4;k++){
    //     for(auto &node : aa_graphs_with_ln_weights["Leu"].nodes[k]){
    //         IndexType index0;
    //         NumType num0;
    //         NucType nuc0;
    //         std::tie(index0,num0,nuc0) = node;
    //         for (auto &node_weight : aa_graphs_with_ln_weights["Leu"].right_edges[node]){
    //             IndexType index_nr;
    //             NumType num_nr;
    //             NucType nuc_nr;
    //             NodeType node_right = std::get<0>(node_weight);
    //             std::tie(index_nr,num_nr,nuc_nr) = node_right;
    //             std::cout << index0 << "," << num0 << "," << GET_ACGU(nuc0) << std::endl;
    //             std::cout << index_nr << "," << num_nr << "," << GET_ACGU(nuc_nr)  << std::endl;
    //             double weight;
    //             NucType nuc1;
    //             NucType nuc2;
    //             NucType nuc3;
    //             std::tie(weight,nuc1,nuc2,nuc3) = best_path_in_one_codon_unit["Leu"][make_tuple(node,node_right)];
    //             std::cout << weight << ","<<  GET_ACGU(nuc1) << ","<<  GET_ACGU(nuc2) << ","<<  GET_ACGU(nuc3) << std::endl;
    //             std::cout <<  ".............."  << std::endl;
    //         }            
    //     }
    // }

    // std::cout <<  "2................................."  << std::endl;

    // for(int k=0;k<3;k++){
    //     for(auto &node : aa_graphs_with_ln_weights["Leu"].nodes[k]){
    //         IndexType index0;
    //         NumType num0;
    //         NucType nuc0;
    //         std::tie(index0,num0,nuc0) = node;
    //         std::cout << index0 << "," << num0 << "," << GET_ACGU(nuc0) << std::endl;
    //         for (auto &node_weight : aa_graphs_with_ln_weights["Leu"].right_edges[node]){
    //             // IndexType index_nr;
    //             // NumType num_nr;
    //             // NucType nuc_nr;
    //             NodeType node_right = std::get<0>(node_weight);
    //             // std::tie(index_nr,num_nr,nuc_nr) = node_right;
    //             // std::cout << index_nr << "," << num_nr << "," << GET_ACGU(nuc_nr)  << std::endl;
    //             for (auto &node_weight2 : aa_graphs_with_ln_weights["Leu"].right_edges[node_right]){
    //                 IndexType index_2;
    //                 NumType num_2;
    //                 NucType nuc_2;
    //                 NodeType node2 = std::get<0>(node_weight2);
    //                 std::tie(index_2,num_2,nuc_2) = node2;
    //                 std::cout << index_2 << "," << num_2 << "," << GET_ACGU(nuc_2)  << std::endl;              
    //                 double weight;
    //                 NucType nuc1;
    //                 NucType nuc2;
    //                 NucType nuc3;
    //                 std::tie(weight,nuc1,nuc2,nuc3) = best_path_in_one_codon_unit["Leu"][make_tuple(node,node2)];
    //                 std::cout << weight << ","<<  GET_ACGU(nuc1) << ","<<  GET_ACGU(nuc2) << ","<<  GET_ACGU(nuc3) << std::endl;
    //                 std::cout <<  ".............."  << std::endl;
    //             }
    //         }
    //     }
    //     std::cout <<  ".............."  << std::endl;
    // }

    // std::cout <<  "3................................."  << std::endl;

    // for(int k=0;k<2;k++){
    //     for(auto &node : aa_graphs_with_ln_weights["Leu"].nodes[k]){
    //         IndexType index0;
    //         NumType num0;
    //         NucType nuc0;
    //         std::tie(index0,num0,nuc0) = node;
    //         std::cout << index0 << "," << num0 << "," << GET_ACGU(nuc0) << std::endl;
    //         for (auto &node_weight1 : aa_graphs_with_ln_weights["Leu"].right_edges[node]){
    //             NodeType node1 = std::get<0>(node_weight1);
    //             for (auto &node_weight2 : aa_graphs_with_ln_weights["Leu"].right_edges[node1]){
    //                 NodeType node2 = std::get<0>(node_weight2);
    //                 for (auto &node_weight3 : aa_graphs_with_ln_weights["Leu"].right_edges[node2]){
    //                     IndexType index_3;
    //                     NumType num_3;
    //                     NucType nuc_3;
    //                     NodeType node3 = std::get<0>(node_weight3);
    //                     std::tie(index_3,num_3,nuc_3) = node3;
    //                     std::cout << index_3 << "," << num_3 << "," << GET_ACGU(nuc_3)  << std::endl;                
    //                     double weight;
    //                     NucType nuc1;
    //                     NucType nuc2;
    //                     NucType nuc3;
    //                     std::tie(weight,nuc1,nuc2,nuc3) = best_path_in_one_codon_unit["Leu"][make_tuple(node,node3)];
    //                     std::cout << weight << ","<<  GET_ACGU(nuc1) << ","<<  GET_ACGU(nuc2) << ","<<  GET_ACGU(nuc3) << std::endl;
    //                     std::cout <<  ".............."  << std::endl;
    //                 }
    //             }
    //         }
    //     }
    //     std::cout <<  ".............."  << std::endl;
    // }



    std::unordered_map<std::string, double> max_path;
    std::unordered_map<std::string, std::string> aa_best_path_in_a_whole_codon;

    for (auto& aa_path_weight : codon.aa_table_) {
        auto& aa = aa_path_weight.first; // char
        for (auto& path_weight : aa_path_weight.second) {
            if (max_path[aa] < path_weight.second) {
                max_path[aa] = path_weight.second;
                aa_best_path_in_a_whole_codon[aa] = path_weight.first;
            }
        }
    }

    aa_graphs_with_ln_weights_ret = aa_graphs_with_ln_weights;
    best_path_in_one_codon_unit_ret = best_path_in_one_codon_unit;
    aa_best_path_in_a_whole_codon_ret = aa_best_path_in_a_whole_codon;
}



template <typename IndexType,
          typename NodeType = tuple<IndexType, NumType, NucType>,
          typename LatticeType = Lattice<IndexType>,
          typename DFAType = DFA<IndexType>>
DFAType get_dfa(unordered_map<string, LatticeType> aa_graphs, vector<string> aa_seq) {
    DFAType dfa = DFAType();
    DFAType dfa1 = DFAType();                 
    NodeType newnode = make_tuple(4 * static_cast<IndexType>(aa_seq.size()), 0,0);
    dfa.add_node(newnode);
    IndexType i = 0;
    IndexType i4;
    string aa;
    LatticeType graph;
    for(auto& item : aa_seq) {
        i4 = i * 4;
        aa = aa_seq[i];
        graph = aa_graphs[aa];
        for (IndexType pos = 0; pos <= 3; pos++) {
            for(auto& node : graph.nodes[pos]) {
                IndexType num = get<1>(node);
                NucType nuc = get<2>(node);
                newnode = make_tuple(i4 + pos,num,nuc);
                dfa.add_node(newnode);
                for (auto& edge : graph.right_edges[node]) {
                    NodeType n2 = get<0>(edge);
                    num = get<1>(n2);
                    NodeType newn2 = make_tuple(i4 + pos + 1, num, get<2>(n2));
                    dfa.add_edge(newnode, newn2, get<1>(edge));
                }
            }
        }
        i++;
    }

    IndexType j=0;
    NodeType endn = make_tuple(3 * static_cast<IndexType>(aa_seq.size()), 0 , 0);
    dfa1.add_node(endn);
    for(i=0; i< aa_seq.size()*4;i=i+4) {
        for (IndexType pos = 1; pos <= 3; pos++) {
            for(auto& node : dfa.nodes[pos+i]) {
                IndexType num = get<1>(node);
                NucType nuc = get<2>(node);
                //std::cout <<  nuc << std::endl;
                newnode = make_tuple(j+pos-1,num,nuc);
                dfa1.add_node(newnode);

                if(pos!=3){
                    for (auto& edge : dfa.right_edges[node]) {
                        NodeType n2 = get<0>(edge);
                        num = get<1>(n2);
                        NodeType newn2 = make_tuple(j + pos, num, get<2>(n2));
                        dfa1.add_edge(newnode, newn2, get<1>(edge));
                    }
                }else{
                    for(auto& edge : dfa.right_edges[node]){
                        if(i+5<4 * static_cast<IndexType>(aa_seq.size())){
                            for(auto& n2 : dfa.nodes[i+5]){
                                num = get<1>(n2);
                                NodeType newn2 = make_tuple(j + pos, num, get<2>(n2));
                                dfa1.add_edge(newnode, newn2, get<1>(edge));
                            }
                        }else{
                            dfa1.add_edge(newnode, endn, get<1>(edge));
                        }
                        
                    }
                }               
            }
        }
        j=j+3;
    }


#ifdef is_verbose
    printf("-----------------DFA------------------------\n");
    for(IndexType pos = 0; pos < 3 * static_cast<IndexType>(aa_seq.size()) + 1; pos++){
        for(auto& node : dfa.nodes[pos]) {
            IndexType p = get<0>(node);
            IndexType num = get<1>(node);
            printf("node, (%d, %d)\n", p, num);
            for(auto &n2 : dfa.auxiliary_right_edges[node]){
                IndexType p2 = get<0>(n2.first);
                IndexType num2 = get<1>(n2.first);
                for(auto nuc : n2.second){
                    printf("              (%d, %d) -(%d,%lf)-> (%d, %d)\n", p, num, get<0>(nuc),get<1>(nuc), p2, num2);
                }
            }
            for(auto &n1 : dfa.auxiliary_left_edges[node]){
                IndexType p1 = get<0>(n1.first); IndexType num1 = get<1>(n1.first);
                for(auto nuc : n1.second){
                    printf("  (%d, %d) <-(%d,%lf)- (%d, %d)\n", p1, num1, get<0>(nuc),get<1>(nuc), p, num);
                }
            }
        }
    }
#endif

    return dfa1;
}

}

