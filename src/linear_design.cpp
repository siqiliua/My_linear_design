#include <iomanip>
#include <string>
#include <vector>
#include <iostream>
#include "beam_cky_parser.h"
#include "beam_cky_parser.cc"
#include "Utils/codon.h"
#include "Utils/common.h"

using namespace LinearDesign;
//using namespace std;

#ifndef CODING_WHEEL
#define CODING_WHEEL "./coding_wheel.txt"
#endif


template <typename ScoreType, typename IndexType>
bool output_result(const DecoderResult<ScoreType, IndexType>& result, 
        const double duration, const double lambda, const bool is_verbose, 
        const Codon& codon, std::string& CODON_TABLE) {

    std::stringstream ss;
    if (is_verbose)
        ss << "Using lambda = " << (lambda / 100.) << "; Using codon frequency table = " << CODON_TABLE << endl;
    ss << "mRNA sequence:  " << result.sequence << endl;
    ss << "mRNA structure: " << result.structure << endl;
    ss << "mRNA folding free energy: " << std::setprecision(2) << fixed << result.score 
                                        << " kcal/mol; mRNA CAI: " << std::setprecision(3) 
                                        << fixed << codon.calc_cai(result.sequence) << endl;
    if (is_verbose)
        ss << "Runtime: " << duration << " seconds" << endl;
    cout << ss.str() << endl;

    return true;
}

void linearDesign(int n, vector<std::string> aa_seq_list){
    // default args
    double lambda = 0.0f;
    bool is_verbose = false;
    std::string CODON_TABLE = "./codon_usage_freq_table_human.csv";

    // parse args
    // if (argc != 4) {
    //     show_usage();
    //     return 1;
    // }else{
    //     lambda = atof(argv[1]);
    //     is_verbose = atoi(argv[2]) == 1;
    //     if (string(argv[3]) != ""){
    //         CODON_TABLE = argv[3];
    //     }
    // } 
    // lambda *= 100.;

    Codon codon(CODON_TABLE);
    //std::cout << codon.aa_table_["A"][0].first << std::endl;

    std::unordered_map<std::string, Lattice<IndexType>> aa_graphs_with_ln_weights;
    std::unordered_map<std::string, std::unordered_map<std::tuple<NodeType, NodeType>, std::tuple<double, NucType, NucType>, std::hash<std::tuple<NodeType, NodeType>>>> best_path_in_one_codon_unit;
    std::unordered_map<std::string, std::string> aa_best_path_in_a_whole_codon;
    prepare_codon_unit_lattice<IndexType>(CODING_WHEEL, codon, aa_graphs_with_ln_weights, best_path_in_one_codon_unit, aa_best_path_in_a_whole_codon, lambda);

    //std::cout << aa_best_path_in_a_whole_codon["A"] << std::endl;
    // for(int i=0;i<4;i++){
    //     int size;
    //     if(i==0){
    //         size=1;
    //     }else if(i==1){
    //         size=2;
    //     }else if(i==2){
    //         size=2;
    //     }else if(i==3){
    //         size=6;
    //     }
         
    //     for(int j=0;j<size;j++){
    //         NodeType node;
    //         IndexType index;
    //         NumType num;
    //         NucType nuc;
    //         node = aa_graphs_with_ln_weights["Leu"].nodes[i][j];
    //         std::tie(index,num,nuc) = node;
    //         std::cout << index << "," << num << "," << nuc << std::endl;

    //         for (auto &node_weight : aa_graphs_with_ln_weights["Leu"].left_edges[node]){
    //             IndexType index_nr;
    //             NumType num_nr;
    //             NucType nuc_nr;
    //             NodeType node_right = std::get<0>(node_weight);
    //             double weight = std::get<1>(node_weight);
    //             std::tie(index_nr,num_nr,nuc_nr) = node_right;
    //             std::cout << index_nr << "," << num_nr << "," << nuc_nr << "," << weight << std::endl;
    //             //
    //             double weight1;
    //             NucType nuc1;
    //             NucType nuc2;
    //             std::tie(weight1,nuc1,nuc2) = best_path_in_one_codon_unit["Leu"][make_tuple(node,node_right)];
    //             std::cout << weight1 << ","<<  nuc1 << ","<<  nuc2 << std::endl;
    //         }
    //          std::cout <<  "........."  << std::endl;

    //     }

    // }
    
    // main loop
    string aa_tri_seq;
    // load input
    for(int i = 0; i < n; i++){
        string& aa_seq = aa_seq_list[i];
        // convert to uppercase
        transform(aa_seq.begin(), aa_seq.end(), aa_seq.begin(), ::toupper);
        aa_tri_seq.clear();
        if (is_verbose)
            cout << "Input protein: " << aa_seq << endl;
        if (!cvt_to_seq(aa_seq, aa_tri_seq)) 
            continue;
        std::cout << aa_tri_seq << std::endl;

        auto protein = split(aa_tri_seq, ' ');
        //std::cout << protein[0] << std::endl;

        //init parser
        BeamCKYParser<ScoreType, IndexType> parser(lambda, is_verbose);

        // parse
        auto system_start = chrono::system_clock::now();
        auto dfa = get_dfa<IndexType>(aa_graphs_with_ln_weights, split(aa_tri_seq, ' '));
        //std::cout << std::get<2>(dfa.nodes[1][0]) << std::endl;

        // for(int i=0;i<31;i++){
        //     for(auto &node : dfa.nodes[i]){
        //         IndexType index;
        //         NumType num;
        //         NucType nuc;
        //         std::tie(index,num,nuc) = node;
        //         std::cout << index << "," << num << "," << nuc << std::endl;
        //         for (auto &node_weight : dfa.auxiliary_right_edges[node]){
        //             IndexType index_nr;
        //             NumType num_nr;
        //             NucType nuc_nr;
        //             NodeType node_left = node_weight.first;
        //             double weight = node_weight.second;
        //             std::tie(index_nr,num_nr,nuc_nr) = node_left;
        //             std::cout << index_nr << "," << num_nr << "," << nuc_nr << "," << weight << std::endl;
        //         }
        //     std::cout <<  "........."  << std::endl;
        //     }
        // }
    
        auto result = parser.parse(dfa, codon, aa_seq, protein, aa_best_path_in_a_whole_codon, best_path_in_one_codon_unit, aa_graphs_with_ln_weights);
        auto system_diff = chrono::system_clock::now() - system_start;
        auto system_duration = chrono::duration<double>(system_diff).count();  

        // output
        output_result(result, system_duration, lambda, is_verbose, codon, CODON_TABLE);
    }
    
}


int main()
{
    int n = 1;
    std::vector<std::string> A; 
    A.push_back("ML*");
    //A.push_back("MPNTL*");
    linearDesign(n, A);
    return 0;
    
}



