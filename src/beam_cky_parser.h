#pragma once

#include <memory>
#include <string>
#include <cstring>
#include <limits>
#include <vector>
#include <unordered_map>
#include <thread>
#include <mutex>

#include "Utils/network.h"
#include "Utils/codon.h"
#include "Utils/flat.h"


namespace LinearDesign {

namespace detail {

    struct NodeNucIndex {
        LINEAR_DESIGN_INLINE size_t operator()(const NodeNucpair& node_nucpair) const {
            return (node_nucpair.node_first << 6) | (node_nucpair.node_second << 3) | node_nucpair.nucpair;
        }
    };

    struct NodeNucReverseIndex {
        LINEAR_DESIGN_INLINE NodeNucpair operator()(const size_t index) const {
            NodeNucpair node_nucpair = {IndexType(index >> 6), NumType((index & 0x3f) >> 3), NucPairType(index & 0x7)};
            return node_nucpair;
        }
    };

    struct NodeIndex {
        LINEAR_DESIGN_INLINE size_t operator()(const NodeType& node) const {
            return (std::get<0>(node) << 6) | (std::get<1>(node) << 3 )| std::get<2>(node);
        }
    };

    // struct NodeNucReverseIndex2 {
    //     LINEAR_DESIGN_INLINE NodeType operator()(const size_t index) const {
    //         return {index >> 1, (index & 0x1)};
    //     }
    // };

} /* detail */    

template <typename IndexType,
          typename NucType = IndexType,
          typename NodeType = std::tuple<IndexType, NumType, NucType>,
          typename DFAType = DFA<IndexType>>
string get_nuc_from_dfa_cai(DFAType& dfa, const NodeType& start_node, const NodeType& end_node,
        const std::vector<std::string>& protein, std::unordered_map<std::string, std::unordered_map<std::tuple<NodeType, NodeType>, 
        std::tuple<double, NucType, NucType, NucType>, std::hash<std::tuple<NodeType, NodeType>>>>&
        best_path_in_one_codon_unit, std::unordered_map<std::string, std::string>& aa_best_path_in_a_whole_codon) {

    IndexType s_index = std::get<0>(start_node);
    IndexType t_index = std::get<0>(end_node);

    if (s_index >= t_index)
        return "";

    auto aa_left = protein[s_index / 3]; // tri letter
    auto aa_right = protein[t_index / 3];
    auto start_node_re_index = make_tuple(s_index % 3+1, std::get<1>(start_node), std::get<2>(start_node));
    auto end_node_re_index = make_tuple(t_index % 3+1, std::get<1>(end_node), std::get<2>(end_node));
    if (t_index - s_index < 3) {
        if (s_index / 3 == t_index / 3) {
            std::string temp_seq = "";
            auto& nucs = best_path_in_one_codon_unit[aa_left][make_tuple(start_node_re_index, end_node_re_index)];
            temp_seq.append(1, GET_ACGU(std::get<1>(nucs)));
            if (std::get<2>(nucs) != k_void_nuc) 
                temp_seq.append(1, GET_ACGU(std::get<2>(nucs)));

            if (temp_seq.length() != std::get<0>(end_node) - std::get<0>(start_node)) {
                assert(false);
            }
            return temp_seq;
        } else {
            std::string temp_left = "";
            std::string temp_right = "";
            if (s_index % 3 != 0) {
                auto& nucs = best_path_in_one_codon_unit[aa_left][make_tuple(start_node_re_index, make_tuple(0, 0, 0))];
                temp_left.append(1, GET_ACGU(std::get<1>(nucs)));
                if (std::get<2>(nucs) != k_void_nuc) 
                    temp_left.append(1, GET_ACGU(std::get<2>(nucs)));
            }

            if (t_index % 3 != 0) {
                auto& nucs = best_path_in_one_codon_unit[aa_right][make_tuple(make_tuple(0, 0, 0), end_node_re_index)];                
                if (GET_ACGU(std::get<1>(nucs)) != 'X') 
                    temp_right.append(1, GET_ACGU(std::get<1>(nucs)));
                if (std::get<2>(nucs) != k_void_nuc)
                    temp_right.append(1, GET_ACGU(std::get<2>(nucs)));
            }

            assert((temp_left + temp_right).length() == std::get<0>(end_node) - std::get<0>(start_node));

            return temp_left + temp_right;
        }

    } else {
        std::string temp_left = "";
        std::string temp_mid = "";
        std::string temp_right = "";

        if (s_index % 3 != 0) {
            auto& nucs = best_path_in_one_codon_unit[aa_left][make_tuple(start_node_re_index, make_tuple(0, 0, 0))];
            temp_left.append(1, GET_ACGU(std::get<1>(nucs)));
            if (std::get<2>(nucs) != k_void_nuc) 
                temp_left.append(1, GET_ACGU(std::get<2>(nucs)));
        }

        IndexType protein_start_index = s_index / 3;
        if (s_index % 3 != 0)
            protein_start_index++;

        IndexType protein_end_index = t_index / 3;

        if (protein_start_index != protein_end_index) {
            for (IndexType protein_index = protein_start_index; protein_index < protein_end_index; ++protein_index) {
                
                std::string nucs;
                auto aa_tri = protein[protein_index];
                if (k_map_3_1.count(aa_tri)) {
                    nucs = aa_best_path_in_a_whole_codon[std::string(1, k_map_3_1[aa_tri])];
                } else if (aa_best_path_in_a_whole_codon.count(aa_tri)) {
                    nucs = aa_best_path_in_a_whole_codon[aa_tri];
                } else {
                    assert(false);
                }

                for (auto nuc : nucs) {
                    temp_mid.append(1, nuc);
                }
            }
        }

        if (t_index % 3 != 0) {
            auto& nucs = best_path_in_one_codon_unit[aa_right][make_tuple(make_tuple(0,0,0), end_node_re_index)];
            temp_right.append(1, GET_ACGU(std::get<2>(nucs)));
            if (std::get<3>(nucs) != k_void_nuc)
                temp_right.append(1, GET_ACGU(std::get<3>(nucs)));
        }

        assert((temp_left + temp_mid + temp_right).length() == std::get<0>(end_node) - std::get<0>(start_node));

        return temp_left + temp_mid + temp_right;
    }
}

template <typename ScoreType,
          typename IndexType,
          typename NodeType = tuple<IndexType, NumType, NucType>>
class BeamCKYParser {
public:
    using State_t = State<ScoreType>;
    using DFA_t = DFA<IndexType>;
    using ScoreInnerDate_t = ScoreInnerDate<ScoreType, IndexType>;
    using NextPair_t = vector<unordered_map<NodeType, vector<tuple<NodeType, NucType, double>>, std::hash<NodeType>>>;
    using NextPairSet_t = vector<unordered_map<NodeType, set<tuple<NodeType, NucType, double>>, std::hash<NodeType>>>;
    using PrefixScore_t = unordered_map<NodeType, ScoreType, std::hash<NodeType>>;
    using BestX_t_CAI = Flat<NodeType, Flat<NodeNucpair, State_t, detail::NodeNucIndex>, detail::NodeIndex>;
    using BestM_t_CAI = Flat<NodeType, Flat<NodeType, State_t, detail::NodeIndex>, detail::NodeIndex>;
    using BestC_t_CAI = Flat<NodeType, State_t, detail::NodeIndex>;
    using Broken_codon_t_CAI = Flat<NodeType, Flat<NodeType, double, detail::NodeIndex>, detail::NodeIndex>;

    BeamCKYParser(const double lambda_value, const bool verbose);

    DecoderResult<double, IndexType> parse(DFA_t& dfa, 
        Codon& codon, 
        std::string& aa_seq, 
        std::vector<std::string>& p, 
        std::unordered_map<std::string, std::string>& aa_best_in_codon,
        std::unordered_map<std::string, std::unordered_map<std::tuple<NodeType, NodeType>, 
        std::tuple<double, NucType, NucType, NucType>, std::hash<std::tuple<NodeType, NodeType>>>>& best_path_in_one_codon,
        std::unordered_map<string, Lattice<IndexType>>& aa_graphs_with_ln_weights);

private:
    
    template <IndexType j_num>
    void hairpin_beam(IndexType j, DFA_t& dfa);
    
    template <IndexType j_num>
    void Multi_beam(IndexType j, DFA_t& dfa);
    
    template <IndexType j_num>
    void P_beam(IndexType j, DFA_t& dfa);
    
    template <IndexType j_num>
    void M2_beam(IndexType j, DFA_t& dfa);
    
    template <IndexType j_num>
    void M_beam(IndexType j, DFA_t& dfa);
    
    template <IndexType j_num>
    void C_beam(IndexType j, DFA_t& dfa);

    void update_if_better(State_t &state, const ScoreType newscore, const double cai_score) {
        if (state.score + state.cai_score < newscore + cai_score) {
            state.score = newscore;
            state.cai_score = cai_score;
            double epsilon = 1e-8;
        }
    }

    void update_if_better(State_t &state, const ScoreType newscore, const double cai_score, const NodeType pre_node, const double pre_left_cai) {
        if (state.score + state.cai_score < newscore + cai_score) {
            state.score = newscore;
            state.cai_score = cai_score;
            state.pre_node = pre_node;
            state.pre_left_cai = pre_left_cai;
        }
    }


    void get_next_pair(DFA_t& dfa);
    void get_next_pair_set();

    void get_prev_pair(DFA_t& dfa);
    void get_prev_pair_set();

    void preprocess(DFA_t& dfa);    

    BacktraceResult backtrace(DFA_t& dfa, const State_t& state, NodeType end_node);

    ScoreType quickselect_partition(std::vector<ScoreInnerDate_t>& scores,
        ScoreType lower, ScoreType upper);

    ScoreType quickselect(std::vector<ScoreInnerDate_t>& scores,
        const ScoreType lower, const ScoreType upper, const IndexType k);

    double get_broken_codon_score(const NodeType& start_node, const NodeType& end_node);

    double lambda;
    bool is_verbose;

    IndexType seq_length; 

    BestX_t_CAI bestH, bestP, bestMulti;
    BestM_t_CAI bestM2, bestM, bestM_P; // hzhang: bestM_P
    BestC_t_CAI bestC;

    detail::NodeNucReverseIndex reverse_index;
    //detail::NodeNucReverseIndex2 reverse_index2;

    NextPair_t next_pair;
    NextPairSet_t next_pair_set;

    NextPair_t prev_pair;
    NextPairSet_t prev_pair_set;

    NextPair_t next_list;
    NextPair_t prev_list;

    vector<vector<vector<ScoreType>>> bulge_score;
    vector<vector<ScoreType>> stacking_score;

    std::unordered_map<string, Lattice<IndexType>> aa_graphs_with_ln_w;

    std::vector<std::string> protein;
    std::unordered_map<std::string, std::string> aa_best_path_in_a_whole_codon;
    std::unordered_map<std::string, 
                       std::unordered_map<std::tuple<NodeType, NodeType>, std::tuple<FinalScoreType, NucType, NucType, NucType>, 
                       std::hash<std::tuple<NodeType, NodeType>>>> best_path_in_one_codon_unit;

    Broken_codon_t_CAI get_broken_codon_score_map;

#ifdef SPECIAL_HP
    unordered_map<NodeType, unordered_map<NodeType, unordered_map<int8_t, vector<tuple<string, ScoreType, FinalScoreType>>>, std::hash<NodeType>>, std::hash<NodeType>> hairpin_seq_score_cai;
    void special_hp(DFA_t& dfa, int8_t hairpin_length);
#endif
};

}
