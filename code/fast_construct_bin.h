#pragma once

#include <vector>
#include <lemon/list_graph.h>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <queue>
#include <cmath>
#include <utility>
#include <tuple>
#include <fstream>
#include <string>

// @brief Given a number of bins and the higher bound of elements per bin, gather sequences, such that similar sequences are within the same bin and dissimilar sequences are in different bins.
// @param labMaps : The clustering for the graph used.
// @param level_clusters : For every level and sequence, a pointer to the corresponding cluster in the higher level is returned.
// @param bins : the number of bins used
// @param t_max : The maximum number of elements per bin.
template <typename Hasher>
std::vector<std::vector<size_t>> binning_base(const std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, Hasher>>& labMaps,
                                        const std::vector<std::unordered_map<size_t,const std::vector<size_t>*>>& level_clusters,
                                        const size_t bins, const size_t t_max);
#include "templates/fast_construct_binning_base.tpp"

// @brief Given a number of bins and the higher bound of elements per bin, gather sequences, such that similar sequences are within the same bin and dissimilar sequences are in different bins.
// @param labMaps : The clustering for the graph used.
// @param level_clusters : For every level and sequence, a pointer to the corresponding cluster in the higher level is returned.
// @param fracmin_sketches : The Set of Fracmin sketches corresponding to every sequence per ID
// @param s : The Fraction determined to generate the Fracmin Sketch (used for estimating union size);
// @param bins : the number of bins used
// @param t_max : The maximum number of elements per bin.
// @param f: When splitting a singular sequence, determines how big a bin can be, i.e. |Bin| = t_max * f
template <typename Hasher>
std::tuple<std::vector<std::vector<size_t>>, std::tuple<size_t,size_t,size_t>, std::vector<size_t>, bool> binning(const std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, Hasher>>& labMaps, 
                                            const std::vector<std::unordered_map<size_t,const std::vector<size_t>*>>& level_clusters,
                                            const std::vector<std::vector<std::uint64_t>>& fracmin_sketches,
                                            const double s, const size_t bins, const size_t t_max, const std::vector<double>& fcorrs,
                                            std::unordered_map<const std::vector<size_t>*, size_t>* top_size_cache = nullptr);
#include "templates/fast_construct_binning.tpp"

// @brief In order to Bin Merge Bins again, we want to filter the LSH Tree in order to easily access the new bins again.
// @param filtered_labMaps : This labMap represents the intersect of the relevant sequences (i.e. MergeBins) and the LSH tree. Sequences not in the Merge Bin will not be contained
// @param filtered_level_clusters : Like level_clusters but contains only the sequences present in the Merge Bin
// @param filtered_cluster_storage : Since the pointers would die if not saved, the new clusters will be stored as well.
template <typename Hasher>
struct LSH_Filtered{
    std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, Hasher>> filtered_labMaps;

    std::vector<std::unordered_map<size_t, const std::vector<size_t>*>> filtered_level_clusters;

    std::vector<std::vector<std::vector<size_t>>> filtered_cluster_storage;

};

// @brief Filter the LSH Forest such that only the relevant sequences are represented in a cluster
// @param labMaps : The Map that represents the Forest
// @param level_clusters : Maps for every Level the sequence to its cluster
// @param relevant_seqs : A subset of all sequences saying which are passing through the filter.
template <typename Hasher>
LSH_Filtered<Hasher> filter_LSH(const std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, Hasher>>& labMaps,
                                const std::vector<std::unordered_map<size_t, const std::vector<size_t>*>>& level_clusters,
                                const std::vector<size_t>& relevant_seqs);
#include "templates/fast_construct_filter_LSH.tpp"

// @brief Given a LSH Forest and the relevamt sequences, perform the binning Mechanism on ONLY these sequences
// @param labMaps : The clustering for the graph used.
// @param level_clusters : For every level and sequence, a pointer to the corresponding cluster in the higher level is returned.
// @param fracmin_sketches : The Set of Fracmin sketches corresponding to every sequence per ID
// @param relevant_seqs : The sequences which should be filtered.
// @param s : The Fraction determined to generate the Fracmin Sketch (used for estimating union size);
// @param bins : the number of bins used
// @param t_max : The maximum number of elements per bin.
// @param f: When splitting a singular sequence, determines how big a bin can be, i.e. |Bin| = t_max * f
template <typename Hasher>
std::tuple<std::vector<std::vector<size_t>>, std::tuple<size_t,size_t,size_t>, std::vector<size_t>, bool> binning_given_seqs(const std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, Hasher>>& labMaps, 
                                            const std::vector<std::unordered_map<size_t,const std::vector<size_t>*>>& level_clusters,
                                            const std::vector<std::vector<std::uint64_t>>& fracmin_sketches,
                                            const std::vector<size_t>& relevant_seqs,
                                            const double s, const size_t bins, const size_t t_max, const std::vector<double>& fcorrs,
                                            std::unordered_map<const std::vector<size_t>*, size_t>* top_size_cache = nullptr);
#include "templates/fast_construct_binning_given_seqs.tpp"


// @brief Given a number of bins and the higher bound of elements per bin, gather sequences, such that similar sequences are within the same bin and dissimilar sequences are in different bins.
// @param labMaps : The clustering for the graph used.
// @param level_clusters : For every level and sequence, a pointer to the corresponding cluster in the higher level is returned.
// @param fracmin_sketches : The Set of Fracmin sketches corresponding to every sequence per ID
// @param s : The Fraction determined to generate the Fracmin Sketch (used for estimating union size);
// @param bins : the number of bins used
// @param t_max : The maximum number of elements per bin.
// @param f: When splitting a singular sequence, determines how big a bin can be, i.e. |Bin| = t_max * f
// @note this is the core implementation. Calls above are just wrappers.
template <typename Hasher>
std::tuple<std::vector<std::vector<size_t>>, std::tuple<size_t,size_t,size_t>, std::vector<size_t>, bool> binning_core(const std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, Hasher>>& labMaps, 
                                            const std::vector<std::unordered_map<size_t,const std::vector<size_t>*>>& level_clusters,
                                            const std::vector<std::vector<std::uint64_t>>& fracmin_sketches,
                                            const double s, const size_t bins, const size_t t_max, const std::vector<double>& fcorrs,
                                            std::unordered_map<const std::vector<size_t>*, size_t>* top_size_cache = nullptr);
#include "templates/fast_construct_binning_core.tpp"

// @brief given the amount of Hash functions used and the wanted false positive rate (fpr), compute for 100 bin sizes the correction factor.
// @param fpr : The desired false positive rate.
// @param h : The amount of Hash functions that will be used.
std::vector<double> compute_fcorrs(double fpr, size_t h);

// @brief given the trackfill information of an IBF and where merging starts, compute the estimated merge size of each bin.
// @param trackfill : The amount of elements in each bin
// @param merge_start : Index, where merging starts (inclusive)
size_t merge_average(const std::vector<size_t>& track_fill, const size_t merge_start);

// @brief given the trackfill information of an IBF and where splitting starts and ends, compute the estimated splitting size of each bin.
// @param track_fill : The amount of elements in each bin
// @param split_start : Index, where splitting starts (inclusive)
// @param split_end : Index, where splitting ends (exclusive)
size_t splitting_average(const std::vector<size_t>& track_fill, const size_t split_start, const size_t split_end);

// @brief Given every parameter needed, construct the HIBF completely. The bin size for every IBF is refined p+1 times. 
// @param signatures : contains both the One Permutation Hash and the Fracmin Hash signatures
// @param levels : parameters used for the LSH clustering
// @param s : The Fraction determined to generate the Fracmin Sketch (used for estimating union size);
// @param bins : the number of bins used
// @param fpr : The false positive rate aimed for. It will be used to compute fcorr.
// @param h : The amount of hash functions that will be used. Used for computing fcorr.
// @param p : The refinement parameter for the bin size. The IBF will be constructed p+1 times. If p == 0, it will just be (Union + Sum)/(2*bins*s)
// @param max_level : Limits how many the HIBF is allowed to have.
using IBF = std::vector<std::vector<size_t>>;
template <typename Hasher>
std::tuple<std::vector<std::vector<IBF>>, std::vector<std::vector<std::tuple<size_t,size_t,size_t>>>, std::unordered_map<size_t, std::vector<std::tuple<size_t,size_t,size_t>>>, std::vector<std::vector<size_t>>, std::vector<std::vector<std::pair<size_t,size_t>>>> generate_hibf(const std::pair<std::vector<std::vector<std::uint64_t>>, std::vector<std::vector<std::uint64_t>>>& signatures,
                                            const std::vector<std::pair<size_t,size_t>>& levels,
                                            const double s, const double fpr, const size_t h, const size_t p);
#include "templates/fast_construct_generate_hibf.tpp"

void write_linkage(std::ostream& out,
                    const std::unordered_map<size_t, std::string>& seq_to_path);

void write_config(std::ofstream& out,
                  std::filesystem::path const& dir_path,
                  std::uint8_t q,
                  std::uint32_t w,
                  double fpr,
                  std::uint8_t hash_funcs,
                  size_t user_bins);

void write_header(std::ostream& out, const std::vector<std::vector<IBF>>& hibf_levels, 
                  const std::vector<std::vector<size_t>>& max_bin_ids, 
                  const std::vector<std::vector<std::pair<size_t,size_t>>>& parents);

void write_content(std::ofstream& out, 
                   const std::unordered_map<size_t, std::vector<std::tuple<size_t,size_t,size_t>>>& seq_layout);