#include "../fast_construct_graph.h"

template <typename Hasher>
std::tuple<std::vector<std::vector<size_t>>, std::tuple<size_t,size_t,size_t>, std::vector<size_t>, bool> binning_given_seqs(const std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, Hasher>>& labMaps, 
                                            const std::vector<std::unordered_map<size_t,const std::vector<size_t>*>>& level_clusters,
                                            const std::vector<std::vector<std::uint64_t>>& fracmin_sketches,
                                            const std::vector<size_t>& relevant_seqs,
                                            const double s, const size_t bins, const size_t t_max, const std::vector<double>& fcorrs,
                                            std::unordered_map<const std::vector<size_t>*, size_t>* top_size_cache){
                                            

    LSH_Filtered<Hasher> filtered = filter_LSH(labMaps, level_clusters, relevant_seqs);
    return binning_core(filtered.filtered_labMaps, filtered.filtered_level_clusters, fracmin_sketches, s, bins, t_max, fcorrs, top_size_cache);
}
