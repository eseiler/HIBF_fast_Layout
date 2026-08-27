#include "../fast_construct_graph.h"
#include <limits>

template <typename Hasher>
std::tuple<std::vector<std::vector<size_t>>, std::tuple<size_t,size_t,size_t>, std::vector<size_t>, bool> binning_core(
    const std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, Hasher>>& labMaps,
    const std::vector<std::unordered_map<size_t, const std::vector<size_t>*>>& level_clusters,
    const std::vector<std::vector<std::uint64_t>>& fracmin_sketches,
    const double s,
    const size_t bins,
    const size_t t_max,
    const std::vector<double>& fcorrs,
    std::unordered_map<const std::vector<size_t>*, size_t>* top_size_cache) {

    size_t deepest_lvl = labMaps.size() - 1;
    size_t par_lvl = (deepest_lvl > 0) ? deepest_lvl - 1 : deepest_lvl;
    size_t bin = 0;
    size_t merge_bins = 0;
    size_t split_bins = 0;
    size_t tot_seqs = 0;

    bool overflow = false;

    std::vector<std::vector<size_t>> res(bins);

    if (bins == 0) return {res, {0,0,0}, {}, false};

    std::vector<size_t> track_fill(bins, 0);

    std::vector<std::unordered_set<std::uint64_t>> bin_sketches(bins);

    std::vector<const std::vector<size_t>*> candidates;
    std::vector<const std::vector<size_t>*> top_clusters;
    std::vector<const std::vector<size_t>*> small_clusters;

    std::unordered_set<size_t> reserved_small_seqs;

    std::vector<const std::vector<std::uint64_t>*> supersketch;

    std::vector<size_t> bin_lvls(bins, deepest_lvl);
    std::unordered_set<size_t> binned;

    std::unordered_map<const std::vector<size_t>*, size_t> parent_size;
    std::unordered_map<const std::vector<size_t>*, size_t> top_sizes;
    std::unordered_map<const std::vector<size_t>*, size_t> fair_bins;
    std::unordered_map<const std::vector<size_t>*, std::unordered_set<size_t>> bins_used_by_cluster;

    std::vector<std::unordered_set<const std::vector<size_t>*>> used_clusters_per_level(deepest_lvl + 1);
    std::vector<std::unordered_map<const std::vector<size_t>*, size_t>> bin_for_cluster(deepest_lvl + 1);
    std::vector<std::unordered_map<const std::vector<size_t>*, size_t>> itteration_progress(deepest_lvl + 1);

    std::queue<size_t> valid_bins;
    std::priority_queue<std::pair<size_t,size_t>, std::vector<std::pair<size_t,size_t>>, std::greater<std::pair<size_t,size_t>>> still_fillable;

    /// @note we declare some helper functions here, which come in handy later.

    auto update_bin = [&](size_t lvl, const std::vector<size_t>* cluster, size_t b) {bin_for_cluster[lvl][cluster] = b;};

    auto try_insert_sequence = [&](size_t seq, size_t b, bool force) {
        if (b >= bins) return false;
        // A bin holding more than one sequence becomes a merge bin, i.e. a whole IBF on the level
        // below, and that child gets the same number of technical bins as this one. Letting a bin
        // grow past that guarantees the child has to merge again.
        if (!force && res[b].size() >= bins) return false;

        std::unordered_set<std::uint64_t>& sketch = bin_sketches[b];
        const std::vector<std::uint64_t>& seq_sketch = fracmin_sketches[seq];
        const bool may_reject = !force && !res[b].empty();

        // Seeding, merging and the fallback all offer the same sequence to many bins in turn, so
        // most calls end in a rejection. Both checks below make a rejection cheap: the union can
        // never be smaller than either side, which rejects a full bin without touching the hash
        // table, and the loop stops as soon as the union is known to exceed t_max instead of
        // probing the whole sketch and throwing the result away. Neither changes which insertions
        // are accepted.
        if (may_reject && static_cast<size_t>(static_cast<double>(std::max(sketch.size(), seq_sketch.size())) / s) > t_max)
            return false;

        std::vector<std::uint64_t> new_elems;
        for (std::uint64_t elem : seq_sketch) {
            if (sketch.count(elem)) continue;
            new_elems.push_back(elem);
            if (may_reject && static_cast<size_t>(static_cast<double>(sketch.size() + new_elems.size()) / s) > t_max)
                return false;
        }

        size_t potential_size = static_cast<size_t>(static_cast<double>(sketch.size() + new_elems.size()) / s);
        if (may_reject && potential_size > t_max) return false;
        for (std::uint64_t new_elem : new_elems) sketch.insert(new_elem);
        res[b].push_back(seq);
        binned.insert(seq);
        track_fill[b] = potential_size;
        return true;
    };

    auto get_top_parent = [&](size_t seq) -> const std::vector<size_t>* {
        auto it = level_clusters[0].find(seq);
        return (it != level_clusters[0].end()) ? it->second : nullptr;
    };

    auto cluster_has_budget_for_bin = [&](const std::vector<size_t>* top_parent, size_t b) -> bool {
        if (!top_parent) return true;
        auto fb_it = fair_bins.find(top_parent);
        if (fb_it == fair_bins.end()) return true;
        auto used_it = bins_used_by_cluster.find(top_parent);
        if (used_it != bins_used_by_cluster.end() && used_it->second.count(b)) return true;
        size_t already_used = (used_it != bins_used_by_cluster.end()) ? used_it->second.size() : 0;
        return already_used < fb_it->second;
    };

    auto register_bin_for_cluster = [&](const std::vector<size_t>* top_parent, size_t b) {
        if (!top_parent) return;
        if (!fair_bins.count(top_parent)) return;
        bins_used_by_cluster[top_parent].insert(b);
    };


    /// @note this function is used when distributing the greater bottom-level clusters across multiple bins.

    auto split_large_component_balanced = [&](const std::vector<size_t>& component) {
        if (component.empty()) return;

        const std::vector<size_t>* top_parent = get_top_parent(component.front());

        size_t already_used = (top_parent && bins_used_by_cluster.count(top_parent)) ? bins_used_by_cluster[top_parent].size() : 0;
        size_t budget = (top_parent && fair_bins.count(top_parent)) ? fair_bins[top_parent] : std::numeric_limits<size_t>::max();
        size_t new_bin_budget = (budget > already_used) ? (budget - already_used) : 0;
        std::vector<size_t> available_bins;
        for (size_t b = bin; b < bins && available_bins.size() < new_bin_budget; ++b)if (res[b].empty()) available_bins.push_back(b);
        if (available_bins.empty()) return;
        std::vector<size_t> seqs_to_distribute;
        seqs_to_distribute.reserve(component.size());
        for (size_t seq : component) if (!binned.count(seq) && !reserved_small_seqs.count(seq)) seqs_to_distribute.push_back(seq);

        size_t rr_index = 0;
        for (size_t seq : seqs_to_distribute) {
            size_t num_attempts = 0;

            while (num_attempts < available_bins.size()) {
                size_t target_bin = available_bins[rr_index % available_bins.size()];
                rr_index += 1;
                num_attempts += 1;
                bool force = (num_attempts == available_bins.size());
                if (try_insert_sequence(seq, target_bin, force)) {
                    register_bin_for_cluster(top_parent, target_bin);
                    for (size_t lvl_ = 0; lvl_ <= deepest_lvl; ++lvl_) {
                        auto iterator = level_clusters[lvl_].find(seq);
                        if (iterator != level_clusters[lvl_].end()) {
                            used_clusters_per_level[lvl_].insert(iterator->second);
                            update_bin(lvl_, iterator->second, target_bin);
                        }
                    }
                    break;
                }
            }
        }
        while (bin < bins && !res[bin].empty()) bin += 1;};

    /// @note this part is the splitting. If a sequence sketch is too big (sketch_size > t_max), it must be split across mutliple bins, scalled with f (from fcorrs)

    for (auto& [component, node] : labMaps[0]) {
        tot_seqs += component.size();
        for (size_t seq : component) {
            if (bin >= bins) return {res, {0,0,0}, track_fill, true};
            size_t est_sketch_size = fracmin_sketches[seq].size()/s;
            if (est_sketch_size >= t_max) {
                split_bins += 1;
                size_t split_bins_ = 100;
                size_t already_split = 0;
                size_t low = 1;
                size_t high = 100;

                while (low <= high) {
                    size_t mid = low + (high - low)/2;
                    double capacity = static_cast<double>(mid) * static_cast<double>(t_max) * fcorrs[mid];
                    if (static_cast<size_t>(capacity) >= est_sketch_size){
                        split_bins_ = mid;
                        if (mid == 0) break;
                        high = mid - 1;
                    } 
                    else low = mid + 1; 
                }
                split_bins_ = (split_bins_ > bins - bin) ? (bins - bin) : split_bins_;
                double f = fcorrs[split_bins_];

                for (std::uint64_t elem : fracmin_sketches[seq]) {
                    std::unordered_set<std::uint64_t>& sketch = bin_sketches[bin + already_split];
                    size_t sketch_size = track_fill[bin + already_split];
                    size_t potential_sketch_size = sketch_size + static_cast<size_t>(static_cast<double>(1) / s);
                    if (potential_sketch_size > t_max * f && already_split + 1 < split_bins_) {
                        res[bin + already_split].push_back(seq);
                        already_split += 1;
                        split_bins += 1;
                    }
                    bin_sketches[bin + already_split].insert(elem);
                    track_fill[bin + already_split] += static_cast<size_t>(static_cast<double>(1) / s);
                }
                binned.insert(seq);
                bin += (already_split + 1);
                res[bin - 1].push_back(seq);
            }
        }
        if (overflow) break;
    }

    if (bin >= bins && binned.size() != tot_seqs) return {res, {0,0,0}, track_fill, true};
    for (auto& [component, node] : labMaps[0]) top_clusters.push_back(&component);
    top_sizes.reserve(top_clusters.size());
    for (const std::vector<size_t>* clust : top_clusters) {
        // The estimated size of a cluster does not depend on t_max, but binning_core is called
        // once per t_max while the refinement converges. The caller can hand in a cache so these
        // unions are computed once per IBF instead of once per t_max.
        size_t estimated_size = 0;
        bool have_size = false;
        if (top_size_cache) {
            auto cached = top_size_cache->find(clust);
            if (cached != top_size_cache->end()) { estimated_size = cached->second; have_size = true; }
        }
        if (!have_size) {
            std::vector<const std::vector<std::uint64_t>*> subset;
            subset.reserve(clust->size());
            for (size_t seq : *clust) subset.push_back(&fracmin_sketches[seq]);
            size_t union_size = get_union_size_ptr(subset);
            estimated_size = static_cast<std::uint64_t>(static_cast<double>(union_size) / s);
            if (top_size_cache) (*top_size_cache)[clust] = estimated_size;
        }
        top_sizes[clust] = estimated_size;
        if (estimated_size < t_max) for (size_t seq : *clust) supersketch.push_back(&fracmin_sketches[seq]);
    }

    /// @note We compute the "Share" for every top-level cluster. We try to estimate how many bins it is allowed to occupy
    /// This is useful for seeding and merging, since we can say if the cluster should get another representative.

    struct FairShare {
        const std::vector<size_t>* clust;
        size_t weight;
        size_t floor_bins;
        double remainder;
    };

    size_t bins_left_for_top = (bins > bin) ? (bins - bin) : 0;
    std::vector<FairShare> shares;
    shares.reserve(top_clusters.size());
    size_t total_weight = 0;

    for (const std::vector<size_t>* clust : top_clusters) {
        size_t remaining = 0;
        for (size_t seq : *clust)if (!binned.count(seq)) remaining += 1;
        size_t weight = (remaining == 0) ? 0 : std::max<size_t>(top_sizes[clust], 1);

        shares.push_back({clust, weight, 0, 0.0});
        total_weight += weight;

        fair_bins[clust] = 0;
    }

    if (total_weight > 0 && bins_left_for_top > 0) {
        size_t distributed = 0;
        for (FairShare& share : shares) {
            double exact = (static_cast<double>(share.weight) / static_cast<double>(total_weight)) * static_cast<double>(bins_left_for_top);

            share.floor_bins = static_cast<size_t>(exact);
            share.remainder = exact - static_cast<double>(share.floor_bins);

            fair_bins[share.clust] = share.floor_bins;
            distributed += share.floor_bins;
        }
        size_t leftover = (bins_left_for_top > distributed) ? (bins_left_for_top - distributed) : 0;
        std::sort(shares.begin(), shares.end(), [](const FairShare& a, const FairShare& b) {
            if (a.remainder != b.remainder) return a.remainder > b.remainder;
            return a.weight > b.weight;
        });
        for (size_t i = 0; i < leftover && i < shares.size(); ++i) fair_bins[shares[i].clust] += 1;
    }
    


    /// @note this is again a bit of setup; we sort the top level clusters and declare the clusters too small for seeding.
    /// We save these clusters for later in the Merging phase. 
    /// We already estimate how many bins these clusters might need.
    std::sort(top_clusters.begin(), top_clusters.end(), [&top_sizes](const std::vector<size_t>* a, const std::vector<size_t>* b) {
        size_t size_a = top_sizes[a];
        size_t size_b = top_sizes[b];
        if (size_a != size_b) return size_a < size_b;
        return a->front() < b->front();
    });

    for (auto it = top_clusters.rbegin(); it != top_clusters.rend(); ++it) {
        const std::vector<size_t>* clust = *it;
        if (top_sizes[clust] < t_max) {
            small_clusters.push_back(clust);
            used_clusters_per_level[0].insert(clust);

            for (size_t seq : *clust) reserved_small_seqs.insert(seq);
        }
    }

    merge_bins = bin;

    size_t allowed_merge = static_cast<size_t>(static_cast<double>(get_union_size_ptr(supersketch)) / (s * static_cast<double>(t_max)));

        /// @note This is the Late Merging
    /// After using Every cluster, we are left with our estimate usage of the Merge Bins.
    /// We start by entering the "biggest" small cluster and itterate over every usable bin. The small clusters left are used in the Fallback Mechanism.
    
    // Bins below merge_bins hold pieces of split sequences. record_bins expects those to stay a
    // contiguous run of one-sequence bins, so nothing else may be placed there.
    size_t seeding_max_bin = std::max(merge_bins, (bins > allowed_merge) ? (bins - allowed_merge) : bins);
    size_t empty_start = seeding_max_bin;
    size_t empty_end = bins;

    if (empty_start < empty_end) {
        size_t available_merge_bins = empty_end - empty_start;
        size_t merge_rr_idx = 0;

        for (const std::vector<size_t>* cluster : small_clusters) {
            for (size_t seq : *cluster) {
                if (binned.count(seq)) continue;
                size_t attempts = 0;
                while (attempts < available_merge_bins) {
                    size_t curr_bin = empty_start + (merge_rr_idx % available_merge_bins);
                    merge_rr_idx++;
                    attempts++;

                    if (cluster_has_budget_for_bin(cluster, curr_bin)) {
                        if (try_insert_sequence(seq, curr_bin, false)) {
                            register_bin_for_cluster(cluster, curr_bin);

                            for (size_t lvl = 0; lvl <= deepest_lvl; lvl++) {
                                auto it = level_clusters[lvl].find(seq);
                                if (it != level_clusters[lvl].end()) {
                                    used_clusters_per_level[lvl].insert(it->second);
                                    update_bin(lvl, it->second, curr_bin);
                                }
                            }
                            break;
                        }
                    }
                }
            }
        }
    }

    /// @note This is the cluster Splitting
    /// If a cluster is too big on the bottom level, we call the helper function above and Round-Robin distribute them.

    while (bin < bins && !res[bin].empty()) bin += 1;

    size_t start = bin;

    for (auto& [component, node] : labMaps[deepest_lvl]) {
        bool all_binned = std::all_of(component.begin(), component.end(), [&](size_t s) {
            return binned.count(s) != 0 || reserved_small_seqs.count(s) != 0;
        });
        if (all_binned)continue;
        std::vector<const std::vector<std::uint64_t>*> subset;
        subset.reserve(component.size());
        for (size_t seq : component) {
            if (!binned.count(seq) && !reserved_small_seqs.count(seq)) {
                subset.push_back(&fracmin_sketches[seq]);
            }
        }
        size_t union_size = get_union_size_ptr(subset);
        size_t estimator = static_cast<std::uint64_t>(static_cast<double>(union_size) / s);
        if (estimator <= t_max) {
            candidates.push_back(&component);
            continue;
        }
        split_large_component_balanced(component);
        candidates.push_back(&component);
    }

    /// @note This is the seeding algorithm
    /// We start by the smallest clusters and generate a top down diversity by always checking if the top level cluster already was used. If all are used, we descend one level.

    parent_size.reserve(candidates.size());

    for (const std::vector<size_t>* cand : candidates) {
        auto it = level_clusters[par_lvl].find(cand->front());
        parent_size[cand] = (it != level_clusters[par_lvl].end()) ? it->second->size() : cand->size();
    }

    std::sort(candidates.begin(), candidates.end(), [&parent_size](const std::vector<size_t>* a, const std::vector<size_t>* b) {
        if (a->size() != b->size()) return a->size() < b->size();
        size_t size_a = parent_size[a];
        size_t size_b = parent_size[b];
        if (size_a != size_b) return size_a > size_b;
        return a->front() < b->front();
    });

    
    size_t seeding_rr_bin = bin;
    // Late merging reserves the tail of the bin range, and splitting consumes the head, so the
    // seeding window can be empty. Subtracting unsigned without the guard wraps to ~2^64 and
    // turns the loop below into an effectively infinite spin.
    size_t num_seeding_bins = (seeding_max_bin > start) ? (seeding_max_bin - start) : 0;

    for (size_t cl = 0; cl < candidates.size(); cl++) {
        const std::vector<size_t>* cluster = candidates[cl];
        size_t seq_idx = 0;
        while (seq_idx < cluster->size() && (binned.count((*cluster)[seq_idx]) || reserved_small_seqs.count((*cluster)[seq_idx]))) seq_idx += 1;

        if (seq_idx >= cluster->size()) continue;
        size_t seed_seq = (*cluster)[seq_idx];
        const std::vector<size_t>* top_parent = get_top_parent(seed_seq);

        size_t attempts = 0;
        while (attempts < num_seeding_bins) {
            size_t empty_search_count = 0;
            while (empty_search_count < num_seeding_bins && !res[seeding_rr_bin].empty()) {
                seeding_rr_bin += 1;
                if (seeding_rr_bin >= seeding_max_bin) seeding_rr_bin = start;
                empty_search_count++;
            }
            if (seeding_rr_bin >= seeding_max_bin) seeding_rr_bin = start;

            if (cluster_has_budget_for_bin(top_parent, seeding_rr_bin)) {
                if (try_insert_sequence(seed_seq, seeding_rr_bin, false)) {
                    register_bin_for_cluster(top_parent, seeding_rr_bin);

                    for (size_t lvl_ = 0; lvl_ <= deepest_lvl; ++lvl_) {
                        auto iterator = level_clusters[lvl_].find(seed_seq);
                        if (iterator != level_clusters[lvl_].end()) {
                            used_clusters_per_level[lvl_].insert(iterator->second);
                            update_bin(lvl_, iterator->second, seeding_rr_bin);
                        }
                    }
                    seeding_rr_bin += 1;
                    if (seeding_rr_bin >= seeding_max_bin) seeding_rr_bin = start;
                    break;
                }
            }
            seeding_rr_bin += 1;
            if (seeding_rr_bin >= seeding_max_bin) seeding_rr_bin = start;
            attempts += 1;
        }
    }


    /// @note This part is the climbing mechanism. We take every representative and climb its cluster upwards by recursively going up one level after finishing. 
    /// This goes to lvl 0 (root), but elements are not taken from the mother clusters.

    for (size_t b = start; b < seeding_max_bin; b++) if (!res[b].empty() && track_fill[b] < t_max)valid_bins.push(b);

    while (!valid_bins.empty()) {
        size_t b = valid_bins.front();
        valid_bins.pop();
        size_t repr = res[b][0];
        size_t curr_lvl = bin_lvls[b];
        if (curr_lvl == 0) continue;
        auto it = level_clusters[curr_lvl].find(repr);
        if (it == level_clusters[curr_lvl].end()) {
            if (curr_lvl > 0) {
                bin_lvls[b] -= 1;
                valid_bins.push(b);
            }
            continue;
        }

        const std::vector<size_t>& cluster = *(it->second);
        update_bin(curr_lvl, it->second, b);
        size_t& start_it = itteration_progress[curr_lvl][it->second];
        for (; start_it < cluster.size(); start_it++) {
            if (track_fill[b] >= t_max) break;

            size_t seq = cluster[start_it];
            if (binned.count(seq) || reserved_small_seqs.count(seq)) continue;
            if (!try_insert_sequence(seq, b, false)) break;  
        }

        if (track_fill[b] < t_max && curr_lvl > 0) {
            bin_lvls[b] -= 1;
            valid_bins.push(b);
        }
    }


    /// @note The Fallback mechanism.
    /// For every not entered sequence, we check whether there is already a cluster once checked with this element. If so, we can enter this element in the bin corresponding to the cluster.
    /// If there is no such cluster, take the emptiest bin and enter the sequence as a "new seed" and mark the cluster/s as seen for later.

    for (size_t b = start; b < bins; b++) if (track_fill[b] < t_max) still_fillable.push({track_fill[b], b});

    for (auto& [seq, _] : level_clusters[0]) {
        if (binned.count(seq)) continue;
        bool entered = false;
        for (size_t lvl = deepest_lvl; lvl > 0; --lvl) {
            auto it = level_clusters[lvl].find(seq);
            if (it == level_clusters[lvl].end()) continue;

            auto iterator = bin_for_cluster[lvl].find(it->second);
            if (iterator == bin_for_cluster[lvl].end()) continue;

            size_t b = iterator->second;
            if (b < start) continue;

            if (!try_insert_sequence(seq, b, false)) continue;
            bool overshoot = track_fill[b] > t_max;
            if (!overshoot) {
                for (size_t lvl_ = 0; lvl_ <= deepest_lvl; lvl_++) {
                    auto it_ = level_clusters[lvl_].find(seq);
                    if (it_ != level_clusters[lvl_].end()) update_bin(lvl_, it_->second, b);
                }
            }
            entered = true;
            break;
        }

        if (!entered) {
            std::vector<std::pair<size_t,size_t>> skipped;
            while (!still_fillable.empty()) {
                auto [fill, b] = still_fillable.top();
                still_fillable.pop();
                if (track_fill[b] >= t_max) continue;
                if (try_insert_sequence(seq, b, false)) {
                    entered = true;
                    bool overshoot = track_fill[b] > t_max;
                    if (!overshoot) {
                        for (size_t lvl_ = deepest_lvl; lvl_ > 0; lvl_--) {
                            auto it_ = level_clusters[lvl_].find(seq);
                            if (it_ != level_clusters[lvl_].end()) update_bin(lvl_, it_->second, b);
                        }
                    }
                    if (track_fill[b] < t_max) still_fillable.push({track_fill[b], b});
                    break;
                }
                skipped.push_back({track_fill[b], b});
            }
            for (auto& s : skipped) still_fillable.push(s);
            if (!entered) {
                // Late merging can claim every bin from `start` on, leaving start == bins. Indexing
                // track_fill with it reads out of bounds and try_insert_sequence rejects b >= bins,
                // which used to drop the sequence from the layout without a trace. Retry over every
                // bin that is not part of the split region; if splitting claimed all of them there
                // is genuinely no room and this t_max has overflowed.
                if (merge_bins >= bins) return {res, {0,0,0}, track_fill, true};
                size_t search_start = (start < bins) ? start : merge_bins;
                size_t best_b = search_start;
                size_t best_fill = track_fill[search_start];

                for (size_t b = search_start + 1; b < bins; b++) {
                    if (track_fill[b] < best_fill) {
                        best_fill = track_fill[b];
                        best_b = b;
                    }
                }
                try_insert_sequence(seq, best_b, true);
                bool overshoot = track_fill[best_b] > t_max;
                if (!overshoot) {
                    for (size_t lvl_ = deepest_lvl; lvl_ > 0; lvl_--) {
                        auto it_ = level_clusters[lvl_].find(seq);
                        if (it_ != level_clusters[lvl_].end()) update_bin(lvl_, it_->second, best_b);
                    }
                }
                if (track_fill[best_b] < t_max) still_fillable.push({track_fill[best_b], best_b});
            }
        }
    }


    /// @note This is preperation for what will be returned.

    std::vector<size_t> permutation(res.size());
    std::iota(permutation.begin(), permutation.end(), 0);
    std::sort(permutation.begin(), permutation.end(), [&res](size_t a, const size_t b) {
        if (res[a].size() != res[b].size()) return res[a].size() < res[b].size();
        if (res[a].empty()) return false;
        if (res[b].empty())return true;
        return res[a].front() < res[b].front();
    });

    std::vector<std::vector<size_t>> sorted_res(res.size());
    std::vector<size_t> sorted_trackfill(res.size());

    for (size_t k = 0; k < permutation.size(); k++) {
        sorted_res[k] = std::move(res[permutation[k]]);
        sorted_trackfill[k] = track_fill[permutation[k]];
    }

    res = std::move(sorted_res);
    track_fill = std::move(sorted_trackfill);

    auto merge_it = std::partition_point(res.begin(), res.end(), [](const std::vector<size_t>& bin) {return bin.size() < 2;});

    auto split_it = std::partition_point(res.begin(), res.end(), [](const std::vector<size_t>& bin) {return bin.size() == 0;});

    size_t merge_start = std::distance(res.begin(), merge_it);
    size_t split_start = std::distance(res.begin(), split_it);

    return {res, std::make_tuple(split_start, split_bins, merge_start), track_fill, overflow};
}