#include "../fast_construct_bin.h"
#include <future>
#include <thread>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <tuple>
#include <vector>
#include <unordered_map>
#include <limits>

template <typename Hasher>
std::tuple<std::vector<std::vector<IBF>>, 
           std::vector<std::vector<std::tuple<size_t,size_t,size_t>>>, 
           std::unordered_map<size_t, std::vector<std::tuple<size_t,size_t,size_t>>>, 
           std::vector<std::vector<size_t>>, 
           std::vector<std::vector<std::pair<size_t,size_t>>>> 
generate_hibf(const std::tuple<std::vector<std::vector<std::uint64_t>>, std::vector<std::vector<std::uint64_t>>, std::unordered_map<size_t, std::string>>& signatures,
              const std::vector<std::pair<size_t,size_t>>& levels,
              const double s, const double fpr, const size_t h, const size_t p, const size_t threads) {
    
    using BinResult = std::tuple<std::vector<std::vector<size_t>>, std::tuple<size_t,size_t,size_t>, std::vector<size_t>, bool>;
    using ChildResult = BinResult;

    size_t max = std::numeric_limits<size_t>::max();
    const std::vector<std::vector<std::uint64_t>>& oph_sigs = std::get<0>(signatures);
    const std::vector<std::vector<std::uint64_t>>& fracmin_sigs = std::get<1>(signatures);
    const std::vector<double> fcorrs = compute_fcorrs(fpr, h);
    
    std::unordered_map<size_t, std::vector<std::tuple<size_t,size_t,size_t>>> seq_layout;
    std::vector<std::vector<size_t>> max_bin_ids;
    std::vector<std::vector<std::pair<size_t,size_t>>> parents;
    
    lemon::ListGraph graph;
    std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, Hasher>> labMaps = generate_all<Hasher>(oph_sigs, levels, graph);
    std::vector<std::unordered_map<size_t, const std::vector<size_t>*>> clusts = get_clusters(labMaps);
    
    // Estimated cost of a layout, in k-mers. An IBF allocates all of its technical bins at the size
    // of the fullest one, and every merge bin turns into a whole additional IBF one level down.
    auto layout_cost = [](const BinResult& r, size_t sub_bins) {
        const std::vector<std::vector<size_t>>& result = std::get<0>(r);
        const std::vector<size_t>& fill = std::get<2>(r);
        const size_t merge_start = std::get<2>(std::get<1>(r));
        size_t max_fill = 0;
        for (size_t f : fill) max_fill = std::max(max_fill, f);
        size_t child_content = 0;
        for (size_t b = merge_start; b < result.size(); b++) if (!result[b].empty()) child_content += fill[b];
        return static_cast<double>(sub_bins) * static_cast<double>(max_fill) + static_cast<double>(child_content);
    };

    // A merge bin holding every sequence the IBF was given would produce a child IBF with the
    // identical input, and generate_hibf would recurse on it forever. Requiring every merge bin to
    // be strictly smaller than the input makes the recursion depth finite by construction.
    auto makes_progress = [](const BinResult& r, size_t n_seqs) {
        const std::vector<std::vector<size_t>>& result = std::get<0>(r);
        const size_t merge_start = std::get<2>(std::get<1>(r));
        for (size_t b = merge_start; b < result.size(); b++) if (result[b].size() >= n_seqs) return false;
        return true;
    };

    // Last resort for when no t_max produced a usable layout: one user bin per technical bin for as
    // long as they fit, round robin afterwards. Always a valid layout, and never puts every sequence
    // into a single bin, so it always makes progress.
    auto fallback_layout = [&](const std::vector<size_t>& seqs, size_t sub_bins) -> BinResult {
        std::vector<std::vector<size_t>> res(sub_bins);
        std::vector<size_t> fill(sub_bins, 0);
        for (size_t i = 0; i < seqs.size(); i++) res[i % sub_bins].push_back(seqs[i]);

        for (size_t b = 0; b < sub_bins; b++) {
            if (res[b].empty()) continue;
            std::vector<const std::vector<std::uint64_t>*> ptrs;
            ptrs.reserve(res[b].size());
            for (size_t seq : res[b]) ptrs.push_back(&fracmin_sigs[seq]);
            fill[b] = static_cast<size_t>(static_cast<double>(get_union_size_ptr(ptrs)) / s);
        }

        // binning_core orders its bins empty first, then split bins, then merge bins; keep that.
        std::vector<size_t> permutation(sub_bins);
        std::iota(permutation.begin(), permutation.end(), 0);
        std::sort(permutation.begin(), permutation.end(), [&res](size_t a, size_t b) {
            if (res[a].size() != res[b].size()) return res[a].size() < res[b].size();
            if (res[a].empty()) return false;
            if (res[b].empty()) return true;
            return res[a].front() < res[b].front();
        });

        std::vector<std::vector<size_t>> sorted_res(sub_bins);
        std::vector<size_t> sorted_fill(sub_bins);
        for (size_t i = 0; i < sub_bins; i++) {
            sorted_res[i] = std::move(res[permutation[i]]);
            sorted_fill[i] = fill[permutation[i]];
        }

        size_t split_start = 0;
        while (split_start < sub_bins && sorted_res[split_start].empty()) split_start += 1;
        size_t merge_start = split_start;
        while (merge_start < sub_bins && sorted_res[merge_start].size() < 2) merge_start += 1;

        return {std::move(sorted_res), std::make_tuple(split_start, size_t{0}, merge_start), std::move(sorted_fill), false};
    };

auto refine_and_bin = [&](const std::vector<size_t>& seqs, size_t sub_bins, size_t lower, size_t upper, bool all_seqs) -> BinResult {
    if (sub_bins == 0) return BinResult{{}, std::make_tuple(size_t{0}, size_t{0}, size_t{0}), {}, false};

    const size_t n_seqs = all_seqs ? fracmin_sigs.size() : seqs.size();

    // Every sequence at least this big gets split across several technical bins, so above this t_max
    // nothing is split at all and binning cannot run out of bins. Keeping it inside the search
    // interval guarantees a non-overflowing t_max is always reachable, whatever the interval derived
    // from the content happens to be.
    size_t max_single = 0;
    if (all_seqs)
        for (const std::vector<std::uint64_t>& sketch : fracmin_sigs)
            max_single = std::max(max_single, static_cast<size_t>(static_cast<double>(sketch.size()) / s));
    else
        for (size_t seq : seqs)
            max_single = std::max(max_single, static_cast<size_t>(static_cast<double>(fracmin_sigs[seq].size()) / s));
    const size_t feasible_t_max = max_single + max_single/16 + 1;

    size_t curr_t_max = (lower + upper)/(2*sub_bins*s);
    size_t curr_upper = std::max(curr_t_max * 2, feasible_t_max);
    size_t curr_lower = 0;
    size_t old_t_max = 0;

    // Filtering the LSH forest down to this IBF's sequences depends only on `seqs`, and the
    // estimated size of a cluster depends only on its contents. Both used to be redone for every
    // t_max the refinement tried; hoisting them here makes the loop cost one binning pass each.
    LSH_Filtered<Hasher> filtered;
    const std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, Hasher>>* use_maps = &labMaps;
    const std::vector<std::unordered_map<size_t, const std::vector<size_t>*>>* use_clusts = &clusts;
    if (!all_seqs) {
        filtered = filter_LSH<Hasher>(labMaps, clusts, seqs);
        use_maps = &filtered.filtered_labMaps;
        use_clusts = &filtered.filtered_level_clusters;
    }
    // Keyed by cluster pointer, so it must not outlive `filtered`.
    std::unordered_map<const std::vector<size_t>*, size_t> top_size_cache;

    auto run = [&](size_t t_max) {
        return binning_core(*use_maps, *use_clusts, fracmin_sigs, s, sub_bins, t_max, fcorrs, &top_size_cache);
    };

    BinResult best;
    bool valid = false;
    double best_cost = 0.0;

    auto consider = [&](const BinResult& r) {
        if (std::get<3>(r)) return;                    // overflowed, so this is not a layout
        if (!makes_progress(r, n_seqs)) return;        // would hand a child its own input
        const double cost = layout_cost(r, sub_bins);
        if (!valid || cost < best_cost) { best = r; best_cost = cost; valid = true; }
    };

    for (size_t it = 0; it <= p;) {
        if (curr_t_max == old_t_max) break;

        BinResult res = run(curr_t_max);
        const bool overflow = std::get<3>(res);
        consider(res);

        // Overflow means t_max was too small to fit the split sequences. Anything else means the
        // layout fits, and a smaller t_max makes every technical bin cheaper, so the search always
        // walks towards the smallest t_max that still yields a layout.
        if (overflow) curr_lower = curr_t_max;
        else          curr_upper = curr_t_max;

        old_t_max = curr_t_max;
        const size_t next_t_max = (curr_lower + curr_upper)/2;
        if (next_t_max == curr_t_max) break;
        curr_t_max = next_t_max;
        if (!overflow) it += 1;   // only count refinement steps that produced a layout
    }

    // Nothing usable yet: try the t_max that cannot overflow by construction, then fall back.
    // Returning binning_core's overflow result is never an option - its ranges are {0,0,0}, which
    // makes every single technical bin a lower level IBF, and its bins hold only the sequences that
    // were placed before it ran out of room.
    if (!valid) consider(run(feasible_t_max));
    if (!valid) {
        if (all_seqs) {
            std::vector<size_t> all_ids(fracmin_sigs.size());
            std::iota(all_ids.begin(), all_ids.end(), 0);
            best = fallback_layout(all_ids, sub_bins);
        } else {
            best = fallback_layout(seqs, sub_bins);
        }
    }

    return best;
};

    auto record_bins = [&](const IBF& result, size_t split_start, size_t merge_start, size_t index) {
        for (size_t b = split_start; b < merge_start;) {
            if (result[b].empty()) { b += 1; continue; }
            size_t seq = result[b][0];
            size_t start = b;
            size_t count = 0;
            while (b < merge_start && !result[b].empty() && result[b][0] == seq) { count += 1; b += 1; }
            seq_layout[seq].push_back({index, start, count});
        }
        for (size_t b = merge_start; b < result.size(); b++) {
            for (size_t seq : result[b]) seq_layout[seq].push_back({index, b, 1});
        }
    };

    auto get_upper_lower = [&](const std::vector<size_t>& seqs, size_t sub_bins) {
        std::vector<const std::vector<std::uint64_t>*> ptrs;
        ptrs.reserve(seqs.size());
        for (size_t seq : seqs) ptrs.push_back(&fracmin_sigs[seq]);
        size_t union_size = get_union_size_ptr(ptrs);
        size_t sum = 0;
        for (size_t seq : seqs) sum += fracmin_sigs[seq].size();
        return std::pair{union_size, sum};
    };

    auto get_sub_bins = [&](size_t N) {
        if (N == 0) return size_t{0};
        double raw = std::sqrt(static_cast<double>(N));
        size_t rounded = static_cast<size_t>(std::ceil(raw / 64.0)) * 64;
        if (rounded == 0) rounded = 64;
        return std::min(rounded, static_cast<size_t>(2000));
    };

    const size_t global_bins = get_sub_bins(fracmin_sigs.size());
    size_t union_size = get_union_size(fracmin_sigs);
    size_t sum_size = 0;
    for (const std::vector<std::uint64_t>& sketch : fracmin_sigs) sum_size += sketch.size();

    std::vector<std::vector<IBF>> hibf_levels;
    std::vector<std::vector<std::tuple<size_t,size_t,size_t>>> ranges;
    
    auto free_level_contents = [](std::vector<IBF>& lvl_ibfs) {
        for (IBF& ibf : lvl_ibfs) {
            for (std::vector<size_t>& bin : ibf) {
                std::vector<size_t>{}.swap(bin);
            }
        }
    };

    std::vector<size_t> dummy = {0};
    auto root = refine_and_bin(dummy, global_bins, union_size, sum_size, true);
    hibf_levels.push_back({std::get<0>(root)});
    auto [split_start, split_bins, merge_start] = std::get<1>(root);
    ranges.push_back({std::get<1>(root)});
    record_bins(std::get<0>(root), split_start, merge_start, 0);
    
    std::vector<size_t> trackfill = std::get<2>(root);
    size_t max_bin_id = trackfill.empty() ? 0 : std::distance(trackfill.begin(), std::max_element(trackfill.begin(), trackfill.end()));
    max_bin_ids.push_back({max_bin_id});
    parents.push_back({{max, max}});

    const size_t max_workers = std::max<size_t>(1, threads);

    for (size_t lvl = 0;; lvl++) {
        if (lvl >= 2) free_level_contents(hibf_levels[lvl - 2]);

        std::vector<std::pair<size_t, size_t>> tasks;
        for (size_t ibf_index = 0; ibf_index < hibf_levels[lvl].size(); ibf_index++) {
            const IBF& ibf = hibf_levels[lvl][ibf_index];
            auto [split_start, split_bins, merge_start] = ranges[lvl][ibf_index];
            for (size_t b = merge_start; b < ibf.size(); b++) {
                if (ibf[b].empty()) continue;
                tasks.push_back({ibf_index, b});
            }
        }

        std::cerr << "[hibf] level " << lvl << ": " << hibf_levels[lvl].size()
                  << " IBF(s), " << tasks.size() << " merge bin(s) to expand\n";

        if (tasks.empty()) break;

        std::vector<IBF> next_lvl;
        std::vector<std::tuple<size_t,size_t,size_t>> next_ranges;
        std::vector<size_t> next_max_bin_ids;
        std::vector<std::pair<size_t,size_t>> next_parents;

        next_lvl.reserve(tasks.size());
        next_ranges.reserve(tasks.size());
        next_max_bin_ids.reserve(tasks.size());
        next_parents.reserve(tasks.size());

        for (size_t chunk_start = 0; chunk_start < tasks.size(); chunk_start += max_workers) {
            size_t chunk_end = std::min(tasks.size(), chunk_start + max_workers);
            std::vector<std::future<ChildResult>> futures;
            futures.reserve(chunk_end - chunk_start);

            for (size_t i = chunk_start; i < chunk_end; i++) {
                size_t ibf_index = tasks[i].first;
                size_t b = tasks[i].second;
                const auto& sub_seqs = hibf_levels[lvl][ibf_index][b];
                auto [sub_lower, sub_higher] = get_upper_lower(sub_seqs, global_bins);

                futures.push_back(
                    std::async(std::launch::async, [&, sub_seqs, sub_lower, sub_higher]() -> ChildResult {
                        return refine_and_bin(sub_seqs, global_bins, sub_lower, sub_higher, false);
                    })
                );
            }

            for (size_t i = chunk_start; i < chunk_end; i++) {
                size_t ibf_index = tasks[i].first;
                size_t b = tasks[i].second;
                ChildResult child = futures[i - chunk_start].get();

                auto [c_split_start, c_split_bins, c_merge_start] = std::get<1>(child);
                const std::vector<size_t>& child_trackfill = std::get<2>(child);
                size_t child_max_bin_id = child_trackfill.empty() ? 0 : std::distance(child_trackfill.begin(), std::max_element(child_trackfill.begin(), child_trackfill.end()));

                next_lvl.push_back(std::move(std::get<0>(child)));
                next_max_bin_ids.push_back(child_max_bin_id);
                next_ranges.push_back(std::get<1>(child));
                next_parents.push_back({ibf_index, b});
                record_bins(next_lvl.back(), c_split_start, c_merge_start, next_lvl.size() - 1);
            }
        }

        hibf_levels.push_back(std::move(next_lvl));
        ranges.push_back(std::move(next_ranges));
        max_bin_ids.push_back(std::move(next_max_bin_ids));
        parents.push_back(std::move(next_parents));
    }

    return {hibf_levels, ranges, seq_layout, max_bin_ids, parents};
}