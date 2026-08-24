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
    
auto refine_and_bin = [&](const std::vector<size_t>& seqs, size_t sub_bins, size_t lower, size_t upper, bool all_seqs) {
    bool valid = false;
    std::tuple<std::vector<std::vector<size_t>>, std::tuple<size_t,size_t,size_t>, std::vector<size_t>, bool> b_res;
    std::tuple<std::vector<std::vector<size_t>>, std::tuple<size_t,size_t,size_t>, std::vector<size_t>, bool> res;
    size_t curr_lower = lower;
    size_t curr_upper = upper;
    size_t curr_t_max = (curr_lower + curr_upper)/(2*sub_bins*s);
    size_t old_t_max = 0;
    size_t best_t_max = curr_t_max;

    curr_upper = curr_t_max * 2;
    curr_lower = 0;

    auto used_bin_count = [](const auto& res_tuple) {
        const auto& ibf = std::get<0>(res_tuple);
        size_t count = 0;
        for (const auto& bin : ibf) if (!bin.empty()) count += 1;
        return count;
    };

    if (all_seqs) {
        std::cerr << "[refine_and_bin] ROOT START: sub_bins = " << sub_bins
                  << ", lower = " << lower << ", upper = " << upper
                  << ", initial t_max = " << curr_t_max << "\n";
    }

    for(size_t it = 0; it <= p;){
        if(curr_t_max == old_t_max) {
            if (all_seqs) {
                std::cerr << "[refine_and_bin] ROOT Convergence reached at it = " << it << " (t_max = " << curr_t_max << ")\n";
            }
            break; 
        }

        res = all_seqs ? binning(labMaps, clusts, fracmin_sigs, s, sub_bins, curr_t_max, fcorrs) 
                       : binning_given_seqs(labMaps, clusts, fracmin_sigs, seqs, s, sub_bins, curr_t_max, fcorrs);

        bool overflow = std::get<3>(res);
        if(overflow){
            if (all_seqs) {
                std::cerr << "[refine_and_bin] ROOT OVERFLOW at t_max = " << curr_t_max << " (iteration not counted)\n";
            }
            curr_lower = curr_t_max;
            if(curr_t_max == 0) break;

            size_t next_t_max = (curr_lower + curr_upper) / 2;
            if (next_t_max == curr_t_max) break; // Intervall lässt sich nicht weiter verkleinern

            curr_t_max = next_t_max;
            continue; 
        }
        if(it == p) break; 

        const auto& [split_start, split_bins_cnt, merge_start] = std::get<1>(res);
        const std::vector<std::vector<size_t>>& result = std::get<0>(res);

        valid = true;
        b_res = res;
        best_t_max = curr_t_max; // Speichere das t_max des korrekten Ergebnisses

        size_t split_bin_amt = merge_start - split_start;
        size_t merge_bin_amt = result.size() - merge_start;

        if(split_bin_amt == 0 && merge_bin_amt == 0) {
            if (all_seqs) {
                std::cerr << "[refine_and_bin] ROOT Empty IBF result, stopping refinement.\n";
            }
            break; 
        }

        const std::vector<size_t>& trackfill = std::get<2>(res);
        size_t split_avg = split_bin_amt ? splitting_average(trackfill, split_start, merge_start) : 0;
        size_t merge_avg = merge_bin_amt ? merge_average(trackfill, merge_start) : 0;

        if (all_seqs) {
            std::cerr << "[refine_and_bin] ROOT [it " << it << "] t_max = " << curr_t_max 
                      << ", used_bins = " << used_bin_count(res)
                      << ", split_avg = " << split_avg 
                      << ", merge_avg = " << merge_avg << "\n";
        }

        if (split_bin_amt == 0 || split_avg < merge_avg) {
            curr_upper = curr_t_max;
        } else {
            curr_lower = curr_t_max;
        }

        if(curr_t_max == 0) break; 
        old_t_max = curr_t_max;    
        curr_t_max = (curr_lower + curr_upper)/2;
        it += 1;                   
    }

    bool final_is_overflow = std::get<3>(res);
    auto final_res = (final_is_overflow && valid) ? b_res : res;
    size_t chosen_t_max = (final_is_overflow && valid) ? best_t_max : curr_t_max;

    if (all_seqs) {
        std::cerr << "[refine_and_bin] ROOT END: chosen t_max = " << chosen_t_max
                  << ", used_bins = " << used_bin_count(final_res)
                  << ", overflow_flag = " << std::get<3>(final_res) << "\n";
    }

    return final_res;
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