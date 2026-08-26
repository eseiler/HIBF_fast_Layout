#include "fast_construct_bin.h"
#include "fast_construct_hashing.h"

std::vector<double> compute_fcorrs(double fpr, size_t h){
    std::vector<double> res(101, 1.0);
    const double log_num = std::log(1.0 - std::pow(fpr, 1.0 / static_cast<double>(h)));

    for(size_t s = 1; s <= 100; s++){
        double pcorr = 1.0 - std::pow(1.0 - fpr, 1.0 / static_cast<double>(s));
        double log_denom = std::log(1.0 - std::pow(pcorr, 1.0 / static_cast<double>(h)));
        res[s] = log_num / log_denom;
    }
    return res;
}

size_t merge_average(const std::vector<size_t>& track_fill, const size_t merge_start){
    size_t sum = 0;
    size_t amount = 0;
    for(size_t b = merge_start; b < track_fill.size(); b++){
        if(track_fill[b] == 0) continue;
        sum += track_fill[b];
        amount += 1;
    }
    return amount ? sum / amount : 0;
}

size_t splitting_average(const std::vector<size_t>& track_fill, const size_t split_start, const size_t split_end){
    size_t sum = 0;
    size_t amount = 0;
    for(size_t b = split_start; b < split_end; b++){
        if(track_fill[b] == 0) continue;
        sum += track_fill[b];
        amount += 1;
    }
    return amount ? sum / amount : 0;
}

void write_linkage(std::ostream& out, const std::unordered_map<size_t, std::string>& seq_to_path){
    out << "@CHOPPER_USER_BINS\n";
    for(size_t i = 0; i < seq_to_path.size(); ++i){
        out << "@" << i << " " << seq_to_path.at(i) << '\n';
    }
    out << "@CHOPPER_USER_BINS_END\n";
}

void write_config(std::ofstream& out,
                  std::filesystem::path const& dir_path,
                  std::uint8_t q,
                  std::uint32_t w,
                  double fpr,
                  std::uint8_t hash_funcs,
                  size_t user_bins){

out << "@CHOPPER_CONFIG\n"
        "@{\n"
        "@    \"chopper_config\": {\n"
        "@        \"version\": 2,\n"
        "@        \"data_file\": {\n"
    <<  "@            \"value0\": \"" << dir_path.string() << "\"\n"
        "@        },\n"
        "@        \"debug\": false, \n"
        "@        \"sketch_directory\": {\n"
        "@            \"value0\": \"\"\n"
        "@        },\n"
    <<  "@        \"k\": " << static_cast<unsigned>(q) << ",\n"
    <<  "@        \"window_size\": " << w << ",\n"
        "@        \"disable_sketch_output\": true,\n"
        "@        \"precomputed_files\": false,\n"
        "@        \"maximum_index_size\": 0,\n"
        "@        \"number_of_partitions\": 1,\n"
        "@        \"output_filename\": {\n"
    <<  "@            \"value0\": \"" << dir_path.string() << "\"\n"
        "@        },\n"
        "@        \"determine_best_tmax\": false,\n"
        "@        \"force_all_binnings\": false\n"
        "@    }\n"
        "@}\n"
        "@CHOPPER_CONFIG_END\n"
        "@HIBF_CONFIG\n"
        "@{\n"
        "@    \"hibf_config\": {\n"
        "@        \"version\": 1,\n"
    <<  "@        \"number_of_user_bins\":" <<  user_bins << ",\n"
    <<  "@        \"number_of_hash_functions\": " << static_cast<unsigned>(hash_funcs) << ",\n"
    <<  "@        \"maximum_fpr\": " <<fpr << ",\n"
    <<  "@        \"relaxed_fpr\": " << fpr << ",\n"
        "@        \"threads\": 1,\n"
        "@        \"sketch_bits\": 12,\n"
        "@        \"tmax\": 0,\n"
        "@        \"alpha\": 1.2,\n"
        "@        \"max_rearrangement_ratio\": 0.5,\n"
        "@        \"disable_estimate_union\": false,\n"
        "@        \"disable_rearrangement\": false\n"
        "@    }\n"
        "@}\n"
        "@HIBF_CONFIG_END\n";
}

void write_header(std::ostream& out, const std::vector<std::vector<IBF>>& hibf_levels, const std::vector<std::vector<size_t>>& max_bin_ids, const std::vector<std::vector<std::pair<size_t,size_t>>>& parents){
    auto merge_bin_label = [](size_t level, size_t index, const std::vector<std::vector<std::pair<size_t,size_t>>>& parents){
        std::vector<size_t> chain;
        while(level > 0){
            auto [par_idx, merge_bin] = parents[level][index];
            chain.push_back(merge_bin);
            index = par_idx;
            level -= 1;
        }
        std::reverse(chain.begin(), chain.end());
        std::string lab;
        for(size_t i = 0; i < chain.size(); i++){
            if(i > 0) lab += ";";
            lab += std::to_string(chain[i]);
        }
        return lab;
    };

    for(size_t level = 0; level < hibf_levels.size(); level++){
        for(size_t ibf_idx = 0; ibf_idx < hibf_levels[level].size(); ibf_idx++){
            size_t const max_bin_id = max_bin_ids[level][ibf_idx];

            if(level == 0) out << "#TOP_LEVEL_IBF fullest_technical_bin_idx:" << max_bin_id << "\n";
            else{
                std::string const label = merge_bin_label(level, ibf_idx, parents);
                out << "#LOWER_LEVEL_IBF_" << label << " fullest_technical_bin_idx:" << max_bin_id << "\n";
            }
        }
    }

    out << "#USER_BIN_IDX\tTECHNICAL_BIN_INDICES\tNUMBER_OF_TECHNICAL_BINS\n";
}

void write_content(std::ofstream& out, const std::unordered_map<size_t, std::vector<std::tuple<size_t,size_t,size_t>>>& seq_layout){
    for(size_t seq = 0; seq < seq_layout.size(); ++seq){
        auto const & entries = seq_layout.at(seq);
        std::string bin_indices;
        std::string number_of_bins;

        for(size_t i = 0; i < entries.size(); i++){
            auto const& [index, start, count] = entries[i];
            if(i > 0) {bin_indices += ";"; number_of_bins += ";";}
            bin_indices += std::to_string(start);
            number_of_bins += std::to_string(count);
        }
        out << seq << "\t" << bin_indices << "\t" << number_of_bins << "\n";
    }
}
