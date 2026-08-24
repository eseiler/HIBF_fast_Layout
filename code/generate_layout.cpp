#include "fast_construct_graph.h"
#include "fast_construct_bin.h"
#include "fast_construct_hashing.h"
#include <lemon/connectivity.h>
#include <fstream>
#include <iostream>
#include <filesystem>

std::vector<std::pair<size_t,size_t>> parse_lvls(const std::string& s){
    std::vector<std::pair<size_t,size_t>> res;
    std::stringstream ss(s);
    std::string tok;
    while(std::getline(ss,tok,',')){
        auto colon = tok.find(':');
        size_t bcount = std::stoul(tok.substr(0,colon));
        size_t bwidth = std::stoul(tok.substr(colon + 1));
        res.push_back({bcount,bwidth});
    }
    return res;
}

int main(int argc, char* argv[]){
    if(argc != 13){
        std::cerr << "Missing parameters: " << "Levels used for LSH + Fraction s for Fracmin. + Refinement count when constructing the HIBF + Max Levels of HIBF\n";
        return 1;
    }

    const std::filesystem::path& dir_path = argv[1];
    const std::uint8_t k = static_cast<std::uint8_t>(std::stoul(argv[2]));
    const std::uint8_t q = static_cast<std::uint8_t>(std::stoul(argv[3]));
    const std::uint32_t w = static_cast<std::uint64_t>(std::stoull(argv[4]));
    const std::uint64_t seed = static_cast<std::uint64_t>(std::stoull(argv[5])) & ((1ULL << 2*q) - 1);
    const double fpr = std::stod(argv[6]);
    std::vector<std::pair<size_t, size_t>> lvls = parse_lvls(argv[7]);
    double s = std::stod(argv[8]);
    size_t refinements = std::stoul(argv[9]);
    const std::uint8_t hash_funcs = static_cast<std::uint8_t>(std::stoul(argv[10]));
    std::ofstream out_path(argv[11]);
    size_t threads = std::stoul(argv[12]);  
    // Generate One Permutation Hashes for each "sequence":
    IntHasher hasher;
    std::tuple<std::vector<std::vector<std::uint64_t>>, std::vector<std::vector<std::uint64_t>>, std::unordered_map<size_t, std::string>> sigs = ophs_fmhs(dir_path, q, k, w, seed, s, hasher, threads);
    auto full_hibf = generate_hibf<standardHasher>(sigs, lvls, s, fpr, hash_funcs, refinements, threads);
    std::unordered_map<size_t, std::string>& seq_to_file = std::get<2>(sigs);

    write_linkage(out_path, seq_to_file);
    write_config(out_path, dir_path, q, w, fpr, hash_funcs, seq_to_file.size());
    write_header(out_path, std::get<0>(full_hibf), std::get<3>(full_hibf), std::get<4>(full_hibf));
    write_content(out_path, std::get<2>(full_hibf));
    out_path.close();

    return 0;
}
