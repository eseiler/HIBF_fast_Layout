#include "../fast_construct_hashing.h"
#include <thread>
#include <atomic>
#include <limits>
#include <cmath>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>

template <typename IntHash>
std::tuple<std::vector<std::vector<std::uint64_t>>, std::vector<std::vector<std::uint64_t>>, std::unordered_map<size_t, std::string>> ophs_fmhs(const std::filesystem::path& dirpath, const std::uint8_t q, const std::uint8_t k, const std::uint32_t w, const std::uint64_t seed, const double s, IntHash&& hashFunc, const size_t threads){

    std::vector<std::filesystem::path> files;
    for(auto const& entry: std::filesystem::directory_iterator{dirpath}){
        if(!entry.is_regular_file()) continue;

        auto const& filepath = entry.path();
        auto ext = filepath.extension().string();

        if(ext != ".fa" && ext != ".fasta" && ext !=".gz" && ext != ".fq" && ext !=".fastq") continue;

        files.push_back(filepath);
    }


    const size_t n_threads = (threads > files.size()) ? files.size() : threads;

    std::vector<std::vector<std::uint64_t>> res_oph(files.size());
    std::vector<std::vector<std::uint64_t>> res_fmh(files.size());
    std::vector<std::string> seq_paths(files.size());

    std::atomic<std::size_t> next_file{0};

    auto worker = [&](){
        seqan3::shape shape = seqan3::ungapped{q};
        auto minimiser_view = seqan3::views::minimiser_hash(shape, seqan3::window_size{w}, seqan3::seed{seed});
        std::uint64_t max = std::numeric_limits<std::uint64_t>::max();
        std::uint64_t thresh = static_cast<std::uint64_t>(static_cast<double>(max*s));

        while(true) {
            const size_t seq_id = next_file.fetch_add(1, std::memory_order_relaxed);
            if(seq_id >= files.size()) break;
            const auto& filepath = files[seq_id];
            auto fin = seqan3::sequence_file_input{filepath};
            auto record = *fin.begin();
            std::vector<std::uint64_t> tmp_oph(std::pow(2,k), max);
            std::vector<std::uint64_t> tmp_fracmin;
            for(auto&& qgram : record.sequence() | minimiser_view){
              std::uint64_t hash = hashFunc(qgram);

              if(hash < thresh) tmp_fracmin.push_back(hash);

              std::uint8_t index = hash & ((1ULL << k) -1);
              std::uint64_t val = hash >> k; // remove the k bits, they'd just be a bias since every value in the bucket would ALWAYS have these bits same.

              if(val < tmp_oph[index]) tmp_oph[index] = val;
            }
          std::sort(tmp_fracmin.begin(), tmp_fracmin.end());
          seq_paths[seq_id] = filepath.string();
          res_fmh[seq_id] = (std::move(tmp_fracmin));
          res_oph[seq_id] = (std::move(tmp_oph));
        }
    };

    std::vector<std::thread> workers;
    workers.reserve(n_threads);

    for(size_t i = 0; i < n_threads; i++){
        workers.emplace_back(worker);
    }
    for(auto& worker : workers) worker.join();


    std::unordered_map<size_t, std::string> seq_to_path;
    seq_to_path.reserve(files.size());

    for(size_t seq_id = 0; seq_id < files.size(); seq_id++){
        seq_to_path.emplace(seq_id, std::move(seq_paths[seq_id]));
    }

 return {std::move(res_oph), std::move(res_fmh), std::move(seq_to_path)};
  }
