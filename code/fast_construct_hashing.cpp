#include "fast_construct_hashing.h"
#include <limits>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>


std::uint64_t xxhash_wrap(std::uint64_t input){ return XXH64(&input, sizeof(input), 0); }

std::tuple<std::vector<std::vector<std::uint64_t>>, std::vector<std::vector<std::uint64_t>>, std::unordered_map<size_t, std::string>> one_permutation_fracmin_hash(const std::vector<std::vector<std::uint64_t>>& hashes, const std::uint8_t k, const double s, std::function<std::uint64_t(std::uint64_t)> hashFunc){
  std::vector<std::vector<std::uint64_t>> res_oph;
  res_oph.reserve(hashes.size());
  std::vector<std::vector<std::uint64_t>> res_fracmin;
  res_fracmin.reserve(hashes.size());
  std::unordered_map<size_t, std::string> seq_to_path;
  size_t seq_id = 0;

  for(const std::vector<std::uint64_t>& hashes_ : hashes){
    std::uint64_t max = std::numeric_limits<std::uint64_t>::max();
    std::vector<std::uint64_t> tmp_oph(std::pow(2,k), max);
    std::vector<std::uint64_t> tmp_fracmin;
    std::uint64_t thresh = static_cast<std::uint64_t>(static_cast<double>(max*s));
    for(std::uint64_t hash : hashes_){
      std::uint64_t hash_ = hashFunc(hash);

      // Fracmin Part
      if(hash_ < thresh) tmp_fracmin.push_back(hash_);

      // OPH Part
      std::uint8_t index = hash_ & ((1ULL << k) -1); //Extract the first k bits of the initial hash value.
      std::uint64_t val = hash_ >> k; // remove the k bits, they'd just be a bias since evry value in the bucket would ALWAYS have these bits same.
      if(val < tmp_oph[index]) tmp_oph[index] = val;
    }
    std::sort(tmp_fracmin.begin(), tmp_fracmin.end());

    std::string test_file = "/path/to/file/" + std::to_string(seq_id) + ".fasta";
    seq_to_path[seq_id] = test_file;
    seq_id += 1;
    res_oph.push_back(tmp_oph);
    res_fracmin.push_back(tmp_fracmin);
  }
  return {std::move(res_oph), std::move(res_fracmin), std::move(seq_to_path)};
}


namespace {

constexpr char sketch_cache_magic[8] = {'F','L','S','K','C','H','0','1'};

template <typename T>
void put(std::ostream& out, const T& value){ out.write(reinterpret_cast<const char*>(&value), sizeof(T)); }

template <typename T>
bool get(std::istream& in, T& value){ return static_cast<bool>(in.read(reinterpret_cast<char*>(&value), sizeof(T))); }

void put_sketches(std::ostream& out, const std::vector<std::vector<std::uint64_t>>& sketches){
    put(out, static_cast<std::uint64_t>(sketches.size()));
    for(const std::vector<std::uint64_t>& sketch : sketches){
        put(out, static_cast<std::uint64_t>(sketch.size()));
        if(!sketch.empty()) out.write(reinterpret_cast<const char*>(sketch.data()), sketch.size() * sizeof(std::uint64_t));
    }
}

bool get_sketches(std::istream& in, std::vector<std::vector<std::uint64_t>>& sketches, std::uint64_t expected){
    std::uint64_t count = 0;
    if(!get(in, count) || count != expected) return false;
    sketches.assign(count, {});
    for(std::vector<std::uint64_t>& sketch : sketches){
        std::uint64_t len = 0;
        if(!get(in, len)) return false;
        sketch.resize(len);
        if(len && !in.read(reinterpret_cast<char*>(sketch.data()), len * sizeof(std::uint64_t))) return false;
    }
    return true;
}

void put_string(std::ostream& out, const std::string& str){
    put(out, static_cast<std::uint64_t>(str.size()));
    if(!str.empty()) out.write(str.data(), str.size());
}

bool get_string(std::istream& in, std::string& str){
    std::uint64_t len = 0;
    if(!get(in, len)) return false;
    str.resize(len);
    if(len && !in.read(str.data(), len)) return false;
    return true;
}

} // namespace

bool write_sketch_cache(const std::filesystem::path& path, const SketchParams& params, const SketchSet& sketches){
    std::ofstream out(path, std::ios::binary);
    if(!out) return false;

    out.write(sketch_cache_magic, sizeof(sketch_cache_magic));
    put(out, params.q);
    put(out, params.k);
    put(out, params.w);
    put(out, params.seed);
    put(out, params.s);
    put_string(out, params.dir_path.string());

    const std::vector<std::vector<std::uint64_t>>& oph = std::get<0>(sketches);
    const std::vector<std::vector<std::uint64_t>>& fmh = std::get<1>(sketches);
    const std::unordered_map<size_t, std::string>& paths = std::get<2>(sketches);

    put(out, static_cast<std::uint64_t>(oph.size()));
    put_sketches(out, oph);
    put_sketches(out, fmh);

    put(out, static_cast<std::uint64_t>(paths.size()));
    for(const auto& [seq_id, file] : paths){
        put(out, static_cast<std::uint64_t>(seq_id));
        put_string(out, file);
    }

    out.flush();
    return static_cast<bool>(out);
}

bool read_sketch_cache(const std::filesystem::path& path, const SketchParams& params, SketchSet& out_sketches){
    if(path.empty()) return false;

    std::ifstream in(path, std::ios::binary);
    if(!in) return false;

    char magic[sizeof(sketch_cache_magic)] = {};
    if(!in.read(magic, sizeof(magic)) || std::memcmp(magic, sketch_cache_magic, sizeof(magic)) != 0){
        std::cerr << "[sketch cache] " << path.string() << " is not a sketch cache, recomputing\n";
        return false;
    }

    std::uint8_t q = 0, k = 0;
    std::uint32_t w = 0;
    std::uint64_t seed = 0;
    double s = 0.0;
    std::string dir;
    if(!get(in, q) || !get(in, k) || !get(in, w) || !get(in, seed) || !get(in, s) || !get_string(in, dir)) return false;

    if(q != params.q || k != params.k || w != params.w || seed != params.seed || s != params.s || dir != params.dir_path.string()){
        std::cerr << "[sketch cache] " << path.string() << " was written with different parameters, recomputing\n";
        return false;
    }

    std::uint64_t count = 0;
    if(!get(in, count)) return false;

    SketchSet loaded;
    if(!get_sketches(in, std::get<0>(loaded), count)) return false;
    if(!get_sketches(in, std::get<1>(loaded), count)) return false;

    std::uint64_t path_count = 0;
    if(!get(in, path_count)) return false;
    std::get<2>(loaded).reserve(path_count);
    for(std::uint64_t i = 0; i < path_count; i++){
        std::uint64_t seq_id = 0;
        std::string file;
        if(!get(in, seq_id) || !get_string(in, file)) return false;
        std::get<2>(loaded).emplace(static_cast<size_t>(seq_id), std::move(file));
    }

    out_sketches = std::move(loaded);
    return true;
}

size_t get_union_size(const std::vector<std::vector<std::uint64_t>>& sketches){
  std::unordered_set<std::uint64_t> elems;
  size_t tots = 0; // Used for reserving space
  for(const std::vector<std::uint64_t>& sketch : sketches) tots += sketch.size();
  elems.reserve(tots); // "Worst Case" every Element is unique.

  for(const std::vector<std::uint64_t>& sketch : sketches) for(std::uint64_t elem : sketch) elems.insert(elem);

  return elems.size();
}

size_t get_union_size_ptr(const std::vector<const std::vector<std::uint64_t>*>& sketches){
  std::unordered_set<std::uint64_t> elems;
  size_t tots = 0; // Used for reserving space
  for(const std::vector<std::uint64_t>* sketch : sketches) tots += sketch->size();
  elems.reserve(tots); // "Worst Case" every Element is unique.

  for(const std::vector<std::uint64_t>* sketch : sketches) for(std::uint64_t elem : *sketch) elems.insert(elem);

  return elems.size();
}

// @note this function is kind of redacted, sind its entire functionality is now in connect_bands. With the advantage, that we dont have to store every band respective. 
//  It is still used in the tests as the loop stays the same and is used in the connect_bands function.
std::vector<std::vector<std::vector<std::uint64_t>>> extract_bands(const std::vector<std::vector<std::uint64_t>>& hashes, const size_t bwidth, const size_t bcount){
  std::vector<std::vector<std::vector<std::uint64_t>>> res;
  for(std::vector<std::uint64_t> oph : hashes){ // itterate over every given hash
    std::vector<std::vector<std::uint64_t>> tmp;
    for(size_t count = 0; count < bcount; count++){ // We extract the Bands from the given signature
      auto start = oph.begin() + count * bwidth;
      auto end = oph.begin() + (count + 1) * bwidth;
      tmp.push_back(std::vector<std::uint64_t>(start,end)); // Work with itterators - efficiency
    }
    res.push_back(tmp);
  }
  return res;
}
