#pragma once

#include <vector>
#include <cstdint>
#include <functional>
#include <filesystem>
#include <utility>
#include <algorithm>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#define XXH_INLINE_ALL
#include "xxhash.h"

// @brief a wrapper around the uint64 hash of xxHash.
// @param input : the integer to hash
std::uint64_t xxhash_wrap(std::uint64_t input);

// @brief Given a set of hashes, compute the one permutation hash and the Fracmin hash of the given collection
// @param hashes : the set of hashes to be hashed
// @param k : the intervall the function is cut in (used for One Permutation hash).
// @param s : The fraction used to compute the Fracmin Hash
// @param hashFunc : the Hash Function to be used
//
// @note in the end, returns a vector of size 2^k (OPH).
// @note this function is rather for debugging purposes, when actually using stuff, refer to @func ophs
// @note the first set of vectors is the OPH, the other one is the Fracmin.
std::tuple<std::vector<std::vector<std::uint64_t>>, std::vector<std::vector<std::uint64_t>>, std::unordered_map<size_t, std::string>> one_permutation_fracmin_hash(const std::vector<std::vector<std::uint64_t>>& hashes,
                                                const std::uint8_t k, const double s, 
                                                std::function<std::uint64_t(std::uint64_t)> hashFunc = xxhash_wrap);


struct IntHasher {
    template <typename T>
    std::uint64_t operator()(const T elem) const {return xxhash_wrap(elem);}
};

// @brief Given a set of sequences, k and a hashfunction, compute the one permutation hash of the given sequences.
// @param filepath : Path to the file containing the sequences to be hashed
// @param q : the size of the Q-Grams applied to the given sequences
// @param k : the Intervall the function is cut in (Used in One Permutation Hash).
// @param w : The window size used for every sequence.
// @param seed : The seed used for computing the windows.
// @param s: The fraction used to compute the Fracmin Hash
// @param threads : Number of Threads used.
//
// @note every vector inside of the returned one has a size of 2^k (OPH).
// @note the first set of vectors is the OPH, the other one is the Fracmin.
template <typename IntHash>
std::tuple<std::vector<std::vector<std::uint64_t>>, std::vector<std::vector<std::uint64_t>>, std::unordered_map<size_t, std::string>> ophs_fmhs(const std::filesystem::path& filepath,
                                              const std::uint8_t q, const std::uint8_t k, const std::uint32_t w, const std::uint64_t seed, const double s, IntHash&& hashFunc, const size_t threads);
#include "templates/ophs_fmhs.tpp"

// @brief What ophs_fmhs returns: the One Permutation Hashes, the Fracmin sketches, and the mapping
//        from sequence id to the file it was read from.
using SketchSet = std::tuple<std::vector<std::vector<std::uint64_t>>,
                             std::vector<std::vector<std::uint64_t>>,
                             std::unordered_map<size_t, std::string>>;

// @brief The parameters the sketches depend on. A cache written with different ones is not reusable.
struct SketchParams {
    std::filesystem::path dir_path;
    std::uint8_t q;
    std::uint8_t k;
    std::uint32_t w;
    std::uint64_t seed;
    double s;
};

// @brief Read sketches back from a cache file written by @func write_sketch_cache.
// @param path : the cache file.
// @param params : the parameters the caller is about to compute with.
// @param out : filled with the cached sketches on success, untouched otherwise.
// @return true when the file exists, is intact, and was written with exactly these parameters.
//
// @note hashing a directory of FASTA files is by far the slowest part of a layout run, so during
//       testing the same sketches can be reused across many layout runs.
bool read_sketch_cache(const std::filesystem::path& path, const SketchParams& params, SketchSet& out);

// @brief Write sketches to a cache file, together with the parameters they were computed with.
// @param path : the cache file, overwritten if it exists.
// @param params : the parameters the sketches were computed with.
// @param sketches : the sketches to store.
// @return true when the file was written completely.
bool write_sketch_cache(const std::filesystem::path& path, const SketchParams& params, const SketchSet& sketches);

// @brief Given a set of (sorted) Fracmin Sketches, compute the size of the union of those sketches,
// @param sketches : The set of Fracmin Sketches.
size_t get_union_size(const std::vector<std::vector<std::uint64_t>>& sketches);

// @brief Given a set of pointers to sketches, compute their Fracmin Union size.
// @param sketches : The set of pointers pointing to corresponding Fracmin Sketches
// @note Computation is nearly exact as @func get_union_size.
size_t get_union_size_ptr(const std::vector<const std::vector<std::uint64_t>*>& sketches);

// @brief Given a set of One-Permutation Hashes, extract the wanted number of bands, given the bandwidth
// @param hashes : The set of One-Permutation Hashes.
// @param bwidth : The size of one band.
// @param bcount : The number of bands to be extracted.
//
// @note for any i : hashes[i].size() > bwidth * bcount.
// @note EveryHash(CollectionOfHashes(SingularBand)) is the nesting
std::vector<std::vector<std::vector<std::uint64_t>>> extract_bands(const std::vector<std::vector<std::uint64_t>>& hashes,
                                                                   const size_t bwidth, const size_t bcount);

// @brief Hash the given Vector. Using XXHash (https://github.com/cyan4973/xxhash)
// @param vec the vector to be hashed
template <typename T>
std::uint64_t hash_vector(const std::vector<T>& vec){
    return XXH64(vec.data(), vec.size() * sizeof(T), 0);
}

// @brief Hash the given vector by returning its size.
// @param vec the vector to be hashed.
// @note this function is debugging purposes and for little examples only!!!
template <typename T>
std::uint64_t hash_vector_debug(const std::vector<T>& vec){
  std::uint64_t res = 0;
  for(std::uint64_t elem: vec){
    res += elem;
  }
  return res;
}