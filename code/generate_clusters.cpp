#include "fast_construct_graph.h"
#include "test_clustering.h"
#include "fast_construct_bin.h"
#include "fast_construct_hashing.h"
#include <lemon/connectivity.h>
#include <fstream>
#include <iostream>
#include <filesystem>

// Helper Function to parse the set of doubles
std::vector<double> parse_doubles(const std::string& s){
    std::vector<double> res;
    std::stringstream ss(s);
    std::string tok;
    while(std::getline(ss,tok,',')) res.push_back(std::stod(tok));
    return res;
}

// Helper Function to parse the set of size_ts
std::vector<size_t> parse_sizes(const std::string& s){
    std::vector<size_t> res;
    std::stringstream ss(s);
    std::string tok;
    while(std::getline(ss,tok,',')) res.push_back(std::stoul(tok));
    return res;
}

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

std::string printout(const size_t seq){
    return "Sequenz " + std::to_string(seq);
}

std::string printout(const std::vector<size_t>& labels){
    std::string res = "Sequenzen ";
    for(size_t lab : labels) res += (std::to_string(lab) += ", ");
    return res;
}

void printgraph(lemon::ListGraph& graph, std::unordered_map<size_t, lemon::ListGraph::Node>& labMap, std::string filepath){
    // Get the connected Components in order to colorize them later on.
    lemon::ListGraph::NodeMap<int> compMap(graph);
    int numComp = lemon::connectedComponents(graph, compMap);

    // In Order to print, we invert the labMap
    std::map<lemon::ListGraph::Node, size_t> nodeMap;
    for(auto& [label, node] : labMap) nodeMap[node] = label;

    // Potential colors, cycle around if there are more than 10 con.comps.
    std::vector<std::string> colors = {
        "red", "blue", "green", "orange", "purple", "cyan", "magenta", "yellow", "brown", "pink"
    };

    std::ofstream dot(filepath);
    dot << "graph G {\n" << "  node [style = filled];\n";
    for(lemon::ListGraph::NodeIt n(graph); n!= lemon::INVALID; ++n){ // Print every node and its label
        std::string color = colors[compMap[n] % colors.size()]; // Get the color, wrap around if more than colors contained.
        dot << "  " << graph.id(n) << " [label=\"" << printout(nodeMap[n]) << "\", color=\"" << color << "\"];\n";
    }
    for(lemon::ListGraph::EdgeIt e(graph); e!= lemon::INVALID; ++e){ // Print every edge and its nodes
        dot << "  " << graph.id(graph.u(e)) << " -- " << graph.id(graph.v(e)) << ";\n";
    }
    dot << "}\n";
    dot.close();
}

template <typename Hasher>
void printgraph(lemon::ListGraph& graph, std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, Hasher>& labMap, std::string filepath){
    // Get the connected Components in order to colorize them later on.
    lemon::ListGraph::NodeMap<int> compMap(graph);
    int numComp = lemon::connectedComponents(graph, compMap);

    // In Order to print, we invert the labMap
    std::map<lemon::ListGraph::Node, std::vector<size_t>> nodeMap;
    for(auto& [label, node] : labMap) nodeMap[node] = label;

    // Potential colors, cycle around if there are more than 10 con.comps.
    std::vector<std::string> colors = {
        "red", "blue", "green", "orange", "purple", "cyan", "magenta", "yellow", "brown", "pink"
    };

    std::ofstream dot(filepath);
    dot << "graph G {\n" << "  node [style = filled];\n";
    for(lemon::ListGraph::NodeIt n(graph); n!= lemon::INVALID; ++n){ // Print every node and its label
        std::string color = colors[compMap[n] % colors.size()]; // Get the color, wrap around if more than colors contained.
        dot << "  " << graph.id(n) << " [label=\"" << printout(nodeMap[n]) << "\", color=\"" << color << "\"];\n";
    }
    for(lemon::ListGraph::EdgeIt e(graph); e!= lemon::INVALID; ++e){ // Print every edge and its nodes
        dot << "  " << graph.id(graph.u(e)) << " -- " << graph.id(graph.v(e)) << ";\n";
    }
    dot << "}\n";
    dot.close();
}

template <typename Hasher>
void printgraph(lemon::ListGraph& graph, std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, Hasher>>& labMaps, std::string filepath){
    // Get the connected Components in order to colorize them later on.
    lemon::ListGraph::NodeMap<int> compMap(graph);
    int numComp = lemon::connectedComponents(graph, compMap);

    // In Order to print, we invert the labMap
    std::map<lemon::ListGraph::Node, std::vector<size_t>> nodeMap;

    size_t cluster_id = 0;
    std::map<lemon::ListGraph::Node, size_t> node_to_id; // Used in Legend later on
    std::vector<lemon::ListGraph::Node> nodes_ordered;
    
    for(auto& labMap : labMaps){ // labMap : unordered_map<vector, Node, Hasher>.
        for(auto& [label, node] : labMap){
            nodeMap[node] = label;

            if(node_to_id.find(node) != node_to_id.end()) continue;

            node_to_id[node] = cluster_id;
            nodes_ordered.push_back(node);
            cluster_id += 1;
        } 
    }
    // Potential colors, cycle around if there are more than 10 con.comps.
    std::vector<std::string> colors = {
        "red", "blue", "green", "orange", "purple", "cyan", "magenta", "yellow", "brown", "pink"
    };

    std::ofstream dot(filepath);
    dot << "graph G {\n" << "  node [style = filled];\n";
    for(lemon::ListGraph::NodeIt n(graph); n!= lemon::INVALID; ++n){ // Print every node and its label
        std::string color = colors[compMap[n] % colors.size()]; // Get the color, wrap around if more than colors contained.
        const std::vector<size_t>& labs = nodeMap[n];
        dot << "  " << graph.id(n) << " [label=\"S" << node_to_id[n] << " (" << labs.size() << ")\"" << ", color=\"" << color << "\"];\n";
    }
    for(lemon::ListGraph::EdgeIt e(graph); e!= lemon::INVALID; ++e){ // Print every edge and its nodes
        dot << "  " << graph.id(graph.u(e)) << " -- " << graph.id(graph.v(e)) << ";\n";
    }

    // Print the legend as well
    dot << " legend [shape = box, style = filled, fillcolor = lightyellow, label=\"";
    for(auto n: nodes_ordered){
        const std::vector<size_t>& labs = nodeMap[n];
        dot << "S" << node_to_id[n] << " (" << labs.size() << "): ";
        for(size_t seq = 0; seq < labs.size(); seq++){
            dot << labs[seq];
            if(seq + 1 < labs.size()) dot << ", ";
        }
        dot << "\\l";
    }
    dot << "\"];\n";

    if(!nodes_ordered.empty()) dot << " " << graph.id(nodes_ordered.back()) << " -- legend [style = invis, minlen = 7]";

    dot << "}\n";
    dot.close();
}

int main(int argc, char* argv[]){
    if(argc != 7){
        std::cerr << "Missing parameters: " << "Amount cluster + cluster size + vector_size + Similarity in cluster + Similarity berween clusters + Fraction s for Fracmin. + Refinement count when constructing the HIBF\n";
        return 1;
    }

    
    size_t vec_size = std::stoul(argv[1]);
    std::vector<double> j_sims = parse_doubles(argv[2]);
    std::vector<size_t> counts = parse_sizes(argv[3]);
    std::vector<std::pair<size_t, size_t>> lvls = parse_lvls(argv[4]);
    double s = std::stod(argv[5]);
    size_t refinements = std::stoul(argv[6]);
    
    std::vector<std::vector<std::uint64_t>> rand_clusts = get_any_cluster(vec_size, j_sims, counts, 0);
    
    lemon::ListGraph graph;
    
    // Generate One Permutation Hashes for each "sequence":
    std::tuple<std::vector<std::vector<std::uint64_t>>, std::vector<std::vector<std::uint64_t>>, std::unordered_map<size_t, std::string>> sigs = one_permutation_fracmin_hash(rand_clusts, 8, s);
    std::vector<std::vector<std::uint64_t>> oph_sigs = std::get<0>(sigs);
    std::vector<std::vector<std::uint64_t>> fracmin_sigs = std::get<1>(sigs);
    
    std::vector<std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, standardHasher>> labMaps = generate_all<standardHasher>(oph_sigs, lvls, graph);
    
    printgraph(graph, labMaps, "../results/All_lvls_analyzed.dot");
    
    std::cout <<"Each level in one go can be seen in ../results/All_lvls_analyzed.dot. \n";
    
    std::vector<std::unordered_map<size_t, const std::vector<size_t>*>> clusts = get_clusters(labMaps);
    
    size_t union_size = get_union_size(fracmin_sigs);
    size_t sum_size = 0;
    for(std::vector<std::uint64_t>& sketch : fracmin_sigs) sum_size += sketch.size();
    size_t t_max = (sum_size + union_size)/(2*64*s);
    std::vector<double> fcorrs = compute_fcorrs(0.01, 4);

    for(double f : fcorrs) std::cout << f << "\n";

    std::vector<std::vector<size_t>> buckets = std::get<0>(binning(labMaps, clusts, fracmin_sigs, s, 64, t_max, fcorrs));

    auto get_pointers = [&](const std::vector<std::size_t>& bin){
        std::vector<const std::vector<std::uint64_t>*> ptrs;
        ptrs.reserve(bin.size());

        for(size_t index : bin) ptrs.push_back(&fracmin_sigs[index]);

        return ptrs;
    };

    std::cout << "t_max was : " << t_max << "\n";
    for(size_t i = 0; i < buckets.size(); i++){
        if(!buckets[i].empty()) std::cout << "Bucket : " << i << " Has fracmin Union size : " << get_union_size_ptr(get_pointers(buckets[i])) << " And estimated Union size of " << static_cast<size_t>(get_union_size_ptr(get_pointers(buckets[i]))/s) << "\n";
    }
    
    lemon::ListGraph buckets_visualized;
    std::unordered_map<std::vector<size_t>, lemon::ListGraph::Node, standardHasher> LabelMap;
    construct_graph(buckets, buckets_visualized, LabelMap);

    printgraph(buckets_visualized, LabelMap, "../results/Buckets_binned.dot");

    std::cout <<"The resulting binning can be seen in ../results/Buckets_binned.dot. \n";

    auto print_ibf = [&](const std::string& filepath, const IBF& ibf){
        std::vector<std::string> colors = {
            "red", "blue", "green", "orange", "purple", "cyan", "magenta", "yellow", "brown", "pink"
        };

        std::ofstream dot(filepath);
        dot << "graph G {\n" << "  node [style = filled];\n";
        for(size_t b = 0; b < ibf.size(); b++){
            std::string color = colors[b % colors.size()];
            dot << "  " << b << " [label=\"" << printout(ibf[b]) << "\", color=\"" << color << "\"];\n";
        }
        dot << "}\n";
        dot.close();
    };

    std::vector<std::uint64_t> big_sig(1000000);
    std::iota(big_sig.begin(),big_sig.end(),1000000);
    std::get<1>(sigs).push_back(big_sig);
    std::get<0>(sigs).push_back(std::vector<std::uint64_t>(256,1));
    std::get<2>(sigs)[625] = "/path/to/file/625.fasta";
    auto full_hibf = generate_hibf<standardHasher>(sigs, lvls, s, 0.01, 4, refinements,1);
    std::unordered_map<size_t, std::string>& seq_to_file = std::get<2>(sigs);


    auto print_hibf = [&](const std::vector<std::vector<IBF>>& hibf_levels_){
        std::filesystem::path root_dir = "../results/HIBF";
        std::filesystem::create_directories(root_dir);

        print_ibf((root_dir / "root_IBF.dot").string(), hibf_levels_[0][0]);

        for(size_t lvl = 1; lvl < hibf_levels_.size(); lvl++){
            std::filesystem::path lvl_dir = root_dir / ("lvl" + std::to_string(lvl));
            std::filesystem::create_directories(lvl_dir);

            for(size_t index = 0; index < hibf_levels_[lvl].size(); index++){
                std::string filename = "IBF." + std::to_string(lvl) + "." + std::to_string(index) + ".dot";
                print_ibf((lvl_dir / filename).string(), hibf_levels_[lvl][index]);
            }
        }
    };

    print_hibf(std::get<0>(full_hibf));
    std::filesystem::path root_dir = "../results/HIBF";
    std::ofstream header_out(root_dir / "Test_Header.txt");
    write_linkage(header_out, seq_to_file);
    write_config(header_out, "/test/path/", 20, 24, 0.01, 4, 1024);
    write_header(header_out, std::get<0>(full_hibf), std::get<3>(full_hibf), std::get<4>(full_hibf));
    write_content(header_out, std::get<2>(full_hibf));
    header_out.close();

    return 0;
}