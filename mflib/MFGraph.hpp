#include <string>
#include <algorithm>

#include <lemon/list_graph.h>
#include <lemon/maps.h>
#include <lemon/lp.h>

#include "MethylRead.hpp"

using namespace lemon;

#ifndef MFGRAPH_H
#define MFGRAPH_H

namespace methylFlow {
    class MFSolver;
    
    class MFGraph {
        friend class MFSolver;
        friend class MFRegionSolver;
        friend class MFCpgSolver;
        
    public:
        MFGraph();
        ~MFGraph();
        
        ListDigraph &get_graph();
        const ListDigraph &get_graph() const;
        
        const std::string &node_name(const ListDigraph::Node &node) const;
        std::string &node_name(const ListDigraph::Node &node);
        
        const int &coverage(const ListDigraph::Node &node) const;
        int &coverage(const ListDigraph::Node &node);
        
        const float &normalized_coverage(const ListDigraph::Node &node) const;
        float &normalized_coverage(const ListDigraph::Node &node);
        
        const MethylRead * read(const ListDigraph::Node &node) const;
        MethylRead * read(const ListDigraph::Node &node);
        
        const float &flow(const ListDigraph::Arc &arc) const;
        float &flow(const ListDigraph::Arc &arc);
        
        const int &effective_length(const ListDigraph::Arc &arc) const;
        int &effective_length(const ListDigraph::Arc &arc);
        
        ListDigraph::Node addNode(const std::string name, const int coverage, MethylRead *read = NULL);
        ListDigraph::Arc addArc(const ListDigraph::Node u, const ListDigraph::Node v, const int length);
        
        const ListDigraph::Node & get_source() const;
        const ListDigraph::Node & get_sink() const;
        
        const bool &isNormalized() const;
        
        // return total coverage
        const int total_coverage() const;
        
        // return total normalized coverage
        const float total_normalized_coverage() const;
        
        // return total flow
        const float total_flow() const;
        
        // return expected coverage
        const float expected_coverage(const ListDigraph::Node node, const float scale) const;
        
        // run reading from instream
        int run( std::istream & instream,
                std::ostream & comp_stream,
                std::ostream & patt_stream,
                std::ostream & region_stream,
                std::ostream & cpg_stream,
                std::string chr,
                const int start,
                const int end,
                const bool flag_SAM,
                const float lambda,
                const float scale_mult,
                const float epsilon,
                const bool verbose,
                const bool pctselect );
        
        // tsv file with readid, pos, length, strand (ignored), methylString, subString
        bool processRead(MethylRead *read, const std::string readid, std::list<ListDigraph::Node> *pactiveSet);
        
        // print the graph
        void print_graph();
        
        void print_regions( std::ostream & region_stream,
                           const float scale_mult,
                           const int componentId, std::string chr );
        
        // clear graph, delete pointers to read/region objects
        void clear_graph();

        
    protected:
        ListDigraph mfGraph;
        ListDigraph::NodeMap<std::string> nodeName_map;

        ListDigraph::NodeMap<int> coverage_map;
        ListDigraph::NodeMap<float> normalized_coverage_map;
        ListDigraph::NodeMap<MethylRead *> read_map;
        
        ListDigraph::ArcMap<float> flow_map;
        ListDigraph::ArcMap<int> effectiveLength_map;
        ListDigraph::Node source, sink;
        IterableBoolMap<ListDigraph, ListDigraph::Node> fake;
        ListDigraph::NodeMap<bool> parentless;
        ListDigraph::NodeMap<bool> childless;
        
    private:
        bool is_normalized;
        
        // run on current component
        const int run_component( const int componentID,
                                std::ostream & comp_stream,
                                std::ostream & patt_stream,
                                std::ostream & region_stream,
                                std::ostream & cpg_stream,
                                std::string chr,
                                const bool flag_SAM,
                                const float lambda,
                                const float scale_mult,
                                const float epsilon,
                                const bool verbose,
                                const bool pctselect );
        
        // merge chains in read overlap graph
        void merge_chains();
        
        // helper function to merge nodes
        void merge_nodes(ListDigraph::Arc arc);
        
        // add source and sink to graph
        void add_terminals();
        
        // normalize coverage per position
        void normalize_coverage();
        float calculate_median(std::vector<float> x);
        
        // add nodes for regularization penalty
        void regularize();
        
        // solve wrapper
        int solve(const float lambda, const float length_mult, const float epsilon, const bool verbose, const bool pctselect);
        
        // run decomposition algorithm
        // componentID: used for printing
        int decompose(const int componentID, std::ostream & patt_stream, std::string chr);
    };
    
    inline const ListDigraph &MFGraph::get_graph() const
    {
        return mfGraph;
    }
    
    inline ListDigraph &MFGraph::get_graph()
    {
        return mfGraph;
    }
    
    
    inline const std::string &MFGraph::node_name(const ListDigraph::Node &node) const
    {
        return nodeName_map[node];
    }
    
    inline std::string &MFGraph::node_name(const ListDigraph::Node &node)
    {
        return nodeName_map[node];
    }
    
    
    inline const int &MFGraph::coverage(const ListDigraph::Node &node) const
    {
        return coverage_map[node];
    }
    
    inline int &MFGraph::coverage(const ListDigraph::Node &node)
    {
        return coverage_map[node];
    }
    
    inline const float &MFGraph::normalized_coverage(const ListDigraph::Node &node) const
    {
        return normalized_coverage_map[node];
    }
    
    inline float &MFGraph::normalized_coverage(const ListDigraph::Node &node)
    {
        return normalized_coverage_map[node];
    }
    
    inline const MethylRead * MFGraph::read(const ListDigraph::Node &node) const
    {
        return read_map[node];
    }
    
    inline MethylRead * MFGraph::read(const ListDigraph::Node &node) 
    {
        return read_map[node];
    }
    
    inline const float &MFGraph::flow(const ListDigraph::Arc &arc) const
    {
        return flow_map[arc];
    }
    
    inline float &MFGraph::flow(const ListDigraph::Arc &arc)
    {
        return flow_map[arc];
    }
    
    inline const int &MFGraph::effective_length(const ListDigraph::Arc &arc) const
    {
        return effectiveLength_map[arc];
    }
    
    inline int &MFGraph::effective_length(const ListDigraph::Arc &arc)
    {
        return effectiveLength_map[arc];
    }
    
    inline const ListDigraph::Node & MFGraph::get_source() const
    {
        return source;
    }
    
    inline const ListDigraph::Node & MFGraph::get_sink() const
    {
        return sink;
    }
    
    inline const bool &MFGraph::isNormalized() const
    {
        return is_normalized;
    }
    
    template<typename V>
    struct DijkstraMinMaxOperationTraits {
        typedef V Value;
        
        static Value zero() {
            return static_cast<Value>(0);
        }
        
        static Value plus(const Value& left, const Value& right) {
            return std::max(left, right);
        }
        
        static Value less(const Value& left, const Value& right) {
            return left < right;
        }
    };
    
} // namespace methylFlow

#endif // MFGRAPH_H
