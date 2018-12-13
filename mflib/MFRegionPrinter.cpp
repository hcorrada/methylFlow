#include <iostream>

#include "MFGraph.hpp"
#include "MFRegionPrinter.hpp"
#include "MethylRead.hpp"
#include "lemon/maps.h"

using namespace lemon;

namespace methylFlow {
    
    MFRegionPrinter::MFRegionPrinter( MFGraph * g,
                                     std::ostream * ostream,
                                     const int cid,
                                      const float scale, std::string chr , const bool gr) : mfGraph(g),
    outstream(ostream),
    componentID(cid),
    scale_mult(scale),
                                                                                                    chromosome(chr),
                                                                                                    graph_only(gr)
    {
    }
    
    MFRegionPrinter::~MFRegionPrinter()
    {
    }
    
    std::ostream & MFRegionPrinter::getstream()
    {
        return *outstream;
    }
    
    void MFRegionPrinter::reach(const ListDigraph::Node & node) {
        IdMap<ListDigraph, ListDigraph::Node> idmap(mfGraph->get_graph());
        
#ifndef NDEBUG
        std::cout << "Region Printer: printing node " << mfGraph->node_name(node) << std::endl;
#endif
        if (node == mfGraph->get_source() || node == mfGraph->get_sink()) return;
        MethylRead * read = mfGraph->read(node);
        if (!read) return;
        
        getstream() << chromosome << "\t"<< read->start() << "\t" << read->end();
        getstream() << "\t" << componentID << "\t" << mfGraph->node_name(node);
        getstream() << "\t" << mfGraph->coverage(node);
        getstream() << "\t" << (mfGraph->isNormalized() ? mfGraph->normalized_coverage(node) : 0.);

        if (!graph_only) {
          getstream() << "\t" << mfGraph->expected_coverage(node, scale_mult);
        }

        getstream() << "\t" << read->getMethString() << std::endl;
    }
} // namespace methylFlow
