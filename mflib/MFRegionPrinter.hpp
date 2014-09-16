#include <lemon/bfs.h>

using namespace lemon;

#ifndef MFREGIONPRINTER_H
#define MFREGIONPRINTER_H

namespace methylFlow {
    class MFGraph;
    
    class MFRegionPrinter : public BfsVisitor<ListDigraph> {
        
        friend class MFGraph;
        
    public:
        MFRegionPrinter(MFGraph * g, std::ostream * ostream, const int cid, const float scale_mult, int chr);
        ~MFRegionPrinter();
        std::ostream & getstream();
        void reach (ListDigraphBase::Node const &node);
    protected:
        MFGraph *mfGraph;
        std::ostream * outstream;
        int componentID;
        float scale_mult;
        int chromosome;
    };
    
} // namespace methylFlow

#endif // MFREGIONPRINTER_H
