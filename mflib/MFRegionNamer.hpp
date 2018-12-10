#include <lemon/bfs.h>

using namespace lemon;

#ifndef MFREGIONNAMER_H
#define MFREGIONNAMER_H

namespace methylFlow {
  class MFGraph;

  class MFRegionNamer : public BfsVisitor<ListDigraph> {
    friend class MFGraph;

  public:
    MFRegionNamer(MFGraph * g);
    ~MFRegionNamer();
    void reach (ListDigraphBase::Node const &node);

  protected:
    MFGraph *mfGraph;
    int curId;
    std::stringstream buffer; 

  };
}
#endif // MFREGIONNAMER_H
