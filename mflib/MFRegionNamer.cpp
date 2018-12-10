#include "MFGraph.hpp"
#include "MFRegionNamer.hpp"

using namespace lemon;

namespace methylFlow {
  MFRegionNamer::MFRegionNamer ( MFGraph * g) : mfGraph(g), curId(0)
  {
  }

  MFRegionNamer::~MFRegionNamer()
  {
  }

  void MFRegionNamer::reach(const ListDigraph::Node & node) {
    curId++;
    buffer.str("");
    buffer << curId;
    mfGraph->rename_node(node, buffer.str());
  }
}
