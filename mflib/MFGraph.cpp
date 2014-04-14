#include <iostream>
#include <stack>

#include <lemon/bfs.h>
#include <lemon/path.h>

#include "MFGraph.hpp"
#include "MFRegionPrinter.hpp"

namespace methylFlow {

  MFGraph::MFGraph() : mfGraph(), nodeName_map(mfGraph), coverage_map(mfGraph), normalized_coverage_map(mfGraph), read_map(mfGraph),
		     flow_map(mfGraph), effectiveLength_map(mfGraph),
		     source(), sink(), fake(mfGraph, false),
		       parentless(mfGraph, false), childless(mfGraph, false), is_normalized(false)
{
}

MFGraph::~MFGraph()
{
}

  void MFGraph::clear_graph()
  {
    std::vector<ListDigraph::Node> nodes;
    for (ListDigraph::NodeIt n(mfGraph); n != INVALID; ++n) {
      nodes.push_back(n);
    }

    for (std::vector<ListDigraph::Node>::iterator it = nodes.begin(); it != nodes.end(); ++it) {
      if (read_map[*it]) delete read_map[*it];
      mfGraph.erase(*it);
    }

  }

  ListDigraph::Node MFGraph::addNode(const std::string name, const int coverage, MethylRead *read)
{
  ListDigraph::Node n = mfGraph.addNode();
  nodeName_map[n] = name;
  coverage_map[n] = coverage;
  if (is_normalized) {
    normalized_coverage_map[n] = (float) coverage;
  }

  if (read) read->node = n;
  read_map[n] = read;

  return n;
}

ListDigraph::Arc MFGraph::addArc(const ListDigraph::Node u, const ListDigraph::Node v, const int length)
{
  ListDigraph::Arc arc = mfGraph.addArc(u, v);
  effectiveLength_map[arc] = (float) length;
  return arc;
}

void MFGraph::print_graph()
{
  MethylRead *read;
  for(ListDigraph::NodeIt n(mfGraph); n != INVALID; ++n) {
    std::cout << "Node: " << nodeName_map[n] << " cov: " << coverage_map[n];

    if (is_normalized) {
      std::cout << " normalized_cov: " << normalized_coverage_map[n];
    }

    if (read_map[n]) {
      read = read_map[n];
      std::cout << " " << read->getString();
    } 
    std::cout << std::endl;
  }
  std::cout << std::endl;

  for (ListDigraph::ArcIt arc(mfGraph); arc != INVALID; ++arc) {
    std::cout << "Arc: " << nodeName_map[mfGraph.source(arc)];
    std::cout << " -> " << nodeName_map[mfGraph.target(arc)];
    std::cout << " len: " << effectiveLength_map[arc];
    std::cout << " flow: " << flow_map[arc] << std::endl;
  }
  std::cout << std::endl;
}

  const float MFGraph::expected_coverage(const ListDigraph::Node node, const float scale)  const {
    float out = 0.;
    for (ListDigraph::OutArcIt arc(mfGraph, node); arc != INVALID; ++arc) {
      out += flow_map[arc] * effectiveLength_map[arc] / scale;
    }
    
    return out;
  }

  void MFGraph::print_regions( std::ostream & region_stream,
			       const float scale_mult, const int componentId )
  {
    MFRegionPrinter regionPrinter(this, &region_stream, componentId, scale_mult);
    BfsVisit<ListDigraph, MFRegionPrinter, BfsVisitDefaultTraits<ListDigraph> > bfs(mfGraph, regionPrinter);

    bfs.run(source);
  }

  // assumes file is sorted by position
  int MFGraph::run( std::istream & instream,
		    std::ostream & comp_stream,
		    std::ostream & patt_stream,
		    std::ostream & region_stream,
		    const float lambda,
		    const float scale_mult,
		    const bool verbose )
{
  std::string readid, rStrand, methStr, substStr;
  std::string input;

  int rPos, rLen;
  std::list<ListDigraph::Node> activeSet;
  
  ListDigraph::Node node;
  int rightMostPos = 0;

  const int READ_LIMIT = 100;
  #ifndef NDEBUG
  bool check_count = true;
  #else
  bool check_count = false;
  #endif

  int count = 0;
  int componentCount = 0;

  if (verbose) {
    std::cout << "[methylFlow] Reading from file " << std::endl; 
  }

  while (!check_count || count < READ_LIMIT) {
    std::getline( instream, input );
    if( !instream ) break; // checks end of file
    count++;

    // parse tab-separated line
    std::istringstream buffer(input);
    buffer >> readid >> rPos >> rLen >> rStrand >> methStr >> substStr;
    if ( !buffer || !buffer.eof() ) {
      std::cerr << "[methylFlow] Error parsing input" << std::endl;
      return -1;
    }

    // construct object with read info
     MethylRead * m = new MethylRead(rPos, rLen);
     if (methStr != "*")
       m->parseMethyl(methStr);

     #ifndef NDEBUG
     std::cout << "read: " << readid << " " << m->getString() << std::endl;
     #endif

     // does this read start after the rightMost end position?
     if (m->start() > rightMostPos) {
       // clear active reads if necessary
       if (!activeSet.empty())
	 activeSet.clear();

       // process this connected component
       if (count > 1) {
	 componentCount++;
	 if (verbose) {
	   std::cout << "[methylFlow] Processing component " << componentCount << std::endl;
	 }
	 
	 run_component( componentCount, 
		       comp_stream, 
			patt_stream, 
			region_stream,
			lambda, 
			scale_mult,
			verbose );
	 clear_graph();
       }
     }

     // if no reads in active set, add the node to the graph
     if (activeSet.empty()) {
       node = addNode(readid, 1, m);
       activeSet.push_front(node);

       // update the right-most position
       rightMostPos = m->end();
       continue;
     }

     // add read to graph
     if (processRead(m, readid, &activeSet) && m->end() > rightMostPos) {
       rightMostPos = m->end();
     }
  }

  componentCount++;
  if (verbose) {
    std::cout << "[methylFlow] Processing last component " << componentCount << std::endl;
  }
  run_component( componentCount, 
		 comp_stream, 
		 patt_stream, 
		 region_stream,
		 lambda,
		 scale_mult,
		 verbose );
  clear_graph();
  return 0;
}

  bool MFGraph::processRead(MethylRead *read, const std::string readid, std::list<ListDigraph::Node> *pactiveSet)
  {
    #ifndef NDEBUG
    std::cout << "processing read " << readid << std::endl;
    std::cout << "current active set: ";
    for (std::list<ListDigraph::Node>::iterator it = pactiveSet->begin(); it != pactiveSet->end(); ++it) {
      std::cout << nodeName_map[*it] << " ";
    }
    std::cout << std::endl;
    #endif

    // search for identical/superread/subread from left side of active set
    for (std::list<ListDigraph::Node>::iterator it = pactiveSet->begin(); it != pactiveSet->end(); ) {
      ListDigraph::Node active_node = *it;
      if (read_map[active_node]->start() > read->start()) {
	// we're done here
	break;
      }

      ReadComparison cmp = read_map[active_node]->compare(read);

      #ifndef NDEBUG
      std::cout << "checking for identical " << nodeName_map[active_node] << std::endl;
      #endif

      bool erased = false;
      switch(cmp) {
      case IDENTICAL:
      case SUBREAD:
	// increase coverage of corresponding node
	// stop looking through active set
	coverage_map[active_node] += 1;
	#ifndef NDEBUG
	std::cout << "found!" << std::endl;
	#endif
	return false;
      case SUPERREAD:
	// replace read info of corresponding node
	// stop looking through active set
	coverage_map[active_node] += 1;
	if (read_map[active_node]) delete read_map[active_node];
	read_map[active_node] = read;
	#ifndef NDEBUG
	std::cout << "found!" << std::endl;
	#endif
	return false;
      case METHOVERLAP:
      case OVERLAP:
	// keep going
	break;
      case NONE:
	// remove node from active set since no more possible overlaps to be found
	#ifndef NDEBUG
	std::cout << "removing from active set " << nodeName_map[active_node] << std::endl;
	#endif
	it = pactiveSet->erase(it);
	erased = true;
      default:
	break;
      }
      if (!erased) ++it;
    }


    #ifndef NDEBUG
    std::cout << "new node added" << std::endl;
    #endif

    // didn't find identical/superread/subread
    // so we need a new node
    ListDigraph::Node new_node = addNode(readid, 1, read);

    // now we search for all consistent overlaps from the right
    // if a consistent overlap is found, an arc is added
    // all reachable nodes are marked to avoid extra arc
    ListDigraph::NodeMap<bool> reachable(mfGraph, false);
    for (std::list<ListDigraph::Node>::reverse_iterator rit = pactiveSet->rbegin(); rit != pactiveSet->rend(); ++rit) {
      ListDigraph::Node active_node = *rit;
      #ifndef NDEBUG
      std::cout << "comparing node " << nodeName_map[active_node] << std::endl;
      #endif

//      if (read_map[active_node]->start() >= read->start()) {
//	// we're done
//	#ifndef NDEBUG
//	std::cout << "got to left-most point" << std::endl;
//	#endif
//	break;
//      }
//
      if (reachable[active_node]) {
	// this node is reachable so ignore
	#ifndef NDEBUG
	std::cout << "reachable" << std::endl;
	#endif
	continue;
      }

      ReadComparison cmp = read_map[active_node]->compare(read);
      switch(cmp) {
      case METHOVERLAP:
	#ifndef NDEBUG
	std::cout << "consistent overlap" << std::endl;
	#endif

	// add arc since this is a consistent overlap to a node that is not reached 
	addArc(active_node, new_node, read->start() - read_map[active_node]->start());
	
	// now mark all new reachable nodes
	// scanning active set from left
	for (std::list<ListDigraph::Node>::iterator it = pactiveSet->begin(); &*it != &*rit; ++it) {
	  ListDigraph::Node search_node = *it;
	  
	  #ifndef NDEBUG
	  std::cout << "is this node reachable" << nodeName_map[search_node] << std::endl;
	  #endif

	  if (reachable[search_node]) {
	    // this node is already reachable so skip
	    #ifndef NDEBUG
	    std::cout << "it is, skip" << std::endl;
	    #endif
	    continue;
	  }

	  // see if we can reach this node
	  #ifndef NDEBUG
	  std::cout << "try bfs" << std::endl;
	  #endif

	  Bfs<ListDigraph> bfs(mfGraph);
	  reachable[search_node] = bfs.run(search_node, new_node);
	  if (reachable[search_node]) {
	    #ifndef NDEBUG
	    std::cout << "reached! mark all the rest" << std::endl;
	    #endif

	    // we reached it, so
	    // now mark all nodes in the path
	    Path<ListDigraph> path = bfs.path(new_node);
	    for (PathNodeIt<Path<ListDigraph> > reachable_node(mfGraph, path); reachable_node != INVALID; ++reachable_node) {
	      #ifndef NDEBUG
	      std::cout << nodeName_map[reachable_node] << " ";
	      #endif

	      reachable[reachable_node] = true;
	    }
	    #ifndef NDEBUG
	    std::cout << std::endl;
	    #endif
	  }
	}
	break;
      case OVERLAP:
	#ifndef NDEBUG
	std::cout << "inconsistent overlap, keep going" << std::endl;
	#endif
	// do not do anything
	break;
      case IDENTICAL:
      case SUBREAD:
      case SUPERREAD:
      case NONE:
	// we should never be here
	std::cout << "Error adding read" << std::endl;
	if (read_map[new_node]) delete read_map[new_node];
	mfGraph.erase(new_node);
	return false;
      }
    }
    
    // since we're here we added a node
    // now add it to the active set
    pactiveSet->push_back(new_node);

    #ifndef NDEBUG
    std::cout << std::endl << std::endl;
    #endif
    return true;
  }

  void MFGraph::merge_nodes(ListDigraph::Arc arc)
  {
    ListDigraph::Node v = mfGraph.source(arc);
    ListDigraph::Node u = mfGraph.target(arc);


    coverage_map[v] += coverage_map[u];
    normalized_coverage_map[v] += normalized_coverage_map[u];
    read_map[v]->merge(read_map[u]);
    if (read_map[u]) delete read_map[u];
    mfGraph.contract(v, u);
  }

  void MFGraph::merge_chains()
  {
    std::stack<ListDigraph::Node> stack;
    stack.push(source);
    bool doagain = false;
    ListDigraph::NodeMap<bool> reached(mfGraph, false);

    while (!stack.empty()) {
      // grab and pop next node to process
      ListDigraph::Node curNode = stack.top();
      stack.pop();
      
      // check if out degree == 1 (and not a fake node)
      if (countOutArcs(mfGraph, curNode) == 1 && !fake[curNode]) {
	// grab the target of single arc
	ListDigraph::OutArcIt arc(mfGraph, curNode);
	ListDigraph::Node otherNode = mfGraph.target(arc);

	// check if target has single parent and is not fake
	if (countInArcs(mfGraph, otherNode) == 1 && !fake[otherNode]) {
	  // we can merge these nodes
	  doagain = true;
#ifndef NDEBUG
	  std::cout << "Merging nodes: " << nodeName_map[curNode] << " " << nodeName_map[otherNode] << std::endl;
#endif

	  // add coverage to current node
	  coverage_map[curNode] += coverage_map[otherNode];
	  normalized_coverage_map[curNode] += normalized_coverage_map[otherNode];

	  // merge methylation patterns
	  read_map[curNode]->merge(read_map[otherNode]);

	  // delete read object for target node
	  if (read_map[otherNode]) delete read_map[otherNode];

	  // connect current node to children of target node
	  for (ListDigraph::OutArcIt arc2(mfGraph, otherNode); arc2 != INVALID; ++arc2) {
	    addArc(curNode, mfGraph.target(arc2), 1);
	  }
	  // remove the target node
	  mfGraph.erase(otherNode);
	}
      }

      if (doagain) {
	// push current node to stack again
	stack.push(curNode);
	doagain = false;
      } else {
	reached[curNode] = true;
	// push children of current node to stack
	for (ListDigraph::OutArcIt arc(mfGraph, curNode); arc != INVALID; ++arc) {
	  ListDigraph::Node otherNode = mfGraph.target(arc);
	  if (!reached[otherNode]) stack.push(otherNode);
	}
      }
    }

    // fix the effective lengths
    for (ListDigraph::ArcIt arc(mfGraph); arc != INVALID; ++arc) {
      ListDigraph::Node source = mfGraph.source(arc);
      ListDigraph::Node target = mfGraph.target(arc);

      if (read_map[source] && read_map[target]) {
	effectiveLength_map[arc] = read_map[target]->start() - read_map[source]->start();
      }
    }

    // remove source and sink
    mfGraph.erase(source);
    mfGraph.erase(sink);
  }

  const int MFGraph::total_coverage() const
  {
    int out = 0;
    for (ListDigraph::NodeIt n(mfGraph); n != INVALID; ++n) {
      out += coverage_map[n];
    }
    return out;
  }

  const float MFGraph::total_normalized_coverage() const
  {
    float out = 0.;
    for (ListDigraph::NodeIt n(mfGraph); n != INVALID; ++n) {
      out += normalized_coverage_map[n];
    }
    return out;
  }
  
  const float MFGraph::total_flow() const
  {
    float total_flow = 0.;
    for (ListDigraph::InArcIt arc(mfGraph, sink); arc != INVALID; ++arc) {
      total_flow += flow_map[arc];
    }
    return total_flow;
  }

  const int MFGraph::run_component( const int componentID,
				    std::ostream & comp_stream,
				    std::ostream & patt_stream,
				    std::ostream & region_stream,
				    const float lambda,
				    const float scale_mult, 
				    const bool verbose ) 
  {
    #ifndef NDEBUG
    std::cout << "starting run_component" << std::endl;
    print_graph();
    #endif

    normalize_coverage();
    #ifndef NDEBUG
    std::cout << "normaliza complete" << std::endl;
    print_graph();
    #endif

    preprocess();
    #ifndef NDEBUG
    std::cout << "preprocess complete" << std::endl;
    print_graph();
    #endif

    merge_chains();
    #ifndef NDEBUG
    std::cout << "merge complete" << std::endl;
    print_graph();
    #endif
  
    if (verbose) {
      std::cout << "[methylFlow] Component " << componentID << " regions created" << std::endl;
    }
    // solve
    preprocess();
    solve( lambda, scale_mult );

    if (verbose) {
      std::cout << "[methylFlow] Component " << componentID << " estimation complete. Writing regions to file." << std::endl;
    }
    print_regions( region_stream, scale_mult, componentID );

    //#ifndef NDEBUG
    print_graph();
    //#endif

    // decompose
    float tflow = total_flow();
    int tcov = total_coverage();
    int npatterns = decompose(componentID, patt_stream);

    if (verbose) {
      std::cout << "[methylFlow] Component " << componentID << " wrote " << npatterns << " patterns to file." << std::endl;
    }

    int start = read(source)->start() + 1;
    int end = read(sink)->end();
    
    comp_stream << start << "\t" << end;
    comp_stream << "\t" << componentID << "\t" << npatterns;
    comp_stream << "\t" << tcov << "\t" << tflow << std::endl;
    std::cout << "[methylFlow] Finished processing component " << componentID << std::endl;
    return 0;
  }

} // namespace MethylFlow

