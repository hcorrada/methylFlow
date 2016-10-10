#include <iostream>
#include <string>
#include <stack>
#include <climits>

#include <time.h>
#include <stdio.h>
//#include <chrono> // for high_resolution_clock
#include <lemon/bfs.h>
#include <lemon/path.h>

#include <sys/time.h>


#include "MFGraph.hpp"
#include "MFRegionPrinter.hpp"
#include "MFCpgEstimator.hpp"

namespace methylFlow {

  MFGraph::MFGraph() : mfGraph(),
                       nodeName_map(mfGraph),
                       coverage_map(mfGraph),
                       normalized_coverage_map(mfGraph),
                       read_map(mfGraph),
                       flow_map(mfGraph), effectiveLength_map(mfGraph),
                       source(), sink(), fake(mfGraph, false),
                       parentless(mfGraph, false), childless(mfGraph, false), is_normalized(false) {
  }

  int findLength(std::string QNAME){
    std::size_t found;
    std::size_t curStringOffset = 0;

    found = QNAME.find("=", curStringOffset);
    return atoi(QNAME.substr( found +1).c_str());
  }

  MFGraph::~MFGraph() {
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

  ListDigraph::Node MFGraph::addNode(const std::string name, const int coverage, MethylRead *read) {
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

  ListDigraph::Arc MFGraph::addArc(const ListDigraph::Node u, const ListDigraph::Node v, const int length) {
    ListDigraph::Arc arc = mfGraph.addArc(u, v);
    effectiveLength_map[arc] = (float) length;
    return arc;
  }

  int MFGraph::get_graph_size() {
      int graphSize  = 0;
      for(ListDigraph::NodeIt n(mfGraph); n != INVALID; ++n) {
          graphSize++;
      }
      return graphSize;
  }
    
  void MFGraph::print_graph() {
      return;
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

  void MFGraph::print_regions(std::ostream & region_stream,
                              const float scale_mult,
                              const int componentId,
                              std::string chr )
  {
    MFRegionPrinter regionPrinter(this, &region_stream, componentId, scale_mult, chr);
    BfsVisit<ListDigraph, MFRegionPrinter, BfsVisitDefaultTraits<ListDigraph> > bfs(mfGraph, regionPrinter);
    bfs.run(source);
  }

  // assumes file is sorted by position
  int MFGraph::run(std::istream & instream,
                   std::ostream & comp_stream,
                   std::ostream & patt_stream,
                   std::ostream & region_stream,
                   std::ostream & cpg_stream,
                   std::string chr,
                   const long start,
                   const long end,
                   const bool flag_SAM,
                   const float lambda,
                   const float scale_mult,
                   const float epsilon,
                   const bool verbose,
                   const bool verboseTime,
                   const bool pctselect) {

    if (verbose) {
      
      std::cout << "start = " << start << std::endl;
      std::cout << "end = " << end << std::endl;
    }
#ifndef NDEBUG
    std::cout << "start = " << start << std::endl;
    std::cout << "end = " << end << std::endl;
#endif

    std::vector<MethylRead*> mVector;
    std::string readid, rStrand, methStr, substStr;
    std::string QNAME, RNAME, CIGAR, RNEXT, SEQ, QUAL, NM, XX, XM, XR, XG;
    int FLAG, POS, MAPQ, PNEXT, TLEN;
    std::string input;
    struct timeval tvalBefore, tvalAfter, tvalStartReadProcess, tvalEndReadProcess;

    int rPos, rLen;
    std::list<ListDigraph::Node> activeSet;

    ListDigraph::Node node;
    int rightMostPos = 0;

    const long READ_LIMIT = LONG_MAX;
#ifndef NDEBUG
    bool check_count = true;
#else
    bool check_count = false;
#endif

    long count = 0;
    long componentCount = 0;


    std::cout << "[methylFlow] Starting methylFlow... " << std::endl;

    // print headers to output files
    comp_stream << "chr\tstart\tend\tcid\tnnode\tnpatterns\ttotal_coverage\ttotal_flow\n";
    patt_stream << "chr\tstart\tend\tcid\tpid\tabundance\tmethylpat\tregions\n";
    region_stream << "chr\tstart\tend\tcid\trid\traw_coverage\tnorm_coverage\texp_coverage\tmethylpat\n";
    cpg_stream << "chr\tpos\tCov\tMeth\n";

    if(flag_SAM){
      while (std::getline(instream, input)){
        if(input.size() && input[0] !='@') break;
        if (verbose) {
          std::cerr << "[methylFlow] 0 Discarding lines start with " << input << std::endl;
        }
      }
    }

    std::string lastChr = "";
      if (verboseTime) {
          gettimeofday (&tvalStartReadProcess, NULL);
      }
    int maxActiveSetSize = 0;

    while (!check_count || count < READ_LIMIT) {
      if (maxActiveSetSize < activeSet.size()) {
        maxActiveSetSize = activeSet.size();
      }
        
      if(count > 0 || !flag_SAM)
        std::getline( instream, input );

      count++;
      if( !instream ) break; // checks end of file

      MethylRead * m;

      if(!flag_SAM){
        // parse tab-separated line
        std::istringstream buffer(input);
        buffer >> readid >> rPos >> rLen >> rStrand >> methStr >> substStr;
        if ( !buffer || !buffer.eof() ) {
          std::cerr << "[methylFlow] Error parsing tsv input" << std::endl;
          return -1;
        }
        if (rPos < start) {
          continue;
        }
        if (rPos + rLen > end) {
          break;
        }
        // construct object with read info
        m = new MethylRead(rPos, rLen);
        if (methStr != "*"){
          m->parseMethyl(methStr);
          mVector.push_back(m);
        }

// #ifndef NDEBUG
//         std::cout << "read: " << readid << " " << m->getString() << std::endl;
// #endif
      }

      else if(flag_SAM){
        //parse SAM format
        std::istringstream buffer(input);
        buffer >> QNAME >> FLAG >> RNAME >> POS >> MAPQ >> CIGAR >> RNEXT >> PNEXT >> TLEN >> SEQ >> QUAL >> NM  >> XX >> XM >> XR >> XG;

// #ifndef NDEBUG
//         std::cout << XX << "\t" << XG << std::endl;
// #endif
        if ( !buffer || !buffer.eof() ) {
          std::cerr << "[methylFlow] Error parsing SAM input" << std::endl;
          return -1;
        }
          if (verboseTime) {
              gettimeofday (&tvalBefore, NULL);
          }

        //parse chr name
        std::string::size_type sz;
        std::string chromosome;
        if (RNAME[0] == 'c' || RNAME[0] == 'C')
          chromosome =  RNAME.substr(3);
        else
          chromosome = RNAME;
        //                chr = atoi(chromosome.c_str());
        chr = chromosome;

// #ifndef NDEBUG
//         std::cout << "chr " << chr << std::endl;
//         std::cout << "str " << XM << std::endl;
//         std::cout << "pos " << POS << std::endl;
// #endif

        rPos = POS;
        rLen = SEQ.length();

        if (rPos < start) {
          continue;
        }
        if (rPos + rLen > end) {
          break;
        }
        //rLen = findLength(QNAME);

        //if (verbose) std::cout << "rLen " << rLen << std::endl;

// #ifndef NDEBUG
//         std::cout << "rLen " << rLen << std::endl;
// #endif
        // construct object with read info
        m = new MethylRead(rPos, rLen);
        if (methStr != "*"){
          m->parseXMtag(XM,CIGAR);

          mVector.push_back(m);

//           #ifndef NDEBUG
//           std::cout << "start = " << m->start() << std::endl;
//           std::cout << "end = " << m->end() << std::endl;
// #endif
        }

// #ifndef NDEBUG
//         std::cout << "read: " << readid << " " << m->getString() << std::endl;
// #endif
      }

      // we assume the input file is sorted

      // does this read start after the rightMost end position?
      if (m->start() > rightMostPos || chr != lastChr) {
        if (verbose) {
          std::cout << "chr " << chr << std::endl;
          std::cout << "[methylFlow] last chr = " << lastChr << ", chr = " << chr << std::endl;
        }

        lastChr = chr;
        // clear active reads if necessary
        if (!activeSet.empty()){
          activeSet.clear();
        }

        // process this connected component
        if (count > 1) {
          componentCount++;

          if (verbose) {
            std::cout << "[methylFlow] Processing component " << componentCount << std::endl;
            std::cout << "[methylFlow] Read number " << count << std::endl;
            std::cout << "[methylFlow] start read " << m->start() << std::endl;
            std::cout << "[methylFlow] rightMostPos " << rightMostPos << std::endl;
          }
            if (verboseTime) {
                gettimeofday (&tvalEndReadProcess, NULL);
                std::cout << "Time in miliseconds process component:\t" << componentCount << "\t" << get_graph_size() << "\t" << ((tvalEndReadProcess.tv_sec - tvalStartReadProcess.tv_sec)*1000  + tvalEndReadProcess.tv_usec/1000) - tvalStartReadProcess.tv_usec/1000 << std::endl;
                gettimeofday (&tvalBefore, NULL);
            }
          //std::cout << "[methylFlow] Processing component " << componentCount << std::endl;
          //std::cout << "[methylFlow] Read number " << get_graph_size() << std::endl;

          // count = 0;
          run_component(componentCount,
                        comp_stream,
                        patt_stream,
                        region_stream,
                        cpg_stream,
                        chr,
                        flag_SAM,
                        lambda,
                        scale_mult,
                        epsilon,
                        verbose,
                        verboseTime,
                        pctselect );
            if (verboseTime) {

                gettimeofday (&tvalAfter, NULL);
                std::cout << "Time in miliseconds run component:\t" << componentCount << "\t" << get_graph_size() << "\t" << ((tvalAfter.tv_sec - tvalBefore.tv_sec)*1000  + tvalAfter.tv_usec/1000) - tvalBefore.tv_usec/1000 << std::endl;
            }
          clear_graph();
            if (verboseTime) {
                gettimeofday (&tvalStartReadProcess, NULL);
        
            }
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
      if (processRead(m, readid, &activeSet, verboseTime) && m->end() > rightMostPos) {
        rightMostPos = m->end();
      }
        if (maxActiveSetSize < activeSet.size()) {
            maxActiveSetSize = activeSet.size();
        }
        
    }

      
    componentCount++;
    //std::cout << "[methylFlow] Processing last component (no. " << componentCount << ")" << std::endl;
    if (verbose) {
      std::cout << "[methylFlow] Read number " << count << std::endl;
    }
      
      if (verboseTime) {
          gettimeofday (&tvalEndReadProcess, NULL);
          std::cout << "Time in miliseconds process last component:\t" << componentCount << "\t" << get_graph_size() << "\t" << ((tvalEndReadProcess.tv_sec - tvalStartReadProcess.tv_sec)*1000  + tvalEndReadProcess.tv_usec/1000) - tvalStartReadProcess.tv_usec/1000 << std::endl;
  
          gettimeofday (&tvalBefore, NULL);
      }
    run_component(componentCount,
                  comp_stream,
                  patt_stream,
                  region_stream,
                  cpg_stream,
                  chr,
                  flag_SAM,
                  lambda,
                  scale_mult,
                  epsilon,
                  verbose,
                  verboseTime,
                  pctselect);
      if (verboseTime) {
          gettimeofday (&tvalAfter, NULL);
          std::cout << "Time in miliseconds run last component:\t" << componentCount << "\t" << get_graph_size() << "\t" << ((tvalAfter.tv_sec - tvalBefore.tv_sec)*1000  + tvalAfter.tv_usec/1000) - tvalBefore.tv_usec/1000 << std::endl;
          std::cout << "max active size:\t" << maxActiveSetSize << std::endl;
      }
    clear_graph();
    return 0;
  }


  bool MFGraph::processRead(MethylRead *read, const std::string readid, std::list<ListDigraph::Node> *pactiveSet, const bool verboseTime)  {
      struct timeval tvalBeforeIndentical, tvalAfterIndentical;
      if (verboseTime) {
          gettimeofday (&tvalBeforeIndentical, NULL);
      }
    // search for identical/superread/subread from left side of active set
    for (std::list<ListDigraph::Node>::iterator it = pactiveSet->begin(); it != pactiveSet->end(); ) {
      ListDigraph::Node active_node = *it;
      if (read_map[active_node]->start() > read->start()) {
        // we're done here
        break;
      }

      ReadComparison cmp = read_map[active_node]->compare(read);

      bool erased = false;
      switch(cmp) {
      case IDENTICAL:
      case SUBREAD:
        // increase coverage of corresponding node
        // stop looking through active set
        coverage_map[active_node] += 1;
        return false;

      case SUPERREAD:
        // replace read info of corresponding node
        // stop looking through active set
        coverage_map[active_node] += 1;
        if (read_map[active_node]) delete read_map[active_node];
        read_map[active_node] = read;
        return false;

      case METHOVERLAP:
      case OVERLAP:
        // keep going
        break;

      case NONE:
        // remove node from active set since no more possible overlaps to be found
        it = pactiveSet->erase(it);
        erased = true;
      default:
        break;
      }
      if (!erased) ++it;
    }
      if (verboseTime) {
          gettimeofday (&tvalAfterIndentical, NULL);
      }
    // didn't find identical/superread/subread
    // so we need a new node
    ListDigraph::Node new_node = addNode(readid, 1, read);

    // now we search for all consistent overlaps from the right
    // if a consistent overlap is found, an arc is added
    // all reachable nodes are marked to avoid extra arc
    ListDigraph::NodeMap<bool> reachable(mfGraph, false);
    for (std::list<ListDigraph::Node>::reverse_iterator rit = pactiveSet->rbegin(); rit != pactiveSet->rend(); ++rit) {
      ListDigraph::Node active_node = *rit;
      if (reachable[active_node]) {
        // this node is reachable so ignore
        continue;
      }

      ReadComparison cmp = read_map[active_node]->compare(read);
      switch(cmp) {
      case METHOVERLAP:
// #ifndef NDEBUG
//         std::cout << "consistent overlap" << std::endl;
// #endif
        // add arc since this is a consistent overlap to a node that is not reached
        addArc(active_node, new_node, read->start() - read_map[active_node]->start());

        // now mark all new reachable nodes
        // scanning active set from left
        for (std::list<ListDigraph::Node>::iterator it = pactiveSet->begin(); &*it != &*rit; ++it) {
          ListDigraph::Node search_node = *it;


          if (reachable[search_node]) {
            // this node is already reachable so skip
            continue;
          }

          // see if we can reach this node
          Bfs<ListDigraph> bfs(mfGraph);
          reachable[search_node] = bfs.run(search_node, new_node);
          if (reachable[search_node]) {
            // we reached it, so
            // now mark all nodes in the path
            Path<ListDigraph> path = bfs.path(new_node);
            for (PathNodeIt<Path<ListDigraph> > reachable_node(mfGraph, path); reachable_node != INVALID; ++reachable_node) {
              reachable[reachable_node] = true;
            }
          }
        }
        break;
      case OVERLAP:
// #ifndef NDEBUG
//         std::cout << "inconsistent overlap, keep going" << std::endl;
// #endif
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

// #ifndef NDEBUG
//     std::cout << std::endl << std::endl;
// #endif
    return true;
  }


  void MFGraph::merge_nodes(ListDigraph::Arc arc) {
    ListDigraph::Node v = mfGraph.source(arc);
    ListDigraph::Node u = mfGraph.target(arc);

    coverage_map[v] += coverage_map[u];
    normalized_coverage_map[v] += normalized_coverage_map[u];
    read_map[v]->merge(read_map[u]);
    if (read_map[u]) delete read_map[u];
    mfGraph.contract(v, u);
  }

  void MFGraph::merge_chains() {
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
    //    mfGraph.erase(source);
    //    mfGraph.erase(sink);
  }

  const int MFGraph::total_coverage() const {
    int out = 0;
    for (ListDigraph::NodeIt n(mfGraph); n != INVALID; ++n) {
      out += coverage_map[n];
    }
    return out;
  }

  const float MFGraph::total_normalized_coverage() const {
    float out = 0.;
    for (ListDigraph::NodeIt n(mfGraph); n != INVALID; ++n) {
      out += normalized_coverage_map[n];
    }
    return out;
  }

  const float MFGraph::total_flow() const {
    float total_flow = 0.;
    for (ListDigraph::InArcIt arc(mfGraph, sink); arc != INVALID; ++arc) {
      total_flow += flow_map[arc];
    }
    return total_flow;
  }

  const int MFGraph::run_component(const int componentID,
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
                                   const bool verboseTime,
                                   const bool pctselect ) {
#ifndef NDEBUG
    std::cout << "starting run_component" << std::endl;
    print_graph();
#endif

    normalize_coverage();
#ifndef NDEBUG
    std::cout << "normaliza complete" << std::endl;
    print_graph();
#endif

    add_terminals();
#ifndef NDEBUG
    std::cout << "preprocess complete" << std::endl;
    print_graph();
#endif

    MFCpgEstimator cpg_estimator(this, &cpg_stream, scale_mult);
    cpg_estimator.computeRaw();
    cpg_estimator.printRaw(chr, false);

    merge_chains();
#ifndef NDEBUG
    std::cout << "merge complete" << std::endl;
    print_graph();
#endif

    if (verbose) {
      std::cout << "[methylFlow] Component " << componentID << " regions created" << std::endl;
    }
      struct timeval tvalBeforeSolve, tvalAfterSolve;
      if (verboseTime) {
          // solve
          gettimeofday (&tvalBeforeSolve, NULL);
      }
    int res = solve( lambda, scale_mult, epsilon, verbose, verboseTime, pctselect );
      if (verboseTime) {
          gettimeofday (&tvalAfterSolve, NULL);
          std::cout << "Time in miliseconds run solve:\t" << componentID << "\t" << get_graph_size() << "\t" << ((tvalAfterSolve.tv_sec - tvalBeforeSolve.tv_sec)*1000  + tvalAfterSolve.tv_usec/1000) - tvalBeforeSolve.tv_usec/1000 << std::endl;
      }
      
    if (res) {
      std::cerr << "[methylFlow] Error solving" << std::endl;
      return res;
    }

    if (verbose) {
      std::cout << "[methylFlow] Component " << componentID << " estimation complete. Writing regions to file." << std::endl;
    }
    print_regions( region_stream, scale_mult, componentID, chr );

#ifndef NDEBUG
    print_graph();
#endif

    // decompose
    float tflow = total_flow();
    int tcov = total_coverage();

#ifndef NDEBUG
    std::cout << "befor decompose" << std::endl;
#endif
      struct timeval tvalBeforeDecompose, tvalAfterDecompose;
      if (verboseTime) {
          gettimeofday (&tvalBeforeDecompose, NULL);
      }
    int npatterns = decompose(componentID, patt_stream, chr);
      if (verboseTime) {
          gettimeofday (&tvalAfterDecompose, NULL);
          std::cout << "Time in miliseconds run decompose:\t" << componentID << "\t" << get_graph_size() << "\t" << ((tvalAfterDecompose.tv_sec - tvalBeforeDecompose.tv_sec)*1000  + tvalAfterDecompose.tv_usec/1000) - tvalBeforeDecompose.tv_usec/1000 << std::endl;
      }
#ifndef NDEBUG
    std::cout << "after decompose" << std::endl;
#endif

    if (verbose) {
      std::cout << "[methylFlow] Component " << componentID << " wrote " << npatterns << " patterns to file." << std::endl;
    }
    // we compute all the offsets from source!
    int start = read(source)->start()+1;
    int end = read(sink)->end();

    comp_stream << chr << "\t" << start << "\t" << end;
    comp_stream << "\t" << componentID << "\t" << get_graph_size() << "\t" << npatterns;
    comp_stream << "\t" << tcov << "\t" << tflow << std::endl;

    if (verbose) {
      std::cout << "[methylFlow] Finished processing component " << componentID << std::endl;
    }
    return 0;
  }
} // namespace MethylFlow
