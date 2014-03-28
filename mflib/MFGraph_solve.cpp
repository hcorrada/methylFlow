#include <iostream>
#include <queue>

#include <lemon/lp.h>
#include <lemon/path.h>
#include <lemon/dijkstra.h>

#include "MFGraph.hpp"

namespace methylFlow {
void MFGraph::preprocess()
{
  int rightMostStart = -1;
  // find childless nodes
  for (ListDigraph::NodeIt n(mfGraph); n != INVALID; ++n) {
    if (countOutArcs(mfGraph, n) == 0) {
      childless[n] = true;
      if (read_map[n] && read_map[n]->start() > rightMostStart) {
	rightMostStart = read_map[n]->start();
      }
    }
  }

  int leftMostStart = rightMostStart;
  // find parentless nodes
  for (ListDigraph::NodeIt n(mfGraph); n != INVALID; ++n) {
    if (countInArcs(mfGraph, n) == 0) {
      parentless[n] = true;
      if (read_map[n] && read_map[n]->start() < leftMostStart) {
	leftMostStart = read_map[n]->start();
      }
    }
  }

  // now add the fake source and sink
  MethylRead *source_read = new MethylRead(leftMostStart - 1, 1);
  source = addNode("s", 0, source_read);
  fake[source] = true;
  for (ListDigraph::NodeIt n(mfGraph); n != INVALID; ++n) {
    if (parentless[n]) {
      addArc(source, n, read_map[n]->start() - source_read->start());
    }
  }

  MethylRead *sink_read = new MethylRead(rightMostStart + 1, 1);
  sink = addNode("t", 0, sink_read);
  fake[sink] = true;
  for (ListDigraph::NodeIt n(mfGraph); n != INVALID; ++n) {
    if (childless[n]) {
      addArc(n, sink, sink_read->start() - read_map[n]->start());
    }
  }
}

  float MFGraph::calculate_median(std::vector<float> x) {
    float median;
    size_t size = x.size();
    
    sort(x.begin(), x.end());
    if (size % 2 == 0) {
      median = (x[size / 2 - 1] + x[size / 2]) / 2.;
    } else {
      median = x[size / 2];
    }

    return median;
  }


  void MFGraph::normalize_coverage()
  {
    // put nodes in a priority queue by start position
    // to ensure we iterate in the proper order
    std::priority_queue<MethylRead*, std::vector<MethylRead*>, CompareReadStarts> position_queue;
    for (ListDigraph::NodeIt n(mfGraph); n != INVALID; ++n) {
      MethylRead *read = read_map[n];
      std::cout << read->start() << std::endl;
      position_queue.push(read);
    }

    std::vector<float> coverage_per_position;
    std::vector<ListDigraph::Node> current_nodes;

    int current_startpos = -1;
    float current_coverage = 0.;

    while (! position_queue.empty() ) {
      MethylRead *read = position_queue.top();
      position_queue.pop();

      ListDigraph::Node n = read->node;

      #ifndef DEBUG
      std::cout << "current node:" << read->start() << " " << coverage_map[n] << std::endl;
      #endif

      if (current_startpos == -1 ) {
	current_startpos = read->start();
	current_coverage = coverage_map[n];
	current_nodes.push_back(n);
	continue;
      }

      if (read->start() < current_startpos) {
	coverage_per_position.push_back(current_coverage);
	
	#ifndef NDEBUG
	std::cout << "new position!" << std::endl;
	std::cout << current_startpos << " " << current_coverage << std::endl;
	#endif

	// divide by position coverage
	for (std::vector<ListDigraph::Node>::iterator it = current_nodes.begin(); it != current_nodes.end(); ++it) {
	  ListDigraph::Node node = *it;
	  normalized_coverage_map[node] = (float) coverage_map[node] / current_coverage;
	  #ifndef NDEBUG
	  std::cout << coverage_map[node] << " " << normalized_coverage_map[node] << std::endl;
	  #endif
	}

	#ifndef NDEBUG
	std::cout << std::endl;
	#endif

	current_nodes.clear();
	current_nodes.push_back(n);
	current_coverage = (float) coverage_map[n];
	current_startpos = read->start();
      } if (read->start() == current_startpos) {
	current_coverage += (float) coverage_map[n];
	current_nodes.push_back(n);
      } else {
	std::cerr << "nodes out of order" << std::endl;
	return;
      }
    }

    // process last position
    coverage_per_position.push_back(current_coverage);
    for (std::vector<ListDigraph::Node>::iterator it = current_nodes.begin(); it != current_nodes.end(); ++it) {
      ListDigraph::Node node = *it;
      normalized_coverage_map[node] = (float) coverage_map[node] / current_coverage;
      #ifndef NDEBUG
      std::cout << coverage_map[node] << " " << normalized_coverage_map[node] << std::endl;
      #endif
    }

    float median_coverage = calculate_median(coverage_per_position);
    #ifndef DEBUG
    std::cout << "median coverage: " << median_coverage << std::endl;
    #endif

    for (ListDigraph::NodeIt n(mfGraph); n != INVALID; ++n) {
      normalized_coverage_map[n] *= median_coverage;
    }
    is_normalized = true;
  }

  int MFGraph::solve(const float lambda, const float length_mult)
{
  Lp lp;
  
  ListDigraph::NodeMap<Lp::Col> alpha(mfGraph);
  ListDigraph::NodeMap<Lp::Col> beta(mfGraph);
  ListDigraph::NodeMap<Lp::Col> nu(mfGraph);

  ListDigraph::ArcMap<Lp::Row> rows(mfGraph);
  ListDigraph::ArcMap<float> scaled_length(mfGraph);

  // scale the lengths
  for (ListDigraph::ArcIt arc(mfGraph); arc != INVALID; ++arc) {
    scaled_length[arc] = effectiveLength_map[arc] / length_mult;
  }

  // add the lambda edges
  int lamcnt = 0;
  for (ListDigraph::InArcIt arc(mfGraph, sink); arc != INVALID; ++arc, ++lamcnt) {
    std::stringstream nodename;
    nodename << "lambda_" << lamcnt;
    ListDigraph::Node newNode = addNode(nodename.str(), 0);
    ListDigraph::Node node = mfGraph.source(arc);
    ListDigraph::Arc newArc = addArc(node, newNode, 0);
    scaled_length[newArc] = scaled_length[arc];
    mfGraph.changeSource(arc, newNode);
    scaled_length[arc] = - lambda;
    childless[node] = false;
    childless[newNode] = true;
  }
  Lp::Expr obj;

  #ifndef NDEBUG
  std::cout << "running solver on " << countNodes(mfGraph) << " nodes" << std::endl;
  #endif

  for (ListDigraph::NodeIt v(mfGraph); v != INVALID; ++v) {
    #ifndef NDEBUG
    std::cout << "Processing node " << nodeName_map[v] << std::endl;
    #endif

    if (fake[v]) continue;

    // add node's alpha variable and bounds
    alpha[v] = lp.addCol();
    lp.colLowerBound(alpha[v], 0.);
    lp.colUpperBound(alpha[v], 1.);

    // add node's beta variable and bounds
    beta[v] = lp.addCol();
    lp.colLowerBound(beta[v], 0.);
    lp.colUpperBound(beta[v], 1.);

    // add node's nu variable
    nu[v] = lp.addCol();

    #ifndef NDEBUG
    std::cout << "LP vars added" << std::endl;
    #endif

    // add node's term in objective 
    obj += normalized_coverage_map[v] * (beta[v] - alpha[v]);

    #ifndef NDEBUG
    std::cout << "obj added" << std::endl;
    #endif
  }

  // bound nu variable for source targets
  for (ListDigraph::OutArcIt arc(mfGraph, source); arc != INVALID; ++arc) {
    ListDigraph::Node v = mfGraph.target(arc);
    rows[arc] = lp.addRow(nu[v] <= 0);
  }
  #ifndef NDEBUG
  std::cout << "nu bounds added" << std::endl;
  #endif

  // add sink constraints
  for (ListDigraph::InArcIt arc(mfGraph, sink); arc != INVALID; ++arc) {
    ListDigraph::Node v = mfGraph.source(arc);
    rows[arc] = lp.addRow(scaled_length[arc] * beta[v] - 
			  scaled_length[arc] * alpha[v] - nu[v] <= 0);
  }
  #ifndef NDEBUG
  std::cout << "sink constraints added" << std::endl;
  #endif

  // add remaining constraints (if not childless)
  for (IterableBoolMap<ListDigraph, ListDigraph::Node>::FalseIt v(fake); v != INVALID; ++v) {
    if (childless[v]) continue;
    
    for (ListDigraph::OutArcIt arc(mfGraph, v); arc != INVALID; ++arc) {
      ListDigraph::Node u = mfGraph.target(arc);
      if (u == INVALID) {
	std::cout << "error getting target from arc" << std::endl;
	return -1;
      }
      rows[arc] = lp.addRow(scaled_length[arc] * beta[v] - 
			    scaled_length[arc] * alpha[v] - nu[v] + nu[u] <= 0);
    }
    #ifndef NDEBUG
    std::cout << "constraints added" << std::endl;
    #endif
  }

  lp.obj(obj);
  lp.max();
  lp.solve();

  #ifndef NDEBUG
  std::cout << "Called solver" << std::endl;
  #endif

  if (lp.primalType() != Lp::OPTIMAL) {
    std::cout << "Did not find optimum" << std::endl;
    return -1;
  }

  // extract flows
  for (ListDigraph::ArcIt arc(mfGraph); arc != INVALID; ++arc) {
    #ifndef NDEBUG
    std::cout << "Extracting flow of arc: " << std::endl;
    std::cout << nodeName_map[mfGraph.source(arc)];
    std::cout << " -> ";
    std::cout << nodeName_map[mfGraph.target(arc)] << std::endl;
    #endif
    Lp::Row row = rows[arc];
    flow_map[arc] = lp.dual(row);
  }
  return 0;
}

  int MFGraph::decompose(const int componentID, std::ostream & patt_stream)
{
  // compute total flow
  float total_flow = this->total_flow();
  #ifndef NDEBUG
  std::cout << "total flow: " << total_flow << std::endl;
  #endif

  int flownum = 0;

  // iterate while residual flow 
  while (total_flow > 0.) {
    flownum++;
    // compute residual flow for each arc
    ListDigraph::ArcMap<int> residual_flow(mfGraph);
    for (ListDigraph::ArcIt arc(mfGraph); arc != INVALID; ++arc) {
      residual_flow[arc] = total_flow - flow_map[arc];
    }

    // run the min-max dijkstra algorithm
    ListDigraph::NodeMap<int> dist(mfGraph);
    Dijkstra<ListDigraph>
      ::SetOperationTraits<DijkstraMinMaxOperationTraits<int> > 
      ::Create dijkstra(mfGraph, residual_flow);
    dijkstra.distMap(dist);
    dijkstra.run(source, sink);

    // get the resulting path and it's flow
    Path<ListDigraph> shortestPath = dijkstra.path(sink);
    float path_flow = total_flow - dijkstra.dist(sink);

    // construct a meth fragment from path here
    // and remove path flow from each arc in path

    MethylRead pattern = MethylRead(*read_map[source]);
    for(Path<ListDigraph>::ArcIt arc(shortestPath); arc != INVALID; ++arc) {
      ListDigraph::Node s = mfGraph.source(arc);
      MethylRead *read = read_map[s];
      if (!read) continue;

      if (s != source) {
	pattern.merge(read_map[s]);
      }
      
      flow_map[arc] -= path_flow;

      // delete arc if no residual flow
      if (flow_map[arc] == 0) {
	mfGraph.erase(arc);
      }
    }

    patt_stream << read_map[source]->start() + 1 << "\t" << read_map[sink]->end();
    patt_stream << "\t" << componentID << "\t" << flownum << "\t" << path_flow;
    patt_stream << "\t" << pattern.getMethString() << std::endl;

    // recompute residual flow
    total_flow -= path_flow;
  }

  // all done
  return flownum;
}


} // namespace methylFlow
