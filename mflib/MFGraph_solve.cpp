#include <iostream>
#include <queue>

#include <lemon/path.h>
#include <lemon/dijkstra.h>

#include "MFGraph.hpp"
#include "MFRegionSolver.hpp"
#include "MFCpgSolver.hpp"

namespace methylFlow {
    void MFGraph::add_terminals()
    {
        int rightMostStart = -1;
        int rightMostEnd = -1;
        
        
        // find childless nodes
        for (ListDigraph::NodeIt n(mfGraph); n != INVALID; ++n) {
            if (countOutArcs(mfGraph, n) == 0) {
                childless[n] = true;
                if (read_map[n] && read_map[n]->start() > rightMostStart) {
                    rightMostStart = read_map[n]->start();
                    rightMostEnd = rightMostStart + read_map[n]->length();
                    //          std::cout << "rightMostEnd" << rightMostEnd << std::endl;
                    
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
       // MethylRead *sink_read = new MethylRead(rightMostStart , rightMostEnd- rightMostStart+1);
        MethylRead *sink_read = new MethylRead(rightMostEnd + 1 , 1);

        sink = addNode("t", 0, sink_read);
        fake[sink] = true;
        for (ListDigraph::NodeIt n(mfGraph); n != INVALID; ++n) {
            if (childless[n]) {
                
                addArc(n, sink, sink_read->start() - read_map[n]->start());
                //addArc(n, sink, 1);
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
            //      std::cout << read->start() << std::endl;
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
            // keep going if node is invalid
            if (!mfGraph.valid(n)) continue;
            
#ifndef NDEBUG
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
#ifndef NDEBUG
        std::cout << "median coverage: " << median_coverage << std::endl;
#endif
        
        for (ListDigraph::NodeIt n(mfGraph); n != INVALID; ++n) {
            normalized_coverage_map[n] *= median_coverage;
        }
        is_normalized = true;
    }
    
    void MFGraph::regularize()
    {
//        
//        int rightMostEnd = -1;
//        // find childless nodes
//        for (ListDigraph::NodeIt n(mfGraph); n != INVALID; ++n) {
//            if (countOutArcs(mfGraph, n) == 0) {
//                childless[n] = true;
//                if (read_map[n] && read_map[n]->start() > rightMostStart) {
//                    rightMostEnd = rightMostStart + read_map[n]->length();
//                    
//                }
//            }
//        }

        // add the lambda edges
        int lamcnt = 0;
        for (ListDigraph::InArcIt arc(mfGraph, sink); arc != INVALID; ++arc, ++lamcnt) {
            std::stringstream nodename;
            nodename << "lambda_" << lamcnt;
            ListDigraph::Node newNode = addNode(nodename.str(), 0);
            ListDigraph::Node node = mfGraph.source(arc);
            
            MethylRead read = MethylRead(*read_map[node]);
            int start = read_map[node]->start();
            int end = read_map[node]->end();
            
            //int newEnd = max(start, min(end, rightMostEnd - rLen));
            
            //ListDigraph::Arc newArc = addArc(node, newNode, 1);
            ListDigraph::Arc newArc = addArc(node, newNode, end - start + 1);
            
            mfGraph.changeSource(arc, newNode);
            
            childless[node] = false;
            childless[newNode] = true;
        }
    }
    
    int MFGraph::solve(const float lambda, const float length_mult, const float epsilon, const bool verbose, const bool pctselect )
    {
        int res;
        
        if (verbose) {
            std::cout << "[methylFlow] Extending graph with regularization nodes" << std::endl;
        }
        regularize();
        
        for (ListDigraph::ArcIt arc(mfGraph); arc != INVALID; ++arc) {
#ifndef NDEBUG
            std::cout << "Extracting flow of arc: " << std::endl;
            std::cout << nodeName_map[ mfGraph.source(arc)];
            std::cout << " -> ";
            std::cout << nodeName_map[ mfGraph.target(arc)] << "\t";
            std::cout << effectiveLength_map[arc] << std::endl;
#endif
        }
        MFSolver *solver;
        if (pctselect) {
            solver = new MFCpgSolver(this, length_mult);
        } else {
            solver =new MFRegionSolver(this);
        }
        
        if (verbose) {
            std::cout << "[methylFlow] Solving optimization problem" << std::endl;
        }
        res = solver->solve(lambda, length_mult, epsilon, verbose);
        
        if (res) {
            std::cerr << "[methylFlow] Error solving for scale =" << length_mult << std::endl;
            return res;
        }
        
        solver->extract_flows();
        return 0;
    }
    
    
    int MFGraph::decompose(const int componentID, std::ostream & patt_stream, std::string chr)
    {
        IdMap<ListDigraph, ListDigraph::Node> idmap(mfGraph);
        
        // compute total flow
        float total_flow = this->total_flow();
#ifndef NDEBUG
        std::cout << "total flow:: " << total_flow << std::endl;
#endif
        
        int flownum = 0;
        
        // iterate while residual flow
        while (total_flow > 0.00001) {
            
            // compute residual flow for each arc
            ListDigraph::ArcMap<float> residual_flow(mfGraph);
            for (ListDigraph::ArcIt arc(mfGraph); arc != INVALID; ++arc) {
                residual_flow[arc] = (total_flow - flow_map[arc]);
                if(flow_map[arc] != 0){
// #ifndef NDEBUG
//                     std::cout << "source: " << mfGraph.id(mfGraph.source(arc))<< ", target: " << mfGraph.id(mfGraph.target(arc)) << ", flow: " << flow_map[arc] <<std::endl;
// #endif
                }
                
            }
#ifndef NDEBUG
            std::cout << "run the min-max dijkstra algorithm " << std::endl;
#endif
            // run the min-max dijkstra algorithm
            ListDigraph::NodeMap<float> dist(mfGraph);
            Dijkstra<ListDigraph, ListDigraph::ArcMap<float> >
            ::SetOperationTraits<DijkstraMinMaxOperationTraits<float> >
            ::Create dijkstra(mfGraph, residual_flow);
            dijkstra.distMap(dist);
            dijkstra.run(source, sink);
#ifndef NDEBUG
            std::cout << "get the resulting path and it's flow " << std::endl;
#endif
            // get the resulting path and it's flow
            Path<ListDigraph> shortestPath = dijkstra.path(sink);
            float path_flow = total_flow - dijkstra.dist(sink);

            // break out if this is not a valid path
            if (path_flow == 0) {
#ifndef NDEBUG
              std::cout << "Found invalid path" << std::endl;
#endif
              break;
            }

            // construct a meth fragment from path here
            // and remove path flow from each arc in path
#ifndef NDEBUG
            std::cout << "dijkstra.dist: " << dijkstra.dist(sink) << ", path_flow:" << path_flow << ", total flow: " << total_flow << " length: " << shortestPath.length() << std::endl;
#endif

            // it's a valid pattern, so increase number of paths found
            flownum++;

            std::stringstream region_list;
            MethylRead* pattern = NULL;
            int start, end;
            //MethylRead pattern = MethylRead(*read_map[source]);
            //int start = read_map[source]->start();
            //int end = read_map[sink]->end();
            
            for(Path<ListDigraph>::ArcIt arc(shortestPath); arc != INVALID; ++arc) {
                ListDigraph::Node s = mfGraph.source(arc);
                ListDigraph::Node t = mfGraph.target(arc);
                if (s == source) {
                    pattern = new MethylRead(*read_map[t]);
                    start = read_map[t]->start();
                    end = read_map[t]->end();
                    // end = read_map[sink]->end();

#ifndef NDEBUG
                    std::cout << "Original pattern = " << pattern->getMethString() << std::endl;
#endif
                    break;
                }
                
            }

#ifndef NDEBUG
            if (!pattern) {
              std::cout << "We should not hit this" << std::endl;
              for (Path<ListDigraph>::ArcIt arc(shortestPath); arc != INVALID; ++arc) {
                ListDigraph::Node s = mfGraph.source(arc);
                ListDigraph::Node t = mfGraph.target(arc);
                std::cout << mfGraph.id(s) << " -> " << mfGraph.id(t) << std::endl;
              }
            }
#endif
             
            for(Path<ListDigraph>::ArcIt arc(shortestPath); arc != INVALID; ++arc) {
// #ifndef NDEBUG
//                 std::cout << " After finding a path, " << "source: " << mfGraph.id(mfGraph.source(arc))<< ", target: " << mfGraph.id(mfGraph.target(arc)) << ", flow: " << flow_map[arc] << " ,arc - pathFlow: " << flow_map[arc] - path_flow << std::endl;
// #endif
                
                ListDigraph::Node s = mfGraph.source(arc);
                ListDigraph::Node t = mfGraph.target(arc);
                
                MethylRead *read = read_map[t];

                // don't print the source node or nodes connected to sink
                if (s != source && t != get_sink()) {
                    region_list << idmap[s];
                    if (!childless[t]) {
                        region_list << ",";
                    }
                }
                
                flow_map[arc] -= path_flow;
                
                // delete arc if no residual flow
                if (flow_map[arc] < 1e-6) {
                    mfGraph.erase(arc);
                }
                
                if (s == source) {
                    continue;
                }
                
                
                if (!read  || t == get_sink()) {
                    continue;
                }
                
             
                if (s != source && t != get_sink()) {
                    pattern->merge(read);
                    end = read->end();
                    //std::cout << "new pattern = " << pattern->getMethString() << std::endl;
                }
                
              
              //  if (childless[t]) {
              //      end = read->end();
              //  }
                
                
            }

#ifndef NDEBUG
            std::cout << "new pattern = " << pattern->getMethString() << std::endl;
#endif

            //Note we add source one nucleotide before every read
            patt_stream << chr << "\t" << start << "\t" << end;
            patt_stream << "\t" << componentID << "\t" << flownum << "\t" << path_flow;
            patt_stream << "\t" << pattern->getMethString() << "\t" << region_list.str() << std::endl;

            delete pattern;
            
            // recompute residual flow
            total_flow -= path_flow;
        }
        
        // all done
        return flownum;
    }
    
    
} // namespace methylFlow
