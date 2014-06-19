// to be run on cbcb server : qsub run.sh -t 3-5 -q xlarge -l mem=24G,walltime=24:00:00 -N Hap



#include <algorithm>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include "evaluation.h"
#include <iostream>
#include <lemon/list_graph.h>
//#include <lemon/lp_base.h>
#include <lemon/lp.h>
#include <lemon/maps.h>
#include <mflib/MethylRead.hpp>
//#include <lemon/bipartite_matching.h>
//#include <lemon/concepts/bpgraph.h>
//lemon::Lp lp;


typedef lemon::ListDigraph Graph;
typedef Graph::Node RedNode;
typedef Graph::Node BlueNode;
typedef Graph::NodeIt RedNodeIt;
typedef Graph::NodeIt BlueNodeIt;
typedef Graph::Node Node;
typedef Graph::Arc Edge;
typedef Graph::ArcIt EdgeIt;
typedef Graph::InArcIt InEdgeIt;
typedef Graph::ArcMap<float> LengthMap;

using namespace methylFlow;

std::map<RedNode, MethylRead*> read_map;

using lemon::INVALID;

Graph g;
LengthMap length(g);

int chr;

std::ifstream truePatternFile, estimatedPatternFile;
std::vector<MethylRead*> trueMethylData,estimatedMethylData;

Graph::Node addNode_Read(MethylRead *read)
{
    Graph::Node n = g.addNode();
    read_map[n] = read;
    return n;
}
void readTruePattern(int start, int end){
	// "s" is for start position, "e" is for end position
	int  s, e, cid, pid;
    float abundance;
	string chr, methylString;
	int i = 0;
	while(!truePatternFile.eof()) {
		i++;
		truePatternFile >> chr >> s >> e >> cid >> pid >> abundance >> methylString;
        cout << "Ture read : "<< chr << " " << s << " " << e <<endl;
		if( s >= start && e <= end){
            MethylRead* m = new MethylRead(s, e-s+1);
            m->parseMethyl(methylString);
            
         //   m->write();
            
            trueMethylData.push_back(m);
		}
	}
    
}

void readEstimatedPattern(int start, int end){
	// "s" is for start position, "e" is for end position
	int  s, e, cid, pid;
    float abundance;
	string  chr, methylString;
    string dummyLine;
    getline(estimatedPatternFile, dummyLine);
	int i = 0;
	while(!estimatedPatternFile.eof()) {
		i++;
		estimatedPatternFile >> chr >> s >> e >> cid >> pid >> abundance >> methylString;
        
        cout << "Estimated read : "<< chr << " " << s << " " << e <<endl;
		if( s >= start &&  e <= end){
            
            MethylRead* m = new MethylRead(s, e-s+1);
            m->parseMethyl(methylString);
      //      m->write();
            estimatedMethylData.push_back(m);
		}
	}
    
}

float cost(Graph::Node u, Graph::Node v) {
    MethylRead* readU = read_map[u];
    MethylRead* readV = read_map[v];
  //  cout << "Distance start" << endl;
    int cost1 =  readU->distance(readV);
    int cost2 = readV->distance(readU);
    return cost1+cost2;
}

void buildGraph() {

    cout << "trueMethylSize = " << trueMethylData.size() << endl;
    cout << "estimatedMethylSize = " << estimatedMethylData.size() << endl;
    for (int i=0; i < trueMethylData.size(); i++) {
        cout << "Add true node " << i << endl;
        addNode_Read(trueMethylData.at(i));
    }
    
    for (int i=0; i < estimatedMethylData.size(); i++) {
        cout << "Add estimated node " << i << endl;
        addNode_Read(estimatedMethylData.at(i));
    }
    
    for (RedNodeIt u(g); u != INVALID; ++u) {
        for (RedNodeIt v(g); v != INVALID; ++v) {
            if (g.id(u) < trueMethylData.size() && g.id(v) >= trueMethylData.size()) {
                cout << "Add arc " << g.id(u) <<" " <<g.id(v) << endl;

                g.addArc(u, v);
            }
        }
    }
    
    
    
    ///////////Set Weights///////////////

    for (EdgeIt e(g); e != INVALID; ++e) {
        RedNode u = g.source(e);
        BlueNode v = g.target(e);
        length[e] = cost(u, v);
        cout << "Cost arc (" << g.id(u) <<", " <<g.id(v) << ") = " << length[e] << endl;
    }
    
}

int main (int argc, char* argv[]) {
    
    
	int start, end;
	//cout << "end main" << endl;
    
	if(argc < 5){
		cout << "Please enter your input" << endl;
		return -1;
	}
	if (argc >= 5){
		start = atoi(argv[3]);
		end = atoi(argv[4]);
		//line = atoi(argv[2]);
	}
	truePatternFile.open(argv[1]);
    estimatedPatternFile.open(argv[2]);

	
	//################     read the input .tsv data to the "line" number
	cout << "reading data " << start << endl;
	readTruePattern(start , end);
    readEstimatedPattern(start , end);

	cout << "build Graph " << start << endl;
    buildGraph();
    
 	cout << "Solve LP " << start << endl;
    
    ////////////LP//////////////////////
    lemon::Lp lp;
    lp.max();
    
    //Graph::ArcMap<lemon::Lp::Col> x;
    std::map<Edge, lemon::Lp::Row> x;
    
    std::map<RedNode, lemon::Lp::Col> alpha;
  //  std::map<Graph::BlueNode, lemon::Lp::Col> beta;
    
  
    lemon::Lp::Expr obj;
    
    for (RedNodeIt r(g); r != INVALID; ++r) {
        alpha[r] = lp.addCol();
        lp.colLowerBound(alpha[r], 0.0);
        obj += alpha[r];
    }
    
 
    for (EdgeIt e(g); e != INVALID; ++e) {
        
        x[e] = lp.addRow();
        RedNode r = g.source(e);
        BlueNode b = g.target(e);
        
        x[e] = lp.addRow(alpha[r] + alpha[b] - length[e] <= 0);
    }
    
    
    
    lp.obj(obj);
    lp.solve();
    
    float v = lp.primal();
    cout << "dual : " << v << endl;

    for (EdgeIt e(g); e != INVALID; ++e) {
        
        lemon::Lp::Row row = x[e];
        float dual = lp.dual(row);
        
        cout << " X(" << g.id(g.source(e)) <<", " << g.id(g.target(e)) << ") = " << dual<< endl;
      //  cout << " W(" << g.id(g.source(e)) <<", " << g.id(g.target(e)) << ") = " << length[e]<< endl;
    }
    
    
    std::cout << std::endl;
}


