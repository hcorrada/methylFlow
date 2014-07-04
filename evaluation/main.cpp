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
std::map<RedNode, int> abundance_map;

std::ofstream evalFile;



using lemon::INVALID;

Graph g;
LengthMap length(g);

int chr,var;

std::ifstream truePatternFile, estimatedPatternFile, mFile, wFile;
std::ofstream weightFile, matchFile;
std::vector<MethylRead*> trueMethylData,estimatedMethylData;
std::vector<float> trueAbundanceData, estimatedAbundanceData;


/////Variables for reading MatchMatrix
int truePatternNum, estimatedPatternNum;
vector<int> abdncTrue, abdncEstimated, idTrue, idEstimated;
std::map<int, float> abdnc_map;
std::map<int, int> matchTrue_map;
std::map<int, int> matchEstimated_map;
std::map<int, float> weight_map;



Graph::Node addNode_Read(MethylRead *read, float abndnc)
{
    Graph::Node n = g.addNode();
    read_map[n] = read;
    abundance_map[n] = abndnc;
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
            
            m->write();
            
            trueMethylData.push_back(m);
            trueAbundanceData.push_back(abundance);
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
		if( s >= start -1  &&  e <= end){
            
            MethylRead* m = new MethylRead(s, e-s+1);
            m->parseMethyl(methylString);
            m->write();
            estimatedMethylData.push_back(m);
            estimatedAbundanceData.push_back(abundance);
		}
	}
    float sum =0;
    for(unsigned int j= 0; j<estimatedAbundanceData.size(); j++){
        sum += estimatedAbundanceData.at(j);
    }
    for(unsigned int j= 0; j<estimatedAbundanceData.size(); j++){
        estimatedAbundanceData.at(j) = estimatedAbundanceData.at(j)*100/sum ;
    }
    
}

float cost(Graph::Node u, Graph::Node v) {
    MethylRead* readU = read_map[u];
    MethylRead* readV = read_map[v];
    int common = 0 ;
    //  cout << "Distance start" << endl;
    int match =  readU->distance(readV, common);
    cout << "match " << match <<  endl;
    
    cout << "common " << common <<  endl;
    
    int mismatch = readU->cpgOffset.size() + readV->cpgOffset.size() - match - common ;
    int totalCpG = readU->cpgOffset.size() + readV->cpgOffset.size() - common;
    cout << "cost " << (float(mismatch) / totalCpG) <<  endl;
    
    return (float(mismatch) / totalCpG) ;
}

void buildGraph() {
    
    cout << "trueMethylSize = " << trueMethylData.size() << endl;
    cout << "estimatedMethylSize = " << estimatedMethylData.size() << endl;
    for (int i=0; i < trueMethylData.size(); i++) {
        cout << "Add true node " << i << endl;
        addNode_Read(trueMethylData.at(i), trueAbundanceData.at(i));
    }
    
    for (int i=0; i < estimatedMethylData.size(); i++) {
        cout << "Add estimated node " << i << endl;
        addNode_Read(estimatedMethylData.at(i), estimatedAbundanceData.at(i));
    }
    
    for (RedNodeIt u(g); u != INVALID; ++u) {
        for (RedNodeIt v(g); v != INVALID; ++v) {
            if (g.id(u) < trueMethylData.size() && g.id(v) >= trueMethylData.size()) {
                cout << "Add arc " << g.id(u) << " " <<g.id(v) << endl;
                
                g.addArc(u, v);
            }
        }
    }
    
    
    
    ///////////Set Weights///////////////
    
    for (EdgeIt e(g); e != INVALID; ++e) {
        RedNode u = g.source(e);
        BlueNode v = g.target(e);
        length[e] = cost(u, v);
        
        cout << "Cost arc (" << g.id(u) <<", " << g.id(v) << ") = " << length[e] << endl;
    }
}

void writeMatchMatrix() {
    weightFile.open("/cbcb/project-scratch/fdorri/Code/methylFlow/testing/weight.txt");
    matchFile.open("/cbcb/project-scratch/fdorri/Code/methylFlow/testing/match.txt");
    
    weightFile << trueMethylData.size() << "\t" << estimatedMethylData.size() << endl;
    matchFile << trueMethylData.size() << "\t" << estimatedMethylData.size() << endl;
    
    for (RedNodeIt u(g); u != INVALID; ++u) {
        if (g.id(u) < trueMethylData.size()) {
            weightFile << g.id(u) << "\t" << abundance_map[u] << endl;;
            matchFile << g.id(u) << "\t" << abundance_map[u] << endl;
        }
    }
    
    
    for (RedNodeIt u(g); u != INVALID; ++u) {
        if (g.id(u) >= trueMethylData.size()) {
            weightFile << g.id(u) << "\t" << abundance_map[u] << endl;;
            matchFile << g.id(u) << "\t" << abundance_map[u] << endl;
        }
    }
    
    
    
    
    int argmin;
    for (RedNodeIt u(g); u != INVALID; ++u) {
        double min = 100000;
        
        if (g.id(u) < trueMethylData.size()) {
            
            for (RedNodeIt v(g); v != INVALID; ++v) {
                if (g.id(v) >= trueMethylData.size()) {
                    float weight = cost(u, v);
                    weightFile << var << "\t" << g.id(u) << "\t" << g.id(v) << "\t" << weight << endl;
                    
                    if(min > weight ){
                        min = cost(u,v);
                        argmin = g.id(v);
                    }
                }
            }
            
            matchFile << var << "\t" << g.id(u) << "\t" << argmin  << "\t" << min << endl;
        }
    }
    matchFile.close();
    weightFile.close();
    
}


void readMatchMatrix() {
    //reading the match file for pruning the matched edges
    int id, var, matchId;
    float abdnc, weight;
    
    mFile.open("/cbcb/project-scratch/fdorri/Code/methylFlow/testing/match.txt");
    mFile >> truePatternNum  >> estimatedPatternNum ;
    for (int i = 0; i < truePatternNum ; i++){
        mFile >> id >> abdnc;
        idTrue.push_back(id);
        abdnc_map[id] = abdnc;
    }
    for (int i = 0; i < estimatedPatternNum ; i++){
        mFile >> id >> abdnc;
        idEstimated.push_back(id);
        abdnc_map[id] = abdnc;
    }
    
    for (int i=0; i< truePatternNum ; i++ ){
        mFile >> var >> id >> matchId >> weight ;
        matchTrue_map[id] = matchId;
        weight_map[id] = weight;
        matchEstimated_map[matchId] = id;
    }
    
}

void computeErrorMatrix(double threshold) {
    float methylCallError = 0;
    float abndncError = 0;
    
    int match = 0;
    
    for (int i = 0 ; i < truePatternNum ; i ++){
        if(weight_map[i] < threshold) {
            abndncError += pow(double(abdnc_map[i] - abdnc_map[matchTrue_map[i]]),2) / 10000;
            methylCallError += weight_map[i];
            match++;
        }
        else{
            abndncError += pow(double(abdnc_map[i]),2) /10000;
            abndncError += pow(double(abdnc_map[matchTrue_map[i]]),2) /10000;
            //methylCallError += weight_map[i];
        }
    }
    
    for(int i = truePatternNum; i < truePatternNum + estimatedPatternNum ; i++){
        if (matchEstimated_map.find(i) == matchEstimated_map.end()) {
            abndncError += pow(double(abdnc_map[i]),2) /10000;
        }
        
    }
    
    
    int TP = match;
    int FN = truePatternNum - match;
    int FP = estimatedPatternNum - match;
    
    methylCallError = methylCallError/match;
    
    //evalFile.open("/cbcb/project-scratch/fdorri/Code/methylFlow/testing/eval.txt", std::ios_base::app);
    
    //evalFile << "abndncError " << "\t"  << "methylCall Error " << "\t" << "TP "<< "\t" << "FN " << "\t" << "FP " << endl;
    evalFile << threshold << "\t" <<abndncError << "\t" << methylCallError << "\t" << TP << "\t" << FN  << "\t"  << FP << std::endl;
    
    //std:cout << abndncError << "\t" << methylCallError << "\t" << TP << "\t" << FN  << "\t"  << FP << std::endl;
    

}

int main (int argc, char* argv[]) {
    
    
	int start, end;
	cout << "start main" << endl;
    
	if(argc < 7){
		cout << "Please enter your input" << endl;
		return -1;
	}
	if (argc >= 7){
		start = atoi(argv[4]);
		end = atoi(argv[5]);
		var = atoi(argv[6]);
        
	}
	truePatternFile.open(argv[1]);
    estimatedPatternFile.open(argv[2]);
    evalFile.open(argv[3],std::ios_base::app);
    
	
	//################     read the input .tsv data to the "line" number
	cout << "reading data " << start << endl;
	readTruePattern(start , end);
    readEstimatedPattern(start , end);
    
	cout << "build Graph " << start << endl;
    buildGraph();
    
    
    //########### write the weight of edges and corresponding matches for every true pattern ########
    /// the ffirst line is #truePattern , #estimated Pattern
    /// then each line is true Pattern id , abundance of pattern
    /// then each line is estimated pattern id , abundance of pattern
    // then the weight information of edges and matching information is seen in the rest of file
    
    
    writeMatchMatrix();
    

    readMatchMatrix();
    
    
    //// computing Error Metric ////////
    for (double thresh=0.05; thresh<0.30; thresh +=0.05)
        computeErrorMatrix(thresh);
    

//}

/*
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
 //lp.colLowerBound(alpha[r], 0.0);
 obj += alpha[r];
 }
 
 
 
 for (EdgeIt e(g); e != INVALID; ++e) {
 
 //   x[e] = lp.addRow();
 RedNode r = g.source(e);
 BlueNode b = g.target(e);
 
 x[e] = lp.addRow(alpha[r] + alpha[b] - length[e] <= 0);
 }
 
 
 
 lp.obj(obj);
 lp.solve();
 
 float v = lp.primal();
 //cout << "dual : " << v << endl;
 
 float abndncError = 0;
 float methylCallError = 0;
 
 int match = 0;
 for (EdgeIt e(g); e != INVALID; ++e) {
 
 lemon::Lp::Row row = x[e];
 float dual = lp.dual(row);
 
 cout << " X(" << g.id(g.source(e)) <<", " << g.id(g.target(e)) << ") = " << dual<< endl;
 cout << " W(" << g.id(g.source(e)) <<", " << g.id(g.target(e)) << ") = " << length[e]<< endl;
 
 if (dual){
 abndncError += pow(double(abundance_map[g.target(e)] - abundance_map[g.source(e)]),2) / 10000;
 methylCallError += length[e];
 match++;
 }
 else{
 abndncError += pow(double(abundance_map[g.target(e)]),2) /10000 + pow(double(abundance_map[g.source(e)]),2) /10000;
 methylCallError +=1;
 }
 }
 int TP = match;
 int FN = trueMethylData.size() - match;
 int FP = estimatedMethylData.size() - match;
 
 methylCallError = methylCallError/(match + FN +FP);
 
 //evalFile.open("/cbcb/project-scratch/fdorri/Code/methylFlow/testing/eval.txt", std::ios_base::app);
 
 //evalFile << "abndncError " << "\t"  << "methylCall Error " << "\t" << "TP "<< "\t" << "FN " << "\t" << "FP " << endl;
 evalFile << abndncError << "\t" << methylCallError << "\t" << TP << "\t" << FN  << "\t"  << FP << std::endl;
 
 //std:cout << abndncError << "\t" << methylCallError << "\t" << TP << "\t" << FN  << "\t"  << FP << std::endl;
 */





}


