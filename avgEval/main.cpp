// to be run on cbcb server : qsub run.sh -t 3-5 -q xlarge -l mem=24G,walltime=24:00:00 -N Hap



#include <set>
#include <algorithm>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <lemon/list_graph.h>
//#include <lemon/lp_base.h>
#include <lemon/lp.h>
#include <lemon/maps.h>
#include <mflib/MethylRead.hpp>
//#include <lemon/bipartite_matching.h>
//#include <lemon/concepts/bpgraph.h>
//lemon::Lp lp;




using namespace std;

float var;
std::string input, output;




//Varibales for computing the avarage error matrices
std::ifstream inputFile;
std::ofstream outputFile;

float methylCallError = 0;
float abndncError = 0;
float TP, FP, FN = 0;

std::set<float> thr;
std::map<float, vector<float> > abdncErr_avg_map;
std::map<float, vector<float> > methylErr_avg_map;
std::map<float, vector<float> > TP_avg_map;
std::map<float, vector<float> > FN_avg_map;
std::map<float, vector<float> > FP_avg_map;


void readData(){
    
    float threshold, abdErr, methErr, TP1, FN1, FP1;
    
    cout << "REad Data start" << endl;

    //std::string fileName = dir;
    //fileName += str;
    std::string line;
    int chek =0;
    thr.clear();
   while (std::getline(inputFile, line)) {
       std::stringstream ss(line);
    //for(int i =0 ; i <10 ; i++){
        ss >> threshold >> abdErr >> methErr >> TP1 >> FN1 >> FP1;
        //cout << threshold << "\t" << abdErr<< "\t" << methErr << "\t"<< TP1 << "\t" << FN1<< "\t" << FP1 << endl;
        if ( thr.find(threshold) == thr.end()){
            thr.insert(threshold);
            abdncErr_avg_map[threshold].clear();
            methylErr_avg_map[threshold].clear();
            TP_avg_map[threshold].clear();
            FN_avg_map[threshold].clear();
            FP_avg_map[threshold].clear();
            
        }
       if( abs(threshold - 0.1) < 0.001 ){
       chek++;
       }

        abdncErr_avg_map[threshold].push_back(abdErr);
        methylErr_avg_map[threshold].push_back(methErr);
        TP_avg_map[threshold].push_back(TP1);
        FN_avg_map[threshold].push_back(FN1);
        FP_avg_map[threshold].push_back(FP1);
       
    }
    cout << "check " << chek << endl;
    
}

float meanFloat(vector<float> vec){
    float sum = 0;
    for (unsigned int i=0; i < vec.size() ; i++) {
        sum += vec[i];
    }
    return sum/vec.size();
}


void computeAverageErrorMatrix(float var, std::ofstream *avgFile){
    
    
    std::set<float>::iterator it;

    for (it=thr.begin(); it!=thr.end(); ++it){
        
        //cout << "it = thr.begin is " << endl;

        //cout << "threshold " << *it << endl;
        abndncError = meanFloat(abdncErr_avg_map[*it]);
        //cout << "abd size " <<abdncErr_avg_map[*it].size()<< endl;
        methylCallError = meanFloat(methylErr_avg_map[*it]);
        //cout << "meth size " <<methylErr_avg_map[*it].size()<< endl;

        TP = meanFloat(TP_avg_map[*it]);
        FP = meanFloat(FP_avg_map[*it]);
        FN = meanFloat(FN_avg_map[*it]);
        
        (*avgFile) <<  var << "\t" << *it << "\t" <<abndncError << "\t" << methylCallError << "\t" << TP << "\t" << FN  << "\t"  << FP << std::endl;
    }
    
}

int main (int argc, char* argv[]) {
    
    std::stringstream buffer;

	cout << "start main" << endl;
    
	if(argc < 4){
		cout << "Please enter your input" << endl;
		return -1;
	}
	
    std::string indirname = argv[1];
    
    buffer.str("");
    buffer << indirname << "/eval.txt";
    inputFile.open( buffer.str().c_str() );
    

    std::string outdirname = argv[2];
    
    buffer.str("");
    buffer << outdirname << "/evalAvg.txt";
    outputFile.open( buffer.str().c_str() ,std::ios_base::app);
    
    
   // inputFile.open(argv[1]);
   // outputFile.open(argv[2], ios::app );
    
    var = atof(argv[3]);
    
	
	//################     read the input .tsv data to the "line" number
	cout << "read Data " << endl;
    readData();
    cout << " Compute average " << endl;
    computeAverageErrorMatrix(var, &outputFile);
    inputFile.close();
    outputFile.close();
    
    
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
    
    return 0;
    
    
    
}


