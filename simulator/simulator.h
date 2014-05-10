#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <string>
#include <sstream>

using namespace std;

struct MethylInfo{
	int pos;
	char type;
	MethylInfo(int p, char t): pos(p), type(t){}
};

struct MethylHap{
    int length;
    vector<MethylInfo> methyl;
};

struct MethylRead{
    int start, length;
    vector<MethylInfo> methyl;
	bool operator < (const MethylRead& m) const
	{ return (this->start < m.start || (this->start == m.start && this->length > m.length));}
};



class simulator {
public:
    
    vector<int> freq;
    vector<int> pos;
    int dnaLength, readLength, HapNum, coverage, error, flag, corrDist;
	vector<MethylHap> methylHapVec;
	
	simulator(){};
    
    ~simulator(){};

    void readData();
    double sigmoid(double x, int corrDist);
    int computeMethylProbability(int corrDist, int dist);
	void buildMethylHap(int length, vector<int> pos, int HapNum, vector<MethylHap>& methylHapVec, int flag);
//    void buildMethylHap(int length, vector<int> pos, int HapNum);

	void selectHP(int readLength, int dnaLength, vector<int> freq, int& hap, int& pos);
	void buildRead(int hap, int pos,int readLength, int error, vector<MethylHap> methylHapVec, MethylRead& read);
    
    void writeMethylRead(MethylRead read, int k);
    void simulate();
	
	
};
