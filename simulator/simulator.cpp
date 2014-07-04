// to be run :
//./Simulate < input.in > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/test0.tsv
// dataFlag = 0  >>> read CpG sites from file
// dataFlag > 0 >>>> dataFlag equals the number of cpg sites
// dataFlag < 0 >>> read the data from rest of the file

//freqFlag = 0 >>> randomly choose the frequency of each pattern
//freqFlag = 1 >>> read the frequency of patterns from rest of the file(second line)
#include "simulator.h"
#include <math.h>
#include <fstream>
#include <iostream>

std::ofstream patternFile;
std::ifstream methylPosFile;

void simulator::readData(){
    
    cin >> chr >> startDNA >> dnaLength >> readLength >> HapNum >> freqFlag >> coverage >> error >> dataFlag >> corrDist;
    if(dataFlag > 0){
        if (freqFlag == 1){
            int x;
            for( int i=0; i < HapNum; i++){
                cin >> x;
                cerr << "freq " << x << endl;
                freq.push_back(x);
            }
        }
        else{
            int sum = 0;
            for( int i=0; i < HapNum; i++){
                freq.push_back(rand() %100);
                sum += freq[i];
            }
            for(int i =0; i<HapNum; i++){
                freq[i] = freq[i]*100/sum;
                //cout << freq[i] << endl;
            }
        }
        
        bool* A = new bool[dnaLength];
        for (int i=0; i<dnaLength; i++)
            A[i] = false;
        
        for(int i =0; i < dataFlag; i++){
            int p = rand() % dnaLength;
            A[p] = true;
        }
        for (int i=0; i< dnaLength; i++) {
            if (A[i])
                pos.push_back(i + startDNA);
        }
        
        
    }
    else if( dataFlag == 0){
        if (freqFlag == 1){
            int x;
            for( int i=0; i < HapNum; i++){
                cin >> x;
                cerr << "freq " << x << endl;
                freq.push_back(x);
            }
        }
        else{
            int sum = 0;
            for( int i=0; i < HapNum; i++){
                freq.push_back(rand() %100);
                sum += freq[i];
            }
            for(int i =0; i<HapNum; i++){
                freq[i] = freq[i]*100/sum;
                //cout << freq[i] << endl;
            }
        }
        
        methylPosFile.open("/cbcb/project-scratch/fdorri/Code/Methylation/methylFlow/1_CpGInfo.txt");
        vector<int> tempPos;
        int s, p;
        char t;
        string dummyLine;
        getline(methylPosFile, dummyLine);
        int i = 0;
        //while(!methylPosFile.eof()) {
        while( i < 100){
            i++;
            methylPosFile >> s >> p >> t;
            //dnaLength = max(dnaLength, s + p );
            //dnaMin = min(dnaMin, s);
            //cerr << "s " << s << " ,p " << p << endl;
            //readLength = max(readLength, p);
            if (s >= startDNA && s+p <= startDNA + dnaLength ){
                tempPos.push_back(s + p);
            }
        }
        sort(tempPos.begin(), tempPos.end());
        pos.push_back(tempPos[0]);
        for(int i=1; i < tempPos.size(); i++){
            if(tempPos[i] != tempPos[i-1]){
                pos.push_back(tempPos[i]);
                cerr << tempPos[i] << endl;
            }
        }
        cerr << "posSize " << pos.size()<< endl;
        cerr << "tempposSize " << tempPos.size()<< endl;
        

        
        
        cerr << dnaLength << endl;
        cerr << readLength << endl;
        
        methylPosFile.close();
    }
    else{
        for(int i=0; i < HapNum; i++){
            int f;
            cin >> f;
            freq.push_back(f);
        }
        int posNum;
        cin >> posNum;
        for(int i=0; i < posNum; i++){
            int p;
            cin >> p;
            pos.push_back(p);
        }
    }
    
    
}
double simulator::sigmoid(double x, int corrDist){
    return 1.0/ (1.0 + exp(-10*(x-corrDist)));
}
int simulator::computeMethylProbability(int corrDist, int dist){
	double r = double(rand())/ RAND_MAX;
	//if( dist < corrDist){
    if ( r > sigmoid(dist, corrDist)){
        return 0;
    }
    else if(rand() % 100 < 50)
        return 0;
    else
        return 1;
    
}


void simulator::buildMethylHap(int dnaLength, vector<int> pos, int HapNum, vector<int> freq, vector<MethylHap>& methylHapVec, int dataFlag){
	
	patternFile.open("/cbcb/project-scratch/fdorri/Code/methylFlow/testing/simPattern.txt");
    char type;
    vector<MethylInfo> hapInfo;
    for(int i=0; i<HapNum; i++){
        hapInfo.clear();
        if(rand() % 100 < 50){
            type = 'M';
        }
        else{
            type ='U';
        }
        hapInfo.push_back(MethylInfo(pos[0], type));
        cerr << "posSize 1 " << pos.size() << endl;
        for( unsigned int j=1; j < pos.size() ; j++){
            int dist = pos[j] - pos[j-1];
            if (computeMethylProbability(corrDist, dist) == 0)
                hapInfo.push_back(MethylInfo(pos[j], type));
            else{
                if (type == 'M'){
                    hapInfo.push_back(MethylInfo(pos[j], 'U'));
                    type = 'U';
                }
                else if (type == 'U'){
                    hapInfo.push_back(MethylInfo(pos[j], 'M'));
                    type = 'M';
                }
            }
            
        }
        MethylHap methylHap;
        methylHap.length = dnaLength;
        methylHap.methyl = hapInfo;
        methylHapVec.push_back(methylHap);
        //patternFile << methylHapVec[i].length << "\t";
        patternFile << chr << "\t" << startDNA << "\t" << dnaLength << "\t" <<"1"  << "\t" << i << "\t" << freq[i] << "\t";
        if( i < HapNum -1){
        //cout << "dnaLength" << dnaLength << endl;
        for(unsigned int j=0; j < methylHapVec[i].methyl.size(); j++){
            patternFile << (methylHapVec[i].methyl[j].pos - startDNA) << ":" << methylHapVec[i].methyl[j].type << ",";
        }
        patternFile << endl ;
        }
        if(i == HapNum -1 ){
            for(unsigned int j=0; j < methylHapVec[i].methyl.size(); j++){
                patternFile << (methylHapVec[i].methyl[j].pos - startDNA) << ":" << methylHapVec[i].methyl[j].type << ",";
            }
        }
        
    }
    patternFile.close();
}

void simulator::selectHP(int readLength, int dnaLength, vector<int> freq, int& hap, int& pos){
    int freqSum = 0;
    //cerr << "freqSize " << freq.size() <<"\t";
    for( unsigned int i =0; i < freq.size(); i++){
        freqSum += freq[i];
        // cerr << "freqHap" << freq[i] <<"\t";
    }
    //cout << endl;
    int randHap = rand() % freqSum;
    //cout << randHap << endl;
    //cout << "freqSum" << freqSum << endl;
    pos = rand() % (dnaLength - readLength - 1) + startDNA;
    int sum = 0;
    for(unsigned int i=0; i < freq.size(); i ++){
        if( randHap < freq[i] + sum){
            hap = i;
            break;
        }
        else{
            sum  += freq[i];
        }
    }
}

void simulator::buildRead(int hap, int pos,int readLength, int error, vector<MethylHap> methylHapVec, MethylRead& read){
    
    read.start = pos;
    read.length = readLength;
    //cout << "hap: " << hap << endl;
    //cout << "size " << methylHapVec[hap].methyl.size() << endl;
    
    for(unsigned int i=0; i< methylHapVec[hap].methyl.size(); i++){
        if( (methylHapVec[hap].methyl[i].pos >= pos) &&  (methylHapVec[hap].methyl[i].pos < pos + readLength)){
            int randErr = rand() % 100 ;
			if ((randErr < error) && (methylHapVec[hap].methyl[i].type =='M')){
				read.methyl.push_back(MethylInfo(methylHapVec[hap].methyl[i].pos,'U'));
			}
			else if ((randErr < error) && (methylHapVec[hap].methyl[i].type =='U')){
				read.methyl.push_back(MethylInfo(methylHapVec[hap].methyl[i].pos,'M'));
			}
			else{
				read.methyl.push_back(MethylInfo(methylHapVec[hap].methyl[i].pos,methylHapVec[hap].methyl[i].type));
			}
        }
        if(methylHapVec[hap].methyl[i].pos > pos + readLength)
            break;
        
    }
}

void simulator::writeMethylRead(MethylRead read, int k){
	
    cout << "read" << k <<  "\t" << read.start  << "\t" << read.length << "\t" << "W" << "\t";
    if(read.methyl.size() == 0)
        cout << "*";

    if(read.methyl.size() > 0){
        //cout << "read" << k <<  "\t" << read.start + 1 << "\t" << read.length << "\t" << "W" << "\t";
        for(unsigned int i=0; i< read.methyl.size()-1; i++){
            cout << read.methyl[i].pos - read.start  << ":" << read.methyl[i].type << ",";
        }
        cout << read.methyl[read.methyl.size()-1].pos - read.start << ":" << read.methyl[read.methyl.size()-1].type;
        //cout << "\t" << "MismatchData" << endl;
    }
    cout << "\t" << "MismatchData" << endl;

}


void simulator::simulate(){
	vector<MethylRead> data;
    readData();
    cerr << "DNA Lenght" << "\t" << dnaLength<< endl;
    buildMethylHap(dnaLength, pos, HapNum, freq, methylHapVec, dataFlag);
    cerr << "hap built" << endl;
    cerr << "max_i" << coverage*dnaLength/readLength << endl;
    
    for(int i=0; i < coverage*dnaLength/readLength; i++){
        int hap=0;
        int pos=0;
        MethylRead read;
        //cerr << "readLength " << readLength << endl;
        
        selectHP(readLength, dnaLength, freq, hap, pos);
        //cerr << "HP Seleced " << i << endl;
        //cerr << "hap " << hap << endl;
        //cerr << "pos " << pos << endl;
        //cerr << "read Length" <<  readLength << endl;
        buildRead(hap, pos, readLength, error, methylHapVec, read);
		data.push_back(read);
		
        //cerr << "read Wrote" << endl;
    }
	sort(data.begin(), data.end());
    cerr << "Read built" << endl;
    cerr << data[30].methyl.size() << endl;
	for(unsigned int i=0; i <data.size(); i++){
        //cout << i;
		writeMethylRead(data[i], i);
	}
}


