// to be run :
//./Simulate < input.in > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/test0.tsv
#include "simulator.h"
#include <fstream>
#include <iostream>

std::ofstream patternFile;


void simulator::readData(){
    
    cin >> dnaLength >> readLength >> HapNum >> coverage >> error >> flag;
    if(flag > 0){
        int sum = 0;
        for( int i=0; i < HapNum; i++){
            freq.push_back(rand() %100);
            sum += freq[i];
        }
        for(int i =0; i<HapNum; i++){
            freq[i] = freq[i]*100/sum;
        }
        
        bool* A = new bool[dnaLength];
        for (int i=0; i<dnaLength; i++)
            A[i] = false;
        
        for(int i =0; i < flag; i++){
            int p = rand() % dnaLength;
            A[p] = true;
        }
        for (int i=0; i< dnaLength; i++) {
            if (A[i])
                pos.push_back(i);
        }
        
        
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

void simulator::buildMethylHap(int dnaLength, vector<int> pos, int HapNum, vector<MethylHap>& methylHapVec){
	
	patternFile.open("/cbcb/project-scratch/fdorri/Code/methylFlow/testing/simPattern.txt");
    char type;
    for(int i=0; i<HapNum; i++){
        vector<MethylInfo> hapInfo;
        //cout << "pos.size " << pos.size() << endl;
        for(unsigned int j=0; j<pos.size(); j++){
            if(rand() % 100 < 50){
                type = 'M';
            }
            else{
                type ='U';
            }
            hapInfo.push_back(MethylInfo(pos[j], type));
        }
        MethylHap methylHap;
        methylHap.length = dnaLength;
        methylHap.methyl = hapInfo;
        methylHapVec.push_back(methylHap);
		patternFile << methylHapVec[i].length << "\t";
		for(unsigned int j=0; j < methylHapVec[i].methyl.size(); j++){
			patternFile << methylHapVec[i].methyl[j].pos   << ":" << methylHapVec[i].methyl[j].type << ",";
		}
		patternFile << "\n" ;
    }
	patternFile.close();
	
}

void simulator::selectHP(int readLength, int dnaLength, vector<int> freq, int& hap, int& pos){
    int freqSum = 0;
    for( unsigned int i =0; i < freq.size(); i++){
        freqSum += freq[i];
        //cout << freq[i] <<"\t";
    }
    //cout << endl;
    int randHap = rand() % freqSum;
    //cout << randHap << endl;
    //cout << "freqSum" << freqSum << endl;
    pos = rand() % (dnaLength - readLength - 1);
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
			if ((rand() % 100 < error) && (methylHapVec[hap].methyl[i].type =='M')){
				read.methyl.push_back(MethylInfo(methylHapVec[hap].methyl[i].pos,'U'));
			}
			if ((rand() % 100 < error) && (methylHapVec[hap].methyl[i].type =='U')){
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
	
    cout << "read" << k <<  "\t" << read.start << "\t" << read.length << "\t" << "W" << "\t";
    if(read.methyl.size()==0)
		cout << "*" ;
    else{
        for(unsigned int i=0; i< read.methyl.size()-1; i++){
            cout << read.methyl[i].pos - read.start  << ":" << read.methyl[i].type << ",";
        }
        cout << read.methyl[read.methyl.size()-1].pos - read.start  << ":" << read.methyl[read.methyl.size()-1].type;
    }
    cout << "\t" << "MismatchData" << endl;
}


void simulator::simulate(){
	vector<MethylRead> data;
    readData();
    //cout << "read Done" << endl;
    buildMethylHap(dnaLength, pos, HapNum, methylHapVec);
    //cout << "hap built" << endl;
    for(int i=0; i < coverage*dnaLength/readLength; i++){
        int hap=0;
        int pos=0;
        MethylRead read;
        selectHP(readLength, dnaLength, freq, hap, pos);
        //cout << "HP Seleced " << i << endl;
        //cout << "hap " << hap << endl;
        //cout << "pos " << pos << endl;
        buildRead(hap, pos, readLength, error, methylHapVec, read);
		data.push_back(read);
		
      //  cout << "read Wrote" << endl;
    }
	sort(data.begin(), data.end());
    //cout << "Read built" << endl;
    //cout << data[30].methyl.size() << endl;
	for(unsigned int i=0; i <data.size(); i++){
        cout << i;
		writeMethylRead(data[i], i);
	}
}


