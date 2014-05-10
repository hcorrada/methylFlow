// to be run :
//./Simulate < input.in > /cbcb/project-scratch/fdorri/Code/methylFlow/testing/test0.tsv
// flag = 0  >>> read CpG sites from file
// flag > 0 >>>> flag equals the number of cpg sites
// flag < 0 >>> read the data from rest of the file
#include "simulator.h"
#include <math.h>
#include <fstream>
#include <iostream>

std::ofstream patternFile;
std::ifstream methylPosFile;

void simulator::readData(){
    
    cin >> dnaLength >> readLength >> HapNum >> coverage >> error >> flag >> corrDist;
    if(flag > 0){
        int sum = 0;
        for( int i=0; i < HapNum; i++){
            freq.push_back(rand() %100);
            sum += freq[i];
        }
        for(int i =0; i<HapNum; i++){
            freq[i] = freq[i]*100/sum;
            //cout << freq[i] << endl;
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
    else if( flag == 0){
        int sum = 0;
        for( int i=0; i < HapNum; i++){
            freq.push_back(rand() %100);
            sum += freq[i];
        }
        for(int i =0; i<HapNum; i++){
            freq[i] = freq[i]*100/sum;
            //cout << freq[i] << endl;
        }
        

        methylPosFile.open("/cbcb/project-scratch/fdorri/Code/Methylation/methylFlow/1_CpGInfo.txt");
        vector<int> posVec;
        int s, p;
        char t;
        string dummyLine;
        getline(methylPosFile, dummyLine);
        int i = 0;
        //while(!methylPosFile.eof()) {
        while( i < 100){
            i++;
            methylPosFile >> s >> p >> t;
            dnaLength = max(dnaLength, s + p );
            //readLength = max(readLength, p);
            pos.push_back(s + p);
        }
        //cout << dnaLength << endl;
        //cout << readLength << endl;

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
    return 1.0/ (1.0 + exp(-(x-corrDist)));
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


void simulator::buildMethylHap(int dnaLength, vector<int> pos, int HapNum, vector<MethylHap>& methylHapVec, int flag){
	
	patternFile.open("simPattern.txt");
    char type;
    vector<MethylInfo> hapInfo;
    if(flag >0 ){
        for(int i=0; i<HapNum; i++){
            //cout << "pos.size " << pos.size() << endl;
            //hapInfo.clear();
            //for(unsigned int j=0; j<pos.size(); j++){
              //  if(rand() % 100 < 50){
                //    type = 'M';
                //}
                //else{
                  //  type ='U';
                //}
                //hapInfo.push_back(MethylInfo(pos[j], type));
            //}
            hapInfo.clear();
            if(rand() % 100 < 50){
                type = 'M';
            }
            else{
                type ='U';
            }
            hapInfo.push_back(MethylInfo(pos[0], type));
            for( unsigned int i=1; i < pos.size() ; i++){
                int dist = pos[i] - pos[i-1];
                if (computeMethylProbability(corrDist, dist) == 0)
                hapInfo.push_back(MethylInfo(pos[i], type));
                else{
                    if (type == 'M')
                    hapInfo.push_back(MethylInfo(pos[i], 'U'));
                    if (type == 'U')
                    hapInfo.push_back(MethylInfo(pos[i], 'M'));
                }
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
    
    else if(flag == 0){
        for(int i=0; i<HapNum; i++){
            hapInfo.clear();
            if(rand() % 100 < 50){
                type = 'M';
            }
            else{
                type ='U';
            }
            hapInfo.push_back(MethylInfo(pos[0], type));
            for( unsigned int i=1; i < pos.size() ; i++){
                int dist = pos[i] - pos[i-1];
                if (computeMethylProbability(corrDist, dist) == 0)
                    hapInfo.push_back(MethylInfo(pos[i], type));
                else{
                    if (type == 'M')
                        hapInfo.push_back(MethylInfo(pos[i], 'U'));
                    if (type == 'U')
                        hapInfo.push_back(MethylInfo(pos[i], 'M'));
                }
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
	
    cout << "read" << k <<  "\t" << read.start + 1 << "\t" << read.length << "\t" << "W" << "\t";
    if(read.methyl.size()==0)
		cout << "*" ;
    else{
        for(unsigned int i=0; i< read.methyl.size()-1; i++){
            cout << read.methyl[i].pos - read.start  + 1 << ":" << read.methyl[i].type << ",";
        }
        cout << read.methyl[read.methyl.size()-1].pos - read.start + 1 << ":" << read.methyl[read.methyl.size()-1].type;
    }
    cout << "\t" << "MismatchData" << endl;
}


void simulator::simulate(){
	vector<MethylRead> data;
    readData();
    //cout << "DNA Lenght" << "\t" << dnaLength<< endl;
    buildMethylHap(dnaLength, pos, HapNum, methylHapVec, flag);
    //cout << "hap built" << endl;
    //cout << "max_i" << coverage*dnaLength/readLength << endl;

    for(int i=0; i < coverage*dnaLength/readLength; i++){
        int hap=0;
        int pos=0;
        MethylRead read;
        //cout << "dnaLength" << dnaLength << endl;

        selectHP(readLength, dnaLength, freq, hap, pos);
        //cout << "HP Seleced " << i << endl;
        //cout << "hap " << hap << endl;
        //cout << "pos " << pos << endl;
        //cout << "read Length" <<  readLength << endl;
        buildRead(hap, pos, readLength, error, methylHapVec, read);
		data.push_back(read);
		
        //cout << "read Wrote" << endl;
    }
	sort(data.begin(), data.end());
    //cout << "Read built" << endl;
    //cout << data[30].methyl.size() << endl;
	for(unsigned int i=0; i <data.size(); i++){
        //cout << i;
		writeMethylRead(data[i], i);
	}
}


