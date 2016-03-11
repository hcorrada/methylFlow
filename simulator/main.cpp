// to be run on cbcb server : qsub run.sh -t 3-5 -q xlarge -l mem=24G,walltime=24:00:00 -N Hap



#include <iostream>
#include <algorithm>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>

#include "simulator.h"

std::ifstream inputFile;

unsigned long mix(unsigned long a, unsigned long b, unsigned long c)
{
    a=a-b;  a=a-c;  a=a^(c >> 13);
    b=b-c;  b=b-a;  b=b^(a << 8);
    c=c-a;  c=c-b;  c=c^(b >> 13);
    a=a-b;  a=a-c;  a=a^(c >> 12);
    b=b-c;  b=b-a;  b=b^(a << 16);
    c=c-a;  c=c-b;  c=c^(b >> 5);
    a=a-b;  a=a-c;  a=a^(c >> 3);
    b=b-c;  b=b-a;  b=b^(a << 10);
    c=c-a;  c=c-b;  c=c^(b >> 15);
    return c;
}

int main (int argc, char* argv[]) {
    cerr << "start simulation" << endl;
    unsigned long seed = mix(clock(), time(NULL), getpid());
    srand(seed);
    cerr << "rand check " << rand()%1000 << endl;


    if(argc < 3){
		cout << "Please enter your input" << endl;
		return -1;
	}

    inputFile.open(argv[1]);
    std::string outdirname = argv[2];
    std::stringstream buffer;

    buffer.str("");
    std::ofstream patternFile;
    buffer << outdirname << "/simPattern.txt";
    patternFile.open( buffer.str().c_str() );
    patternFile << "chr" << "\t" << "startDNA" << "\t" << "dnaLength" << "\t" << "ComponentID"  << "\t" << "PatternID"<< "\t" << "PatternFreq" << "\t" << "methylInfo" << endl;


    buffer.str("");
    std::ofstream shortReadFile;
    buffer << outdirname << "/shortRead.txt";
    shortReadFile.open( buffer.str().c_str() );
    //shortReadFile << "readID" << "\t" << "start"  << "\t" << "length" << "\t" << "W" << "\t" << "methylInfo" << endl;



    //patternFile.open(argv[2]);
   // shortReadFile.open(argv[3],std::ios_base::app);

    simulator * sim= new simulator();
    sim->simulate(inputFile, patternFile, shortReadFile);
    inputFile.close();
    patternFile.close();
    shortReadFile.close();
    delete sim;
	return 1;
}
