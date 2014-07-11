// to be run on cbcb server : qsub run.sh -t 3-5 -q xlarge -l mem=24G,walltime=24:00:00 -N Hap



#include <iostream>
#include <algorithm>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include "simulator.h"
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
    unsigned long seed = mix(clock(), time(NULL), getpid());
    srand(seed);
    cerr << "rand check " << rand()%1000 << endl;
    simulator * sim= new simulator();
    sim->simulate();
    delete sim;
	return 1;
}


