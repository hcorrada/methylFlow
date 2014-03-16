// to be run on cbcb server : qsub run.sh -t 3-5 -q xlarge -l mem=24G,walltime=24:00:00 -N Hap



#include <iostream>
#include <algorithm>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include "simulator.h"

int main (int argc, char* argv[]) {
    
    simulator * sim= new simulator();
    sim->simulate();
    delete sim;
	return 1;
}


