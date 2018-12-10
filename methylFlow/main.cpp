/* -*- mode: C++; indent-tabs-mode: nil; -*-
 *
 * This file is a part of LEMON, a generic C++ optimization library.
 *
 * Copyright (C) 2003-2009
 * Egervary Jeno Kombinatorikus Optimalizalasi Kutatocsoport
 * (Egervary Research Group on Combinatorial Optimization, EGRES).
 *
 * Permission to use, modify and distribute this software is granted
 * provided that this copyright notice appears in all copies. For
 * precise terms see the accompanying LICENSE file.
 *
 * This software is provided "AS IS" with no warranty of any kind,
 * express or implied, and with no claim as to its suitability for any
 * purpose.
 *
 */

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cassert>
#include <cstdlib>
#include <climits>
#include <stdio.h>
#include <sys/time.h>


#include <vector>
#include <lemon/lp.h>
#include <lemon/list_graph.h>
#include <lemon/maps.h>

#include "ezOptionParser.hpp"
#include "mflib/MFGraph.hpp"

using namespace methylFlow;
using namespace ez;

void Usage(ezOptionParser & opt) {
    std::string usage;
    opt.getUsage(usage);
    std::cout << usage << std::endl;
};

//The main entry point
// assume arguments:
//   input file
//   component output file
//   pattern output file
//   region output file
//   lambda
//   scale multiplier
//
// writes in tsv files
//   components: chr start end componentID numPatterns totalCoverage totalFlow<float> (scaled coverage)
//   patterns: chr start end componentID patternID totalFlow<float> methylationPattern (offset<int>:[M|U])
//   regions: chr start end componentID regionID coverage expectedCoverage<float> methylationPattern (offset<int>:[M|U])
int main(int argc, const char **argv)
{
    std::stringstream buffer;

    ezOptionParser opt;
    opt.overview = "MethylFlow: methylation pattern reconstruction";
    opt.syntax = "methylFlow -sam -i reads.sam -o mfoutput [OPTIONS]";
    opt.example = "methylFlow -sam -i reads.sam -o mfoutput -l 10.0 -s 30.0 -e 0.1";

    // help and usage
    opt.add(
            "", //Default
            0, // not required
            0, // no args expected
            0, // no delimiter
            "Display usage instructions.", // help description
            "-h", // flag tokens
            "-help",
            "--help",
            "--usage"
            );

    // input file
    opt.add(
            "", // Default.
            0, // Required (for now, will switch to stdin if missing in future)
            1, // number of args expected
            0, // delimiter, not needed
            "Read input file. Default:Tab-separated format:\nstart length strand methyl<string>(offset<int>[M|U] substitutions<string>(ignored))", // Help description
            "-i", //flag token
            "-in", // flag token
            "--in", // flag token
            "--input" //flag token
            );

    // chr flag
    opt.add(
            "", // Default.
            0, // Required (for now, will switch to stdin if missing in future)
            1, // number of args expected
            0, // delimiter, not needed
            "chr name for tsv files, not required for sam input file", // Help description
            "-chr", //flag token
            "-Chr" // flag token
            );

    //start position
    const float DEFAULT_START = 0;
    buffer.str("");
    buffer << DEFAULT_START;
    opt.add(
            buffer.str().c_str(), //Default
            0, // not required
            1, // number of args expected
            0, // no delimiter
            "Only process reads aligning after given location.", // help description
            "-start",
            "-Start",
            "--start"// flag tokens
            );

    //end position
    const float DEFAULT_END = LONG_MAX;
    buffer.str("");
    buffer << DEFAULT_END;
    opt.add(
            buffer.str().c_str(), //Default
            0, // not required
            1, // number of args expected
            0, // no delimiter
            "Only process reads aligning before given location.", // help description
            "-end",
            "-End",
            "--end"// flag tokens
            );


    // Sam file input
    const bool DEFAULT_FLAG_SAM = false;
    opt.add(
            "", // Default.
            0, // Required (for now, will switch to stdin if missing in future)
            0, // number of args expected
            0, // delimiter, not needed
            "Input file is in SAM format instead of default tab-separated format.", // Help description
            "-sam", // flag token
            "-SAM", // flag token
            "--sam" //flag token
            );


    const char * DEFAULT_OUTDIR = "mfoutput";
    // output directory
    opt.add(
            DEFAULT_OUTDIR, // Default
            0, // Not required, uses default
            1, // number of args
            0, // delimiter
            "Output directory. Directory must exist before running. Files written: components.tsv, patterns.tsv, regions.tsv", // help description
            "-o", // flag tokens
            "-out",
            "--out",
            "--output"
            );

    // lambda
    const float DEFAULT_LAMBDA = -1.0;
    buffer.str("");
    buffer << DEFAULT_LAMBDA;
    opt.add(
            buffer.str().c_str(), //default
            0, // not required uses default
            1, // number of args,
            0, //no delimiter
            "Regularization parameter value.", // help description
            "-l", // flag tokens
            "-lam",
            "-lambda",
            "--lambda"
            );

    // scale parameter
    const float DEFAULT_SCALE = 1.0;
    buffer.str("");
    buffer << DEFAULT_SCALE;
    opt.add(
            buffer.str().c_str(), //default
            0, // not required, uses default
            1, // num args
            0, // no delimiter
            "Scale parameter value.", // help description
            "-s", //flag tokens
            "-S",
            "-scale",
            "--scale"
            );

    // epsilon parameter
    const float DEFAULT_EPSILON = 0.1;
    buffer.str("");
    buffer << DEFAULT_EPSILON;
    opt.add(
            buffer.str().c_str(), // default
            0, // not required, uses default
            1, // num args
            0, // no delimiter
            "Regularization parameter search threshold.", // help description
            "-e", // flag tokens
            "-E",
            "-eps",
            "--eps"
            );

    // verbose option
    const bool DEFAULT_VERBOSE = false;
    buffer.str("");
    opt.add(
            "", // default
            0, // no required uses default
            0, // no args, it's a flag
            0, // no delimiter
            "Verbose option.", // help description
            "-v", // flag tokens
            "-V",
            "-verbose",
            "--verbose"
            );
    
    // verbose for elapsed time option
    const bool DEFAULT_VERBOSE_TIME = false;
    buffer.str("");
    opt.add(
            "", // default
            0, // no required uses default
            0, // no args, it's a flag
            0, // no delimiter
            "VerboseTime option.", // help description
            "-vt", // flag tokens
            "-VT",
            "-verboseTime",
            "--verboseTime"
            );
    
    // use cpg-loss algorithm instead
    // of region-loss
    const bool DEFAULT_PCTSELECT = false;
    opt.add(
            "", // default (don't do it)
            0, // not required uses default
            0, // no args, it's a flag
            0, // no delimiter
            "Use cpg-loss instead of region-loss.", // help description
            "-p", // flag tokens
            "-P",
            "-cpgloss",
            "--cpgloss"
            );

    const bool DEFAULT_GRAPH_ONLY = false;
    opt.add(
            "", // default (don't do it)
            0, // not required uses default
            0, // no args, it's a flag
            0, // no delimiter
            "Print graph info out only, do not solve optimization problem.", // help description
            "-gr", // flag tokens
            "-GR",
            "-graph-only",
            "--graph-only"
            );

    opt.parse(argc, argv);

    if (opt.isSet("-h")) {
        Usage(opt);
        return 0;
    }


    std::vector<std::string> badOptions;
    int i;
    if (!opt.gotRequired(badOptions)) {
        for (i=0; i < badOptions.size(); ++i)
            std::cerr << "ERROR: Missing required for option " << badOptions[i] << ".\n\n";
        Usage(opt);
        return 1;
    }

    if (!opt.gotExpected(badOptions)) {
        for (i=0; i < badOptions.size(); ++i)
            std::cerr << "Error: Got unexpected number of arguments for option " << badOptions[i] << ".\n\n";
        Usage(opt);
        return 1;
    }

    int status = 0;
    std::istream* instream;
    std::string input_filename;
    std::ifstream input;
    if (opt.isSet("-i")) {
        opt.get("-i")->getString(input_filename);
        input.open( input_filename.c_str() );
        instream = &input;
        if (!instream) status = -1;
    } else {
        instream = &std::cin;
    }

    std::string chr;
    if (opt.isSet("-chr")) {
        opt.get("-chr")->getString(chr);
    }


    std::string outdirname = DEFAULT_OUTDIR;
    if (opt.isSet("-o")) {
        opt.get("-o")->getString(outdirname);
    }

    bool flag_SAM = DEFAULT_FLAG_SAM;
    if (opt.isSet("-sam")) {
        flag_SAM = true;
    }

    float lambda = DEFAULT_LAMBDA;
    if (opt.isSet("-l")) {
        opt.get("-l")->getFloat(lambda);
    }

    long start = DEFAULT_START;
    if (opt.isSet("-start")) {
        opt.get("-start")->getLong(start);
    } 

    long end = DEFAULT_END;
    if (opt.isSet("-end")) {
        opt.get("-end")->getLong(end);
    } 


    float scale_mult = DEFAULT_SCALE;
    if (opt.isSet("-s")) {
        opt.get("-s")->getFloat(scale_mult);
    }

    float epsilon = DEFAULT_EPSILON;
    if (opt.isSet("-e")) {
        opt.get("-e")->getFloat(epsilon);
    }

    bool verbose = DEFAULT_VERBOSE;
    if (opt.isSet("-v")) {
        verbose = true;
    }
    
    bool verboseTime = DEFAULT_VERBOSE_TIME;
    if (opt.isSet("-vt")) {
        verboseTime = true;
    }
    
    bool pctselect = DEFAULT_PCTSELECT;
    if (opt.isSet("-p")) {
      pctselect = true;
    }

    bool graph_only = DEFAULT_GRAPH_ONLY;
    if (opt.isSet("-gr")) {
      graph_only = true;
    }

    buffer.str("");
    std::ofstream comp_stream;
    buffer << outdirname << "/components.tsv";
    comp_stream.open( buffer.str().c_str(), std::ofstream::out | std::ofstream::trunc );
    if (!comp_stream) status = -1;

    buffer.str("");
    std::ofstream pattern_stream;
    buffer << outdirname << "/patterns.tsv";
    if (!graph_only) pattern_stream.open( buffer.str().c_str(), std::ofstream::out | std::ofstream::trunc );
    if (!pattern_stream) status = -1;

    buffer.str("");
    std::ofstream region_stream;
    buffer << outdirname << "/regions.tsv";
    region_stream.open( buffer.str().c_str(), std::ofstream::out | std::ofstream::trunc );
    if (!region_stream) status = -1;

    buffer.str("");
    std::ofstream cpg_stream;
    buffer << outdirname << "/cpgs.tsv";
    cpg_stream.open( buffer.str().c_str(), std::ofstream::out | std::ofstream::trunc );
    if (!cpg_stream) status = -1;

    buffer.str("");
    std::ofstream edge_stream;
    buffer << outdirname << "/edges.tsv";
    if (graph_only) edge_stream.open( buffer.str().c_str(), std::ofstream::out | std::ofstream:: trunc );
    if (!edge_stream) status = -1;

    buffer.str("");
    std::ofstream path_stream;
    buffer << outdirname << "/paths.tsv";
    if (graph_only) path_stream.open( buffer.str().c_str(), std::ofstream::out | std::ofstream:: trunc );
    if (!path_stream) status = -1;
    if (status == -1) {
        std::cerr << "[methylFlow] Error opening file." << std::endl;
        return 1;
    }

    MFGraph g;
    struct timeval tvalBefore, tvalAfter;

    gettimeofday (&tvalBefore, NULL);
    status = g.run( *instream,
                   comp_stream,
                   pattern_stream,
                   region_stream,
                   cpg_stream,
                   edge_stream,
                    path_stream,
                   chr,
                   start,
                   end,
                   flag_SAM,
                   lambda,
                   scale_mult,
                   epsilon,
                   verbose,
                   verboseTime,
                   pctselect,
                   graph_only);
    gettimeofday (&tvalAfter, NULL);
    printf("Time in miliseconds total: %ld miliseconds\n", ((tvalAfter.tv_sec - tvalBefore.tv_sec)*1000  + tvalAfter.tv_usec/1000) - tvalBefore.tv_usec/1000);


    // streams are closed when object
    // is destroyed
    return status;
}
