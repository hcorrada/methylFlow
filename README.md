# methylFlow


Cell-specific methylation pattern reconstruction. Currently uses an LP
formulation and solver. The lemon graph library
http://lemon.cs.elte.hu/trac/lemon
and the glpk LP solver
http://www.gnu.org/software/glpk/
are included in this repository. The `ezOptionParser.hpp` is also
included
http://ezoptionparser.sourceforge.net/

## Installation


This project uses `cmake` for building and required at least
version 2.6.

```shell
$ git clone https://github.com/hcorrada/methylFlow.git
$ cd methylFlow
$ git submodule init
$ git submodule update
$ mkdir build && cd build
$ cmake ..
$ make
$ make install
```

## Usage

<pre>
USAGE: methylFlow -i reads.tsv -o mfoutput [OPTIONS]

OPTIONS:

-chr, -Chr ARG                    chr name for tsv files, not required for sam
                                  input file

-cpgloss, -p, -P, --cpgloss       Use cpg-loss instead of region-loss.

-e, -eps, -E, --eps ARG           Regularization parameter search threshold.

-end, -End, --end ARG             Display end instructions.

-h, -help, --help, --usage        Display usage instructions.

-i, -in, --in, --input ARG        Read input file. Default:Tab-separated format:
                                  start length strand
                                  methyl<string>(offset<int>[M|U]
                                  substitutions<string>(ignored))

-l, -lam, -lambda, --lambda ARG   Regularization parameter value.

-o, -out, --out, --output ARG     Output directory. Files written:
                                  components.tsv, patterns.tsv, regions.tsv

-s, -scale, -S, --scale ARG       Scale parameter value.

-sam, -SAM, --sam                 Input file is in SAM format instead of default
                                  tab-separated format.

-start, -Start, --start ARG       Display start instructions.

-v, -verbose, -V, --verbose       Verbose option.

EXAMPLES:

methylFlow -i reads.tsv -o mfoutput -l 10.0 -s 30.0 -e 0.1
</pre>



##Authors

Hector Corrada Bravo <hcorrada@gmail.com>  
Faezeh Dorri  

Center for Bioinformatics and Computational Biology  
University of Maryland  
http://www.cbcb.umd.edu/~hcorrada
