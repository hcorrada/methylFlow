# methylFlow


Cell-specific methylation pattern reconstruction. Currently uses an LP
formulation and solver. The lemon graph library
http://lemon.cs.elte.hu/trac/lemon
and the glpk LP solver
http://www.gnu.org/software/glpk/
are included in this repository. The `ezOptionParser.hpp` is also
included
http://ezoptionparser.sourceforge.net/

The algorithm is described and tested in this publication: [http://bioinformatics.oxfordjournals.org/content/32/11/1618.abstract)(http://bioinformatics.oxfordjournals.org/content/32/11/1618.abstract)

## Installation


This project uses `cmake` for building and requires at least
version 2.6. It also uses `c++11` so use a compiler that supports
this (e.g., g++ >= 4.7 or clang >= 3.4)

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

To compile with DEBUG flags use

```shell
...
$ mkdir build_devel && cd build_devel
$ cmake -DCMAKE_BUILD_TYPE=Debug ..
$ make
...
```
## Usage

<pre>
MethylFlow: methylation pattern reconstruction

USAGE: methylFlow -sam -i reads.sam -o mfoutput [OPTIONS]

OPTIONS:

-chr, -Chr ARG                    chr name for tsv files, not required for sam
                                  input file

-cpgloss, -p, -P, --cpgloss       Use cpg-loss instead of region-loss.

-e, -eps, -E, --eps ARG           Regularization parameter search threshold.

-end, -End, --end ARG             Only process reads aligning before given
                                  location.

-h, -help, --help, --usage        Display usage instructions.

-i, -in, --in, --input ARG        Read input file. Default:Tab-separated format:
                                  start length strand
                                  methyl<string>(offset<int>[M|U]
                                  substitutions<string>(ignored))

-l, -lam, -lambda, --lambda ARG   Regularization parameter value.

-o, -out, --out, --output ARG     Output directory. Directory must exist before
                                  running. Files written: cpgs.tsv, components.tsv,
                                  patterns.tsv, regions.tsv

-s, -scale, -S, --scale ARG       Scale parameter value.

-sam, -SAM, --sam                 Input file is in SAM format instead of default
                                  tab-separated format.

-start, -Start, --start ARG       Only process reads aligning after given
                                  location.

-v, -verbose, -V, --verbose       Verbose option.

EXAMPLES:

methylFlow -sam -i reads.sam -o mfoutput -l 10.0 -s 30.0 -e 0.1
</pre>

### Output

Upon running, the output directory (`mfoutput` in the example above) will contain three files with the following format:

#### `cpgs.txt`

Tab-separated file of coverage and methylation calls per cpg. Columns
 - `chr`: chromosome name
 - `pos`: cpg position
 - `Cov`: number of reads overlapping CpG
 - `Meth`: number of reads indicating CpG is methylated


#### `components.tsv`

Tab-separated file of components found by algorithm. A component is a connected region graph based on overlapping reads. Genomic regions are covered by a single component, thus, cell-specific patterns estimated in a given genomic region are obtained from (one or morei non-overlapping) components that overlap that region. 

Columns:

- `chr`: chromosome name of genomic region covered by connected component
- `start`: starting position of genomic region covered by connected component
- `end`: ending position of genomic region covered by connected component
- `cid`: component id, identifier given to component, used to connect to regions and patterns in other output files
- `npatterns`: number of cell-specific methylation patterns estimated from this connected component.  
- `total_coverage`: total number of reads overlapping this component's genomic region
- `total_flow`: the sum of all estimated abundances (flows) for patterns in this region

#### `patterns.tsv`

Tab-separated file of cell-specific methylation patterns estimated by `methylFlow`. 

Columns:

- `chr`: chromosome name 
- `start`: start position of pattern
- `end`: end position of pattern
- `cid`: component id, corresponds to id of a component in file `components.tsv`
- `pid`: pattern id, identifier given to pattern (unique across patterns within the same component)
- `abundance`: abundance estimated for this pattern
- `methylpat`: comma-separated list of methylation status entries  of cpgs within pattern. Entries are `pos:[M|U]` where position is the location of the CpG from the start of the pattern and `M|U` indicates if the CpG is methylated or unmethylated respectively
- `regions`: comma-separated list of regions included in pattern (see file `regions.tsv`)

#### `regions.tsv`

Tab-separated file of regions that make up the region graph used in the estimation algorithm. Reads are assigned to a region
if they have no disagreement on their methylation pattern. That is, regions contain the longest stretches of overlapping reads with unambiguous methylation patterns.



Columns:

- `chr`: chromosome name
- `start`: start position of region
- `end`: end position of region
- `cid`: component id, corresponds to identified of component in file `components.tsv`
- `rid`: region id, identifier given to region (unique across regions within the same component)
- `raw_coverage`: number of reads assigned to the region
- `norm_coverage`: normalized region coverage
- `exp_coverage`: the sum of abundances of all patterns that include this region
- `methylpat`: methylation pattern of region, given in same format as `patterns.tsv`

##Authors

Hector Corrada Bravo <hcorrada@gmail.com>  
Faezeh Dorri  

Center for Bioinformatics and Computational Biology  
University of Maryland  
http://www.cbcb.umd.edu/~hcorrada
