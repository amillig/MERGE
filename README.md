# MERGE

MERGE represents a method that combines direct coupling analysis and machine learning techniques to predict a protein's fitness from sequence. It requires a binary parameter file outputted by [plmc](https://github.com/debbiemarkslab/plmc/tree/master) and sequence-fitness pairs.

![MERGE_git](https://github.com/amillig/MERGE/assets/58852023/f3da6124-5bee-41a2-b4be-a9c8cd0c4947)

# Prerequisites

#### Download latest version of UniRef100 
1. Download the latest version of UniRef100
`$ wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz`
2. Unzip the file to get a fasta file `$ gzip -d uniref100.fasta.gz`

For further information see https://www.uniprot.org/help/downloads

#### Installing HMMER
1. Download tarball `$ wget http://eddylab.org/software/hmmer/hmmer.tar.gz`
2. Unpack the tarball `$ tar zxf hmmer.tar.gz`
3. Enter directory `$ cd hmmer-3.4`
4. Set installaion path `$ ./configure --prefix /your/install/path`
5. Build HMMER `$ make`
6. Run self tests `$ make check`
7. Install programs and man pages `$ make install`
8.  Add executable to PATH `$ export PATH="/your/install/path/bin:$PATH"`

For further information see http://hmmer.org/documentation.html

#### Installing PLMC
1. Clone plmc repository `$ git clone https://github.com/debbiemarkslab/plmc.git`
2. `$ cd plmc`
3. Build with gcc and OpenMP to enable multicore parallelism `$ make all-openmp`
4. Add executable to PATH `$ export PATH="/your/install/path/plmc/bin:$PATH"`

For further information see https://github.com/debbiemarkslab/plmc

#### MERGE
1. `$ git clone https://github.com/amillig/MERGE.git`
2. `$ cd MERGE`
3. `$ pip install -r requirements.txt` ***AI Requirements file erstellen***



# Workflow

***AI Workflow figure***


Sequence -> jackhmmer -> sto2a2m.py -> a2m -> post_process.py -> PLMC -> MERGE



# Example for the Dataset YAP1_HUMAN_Fields2012-singles

To build a model of the fitness landscape of a protein and explore it *in silico*, the following files are required:

- protein sequence in fasta format (see yap1.fasta)
- variant-fitness-pairs in csv format (see yap1.csv)

### 1. Generate Multiple Sequence Alignment (MSA) Using Jackhmmer

To generate a multiple sequence alignment, the target sequence must be provided in fasta format and the inclusion threshold (--incT) must be set. Inclusion threshold: Start with $0.5 \cdot L$, with $L$ being the sequence length. If the MSA contains too few sequences, increasing incT is recommended. 

Usage: `jackhmmer [-options] <seqfile> <seqdb>`

Options:
    - --incT <x> : consider sequences >= this score threshold as significant
    - --cpu <n>  : number of parallel CPU workers to use for multithreads
    - --noali    : don't output alignments, so output is smaller
    - -A <f>     : save multiple alignment of hits to file <f>

For the YAP1_HUMAN_Fields2012-singles dataset ($L=34$) we choose use an inclusion threshold of 17 (--incT) to generate the alignment:

`$ jackhmmer --incT 17 --cpu 4 --noali -A yap1.sto yap1.fasta uniref100.fasta`


### 2. Post-process the MSA

In a next step, we curate the MSA by 

1. excluding all positions, where the wild type sequence has a gap,
2. excluding all positions that contain more than 30 % gaps,
3. excluding all sequences that contain more than 50 % gaps.

The script `sto2a2m.py` can be found in the git repository MERGE.

`$ python sto2a2m.py -sto yap1.sto`

*yap1.a2m* is generated as an output file.


### 3. Infer Parameters for Potts-Model Using PLMC

Once the a2m file is generated, we proceed with infering the parameters of the statitstical model. 

Usage: `plmc [options] alignmentfile`

Options:
    - -o <f>  : save estimated parameters to file <f> (binary)
    - -n <n>  : maximum number of threads to use in OpenMP
    - -le <x> : set L2 lambda for couplings (e_ij)
    - -lh <x> : set L2 lambda for fields (h_i)
    - -m <n>  : maximum number of iterations
    - -g      : model sequence likelihoods only by coding, non-gapped portions
    - -f <s>  : select only uppercase, non-gapped sites from a focus sequence

For the YAP1_HUMAN_Fields2012-singles dataset we choose
    - 48 threads to use in OpenMP
    - le = 0.2 * (N_sites - 1) = 5.8
    - lh = 0.01.
    - 3500 as maximum iterations
    - yap1 as focus sequence
   
`$ plmc -o yap1.params -n 48 -le 5.8 -lh 0.01 -m 3500 -g -f yap1 yap1.a2m`


### 4. Construct and Explore Model of Fitness Landscape Using MERGE

# References
“Combining evolutionary probability and machine learning enables data-driven protein engineering with minimized experimental effort” by Alexander-Maurice Illig, Niklas E. Siedhoff, Mehdi D. Davari*, and Ulrich Schwaneberg*

# Author
MERGE was developed and written by Alexander-Maurice Illig at RWTH Aachen University.

# Credits
MERGE uses binary parameter files that are generated with [plmc](https://github.com/debbiemarkslab/plmc/tree/master) written by [John Ingraham](https://github.com/jingraham).
