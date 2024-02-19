# MERGE
MERGE represents a method that combines direct coupling analysis and machine learning techniques to predict a protein's fitness from sequence. It requires a binary parameter file outputted by [plmc](https://github.com/debbiemarkslab/plmc/tree/master) and sequence-fitness pairs.

![MERGE_git](https://github.com/amillig/MERGE/assets/58852023/f3da6124-5bee-41a2-b4be-a9c8cd0c4947)

# Usage
In the following, the most important steps for model construction are briefly described. A step-by-step guide is given in [example](https://github.com/amillig/MERGE/tree/main/example). To build a model of the fitness landscape of a protein and explore it *in silico*, the following files are required:
- protein sequence in fasta format
- variant-fitness-pairs in csv format

1. Generate a Multiple Sequence Alignment (MSA) Using Jackhmmer

To generate a multiple sequence alignment, the target sequence must be provided in fasta format and the inclusion threshold (--incT) must be set.

```bash
jackhmmer [-options] <seqfile> <seqdb>
```

2. Post-process the MSA
   
In a next step, the MSA is being post-processed by
- excluding all positions, where the wild type sequence has a gap,
- excluding all positions that contain more than 30 % gaps,
- excluding all sequences that contain more than 50 % gaps.

The script sto2a2m.py can be found [here](https://github.com/amillig/MERGE/tree/main/scripts/sto2a2m.py).
```bash
python sto2a2m.py -sto <stoFile>
```

3. Infer Parameters for Potts-Model Using PLMC
Once the a2m file is generated, the parameters of the statitstical model are inferred.
```bash
plmc [options] alignmentfile
```

4. Construct and Explore the Model of the Fitness Landscape Using MERGE

Finally, a model of the fitness landscape is generated. See the example for details on how to use merge.

# Prerequisites
  ### 1. Download the latest version of UniRef100
  1. Download the latest version of UniRef100
```bash
wget https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz
```

  2. Unzip the file to get a fasta file
```bash
gzip -d uniref100.fasta.gz
```

For further information see [uniprot help](https://www.uniprot.org/help/downloads)

### 2. Installing HMMER
1. Download the tarball
```bash
wget http://eddylab.org/software/hmmer/hmmer.tar.gz
```

2. Unpack the tarball
```bash
tar zxf hmmer.tar.gz
```

3. Enter the directory 'hmmer-3.4'
```bash
cd hmmer-3.4
```

4. Set the installaion path (adjust "/your/install/path" accordingly!)
```bash
./configure --prefix /your/install/path
```

5. Build HMMER
```bash
make
```

6. Run self tests (optional)
```bash
make check
```

7. Install programs and man pages
```bash
make install
```

8. Add executable to PATH for session (adjust "/your/install/path" accordingly!)
```bash
export PATH="/your/install/path/bin:$PATH"
```
or permanently (adjust "/your/install/path" accordingly!)
```bash
echo 'export PATH="/your/install/path/bin:$PATH"' >> ~/.bashrc
```

For further information see [hmmer documentation](http://hmmer.org/documentation.html)

### 3. Installing PLMC
1. Clone the plmc repository 
```bash
git clone https://github.com/debbiemarkslab/plmc.git
```

2. Enter thr plmc directory
```bash
cd plmc
```

3. Build with GCC and OpenMP to enable multicore parallelism
```bash
make all-openmp
```

4. Add executable to PATH for session (adjust "/your/install/path" accordingly!)
```bash
export PATH="/your/install/path/plmc/bin:$PATH"
``` 
or permanently (adjust "/your/install/path" accordingly!)
```bash
echo 'export PATH="/your/install/path/plmc/bin:$PATH"' >> ~/.bashrc
```

For further information see [plmc repository](https://github.com/debbiemarkslab/plmc)


### 4. MERGE
1. Clone the MERGE repository
```bash
git clone https://github.com/amillig/MERGE.git
```

2. Enter the MERGE directory
```bash
cd MERGE
```

3. Install the dependencies
```bash
pip install -r requirements.txt
```

4. Import MERGE as module in Python
```python
import merge
```

# References
“Combining evolutionary probability and machine learning enables data-driven protein engineering with minimized experimental effort” by Alexander-Maurice Illig, Niklas E. Siedhoff, Mehdi D. Davari*, and Ulrich Schwaneberg*

# Author
MERGE was developed and written by [Alexander-Maurice Illig](https://www.researchgate.net/profile/Alexander-M-Illig) at RWTH Aachen University.

# Credits
MERGE uses binary parameter files that are generated with [plmc](https://github.com/debbiemarkslab/plmc) written by John Ingraham.

# License
[License](https://github.com/amillig/MERGE/blob/main/LICENSE)
