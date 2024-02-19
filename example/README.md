# MERGE Tutorial: Example for the Dataset [YAP1 HUMAN Fields2012-singles](https://doi.org/10.1073/pnas.1209751109)

To build a model of the fitness landscape of a protein and explore it in silico, the following files are required:
- protein sequence in fasta format (see [yap1.fasta](https://github.com/amillig/MERGE/blob/main/example/yap1.fasta))
- variant-fitness-pairs in csv format (see [yap1.csv](https://github.com/amillig/MERGE/blob/main/example/yap1.csv))

## 1. Generate A Multiple Sequence Alignment (MSA) Using Jackhmmer

To generate a multiple sequence alignment, the target sequence must be provided in fasta format and the inclusion threshold (--incT) has to be set. It is recommended to start with $\texttt{incT} = 0.5 \cdot L$, with $L$ being the sequence length. If the MSA contains too few sequences, increasing incT is recommended.
```bash
jackhmmer [-options] <seqfile> <seqdb>

Options:
- --incT <x> : consider sequences >= this score threshold as significant
- --cpu <n>  : number of parallel CPU workers to use for multithreads
- --noali    : don't output alignments, so output is smaller
- -A <f>     : save multiple alignment of hits to file <f>
```

For the dataset YAP1 ($L=34$) we set $\texttt{incT} = 17$ to generate the alignment:
```bash
jackhmmer --incT 17 --cpu 4 --noali -A yap1.sto yap1.fasta uniref100.fasta
```

## 2. Post-process the MSA
   
In a next step, the MSA is being post-processed by
- excluding all positions, where the wild type sequence has a gap,
- excluding all positions that contain more than 30 % gaps,
- excluding all sequences that contain more than 50 % gaps.

The script sto2a2m.py can be found [here](https://github.com/amillig/MERGE/tree/main/scripts/sto2a2m.py).
```bash
python sto2a2m.py -sto yap1.sto
```

## 3. Infer parameters for the Potts model using PLMC

Once the a2m file is generated, we proceed with infering the parameters of the statitstical model.
```bash
Usage: plmc [options] alignmentfile

Options:
- -o <f>  : save estimated parameters to file <f> (binary)
- -n <n>  : maximum number of threads to use in OpenMP
- -le <x> : set L2 lambda for couplings (e_ij)
- -lh <x> : set L2 lambda for fields (h_i)
- -m <n>  : maximum number of iterations
- -g      : model sequence likelihoods only by coding, non-gapped portions
- -f <s>  : select only uppercase, non-gapped sites from a focus sequence
```

For the YAP1 dataset we choose
- $48$ threads to use in OpenMP
- $\texttt{le} = 0.2 * (N_\text{sites} - 1) = 5.8$
- $\texttt{lh} = 0.01$
- $3500$ maximum iterations
- yap1 as focus sequence
```bash
plmc -o yap1.params -n 48 -le 5.8 -lh 0.01 -m 3500 -g -f yap1 yap1.a2m
```

## 4. Construct and explore the model of the fitness landscape using MERGE

Run the script [yap1_example.py](https://github.com/amillig/MERGE/blob/main/example/example_yap1.py).
```bash
python yap1_example.py
```

### Explanation of 'yap1_example.py'
First, import the required package
```python
import merge
import numpy as np
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
from scipy.stats import spearmanr
```

Next, user input is required:
- startingPosition : startingPosition != 1 if the sequence in the fasta file does not start at residue 1
- paramsFile       : binary file outputted by PLMC
- csvFile          : ensure that separator is ';' and not ','. The name of the column including the variants should be 'mutant'. Multiple substituted variants should be separated by ','. Variants should be given in the form F3A (PHE at position 3 is subsituted by ALA).
- fitnessColumn    : key of the column including the fitness values
```python
startingPosition = 170
paramsFile = 'yap1.params'
csvFile = 'yap1.csv'
fitnessColumn = 'linear' 
```

In a next step, sequence representations of the wild type [yap1_wt_encoded.npy](https://github.com/amillig/MERGE/blob/main/example/yap1_wt_encoded.npy) and all variants given in [yap1.csv](https://github.com/amillig/MERGE/blob/main/example/yap1.csv) are generated :
```python
encodeCls = merge.Encode(startingPosition, paramsFile)
np.save('yap1_wt_encoded.npy', encodeCls._encode_wt())
data = merge.get_data(fitnessColumn, csvFile, encodeCls)
merge.generate_dataframe(data, 'yap1_encoded.csv')
```

Finally, a model of the fitness landscape is generated using 80 % of the data and evaluated in terms of the coefficient of determination $R^2$ and Spearman's $\rho$.
```python
df = pd.read_csv('yap1_encoded.csv', sep=';')
xWt = np.load('yap1_wt_encoded.npy')
x = df.iloc[:,2:].to_numpy()
deltaE = merge.X_to_deltaE(x, xWt)
y = df.iloc[:,1].to_numpy()   

xTrain, xTest, deltaETrain, deltaETest, yTrain, yTest = train_test_split(x, deltaE, y, train_size=0.8, random_state=42) 
model = merge.CombinedPredictor()
model.train(xTrain, deltaETrain, yTrain)
 
print('R2 score (test set): %.2f '%(r2_score(yTest, model.predict(xTest, deltaETest))))
print('Spearman rho (test set): %.2f '%(spearmanr(yTest, model.predict(xTest, deltaETest))[0]))
```

Optionally, the fitness landscape can be explored to generate improved variants applying Metropolis Hastings algorithm. For this purpose, additional user input is required:
- yWt : Fitness value of the wild type
- maxSubstitutions : Maximum number of subsitutions to introduce during exploration
- factor : Only classify variants as improved if fitness >= factor * yWt
- nWalkers : Number of walkers to generate
- nCores : Number of CPUs to use for parallel execution
```python
yWt = 1.0
maxSubstitutions = 3
factor = 1.1
nWalkers = 96
nCores = 8

explore = merge.Explore(encodeCls, model, yWt, maxSubstitutions=maxSubstitutions, factor=factor)
results = explore.scrape_landscape(nWalkers, nCores=nCores)
df = pd.DataFrame(results, columns=['Variant', 'Predicted Fitness'])
df.to_csv('improved_variants.csv', index=False, sep=';')
```
