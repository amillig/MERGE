# version         v0.1.7
# date            30.01.2024
# author          Alexander-Maurice Illig
# affilation      Institute of Biotechnology, RWTH Aachen
# email           a.illig@biotec.rwth-aachen.de

import numpy as np
import pandas as pd
import multiprocessing

from ._utils import get_single_substitutions
from ._errors import ActiveSiteError

class PLMC:
    """
    Class for extracting data from 'params' file outputted by PLMC.
    Based on https://github.com/debbiemarkslab/plmc/blob/master/scripts/read_params.m
    
    Attributes
    ----------
    paramsFile: str
        Binary parameter file outputed by PLMC.
    """
    def __init__(self, paramsFile:str):
        self._read_paramsFile(paramsFile)
        self.alphabet2index = {aminoAcid:i for i,aminoAcid in enumerate(self.alphabet)}
        self.position2index = {position:i for i,position in enumerate(self.offsetMap)}
        
    def _read_paramsFile(self, paramsFile:str):
        """
        Description
        -----------
        Initializes the PLMC class by extracting information from the binaray 'params' file.
        
        Parameters
        ----------
        paramsFile : str
            Name of the binary 'params' file.
        """
        with open(paramsFile, 'rb') as f:
            self.L,self.q,self.numSeqs,self.numInvalidSeqs,self.numIter = np.fromfile(f, dtype='int32', count=5)
            self.theta,self.lambdaH,self.lambdaJ,self.lambdaGroup,self.nEff = np.fromfile(f, dtype='float32', count=5)
            self.alphabet = np.fromfile(f, dtype='S1', count=self.q).astype('U1')
            self.weights = np.fromfile(f, dtype='float32', count=self.numSeqs+self.numInvalidSeqs)
            self.targetSeq = np.fromfile(f, dtype='S1', count=self.L).astype('U1')
            self.offsetMap = np.fromfile(f, dtype='int32', count=self.L)
            self.fi, = np.fromfile(f, dtype=('float32', (self.L, self.q)), count=1)
            self.hi, = np.fromfile(f, dtype=('float32', (self.L, self.q)), count=1)

            self.fij = np.zeros((self.L, self.L, self.q, self.q), dtype='float32')
            for i in range(self.L-1):
                for j in range(i+1, self.L):
                    self.fij[i,j], = np.fromfile(f, dtype=('float32', (self.q, self.q)), count=1)
                    self.fij[j,i] = self.fij[i,j].T

            self.Jij = np.zeros((self.L, self.L, self.q, self.q), dtype=('float32'))
            for i in range(self.L-1):
                for j in range(i+1, self.L):
                    self.Jij[i,j], = np.fromfile(f, dtype=('float32', (self.q, self.q)), count=1)     
                    self.Jij[j,i] = self.Jij[i,j].T
                    
class Encode(PLMC):
    """
    Class for performing the 'DCA-based encoding'.

    Attributes
    ----------
    startingPosition: int
        Number of leading residue of the fasta sequence used for model construction.
    paramsFile: str
        Binary parameter file outputed by PLMC.
    """

    def __init__(self,startingPosition:int, paramsFile:str):
        self.startingPosition=startingPosition
        super().__init__(paramsFile) # inherit functions and variables from class 'PLMC'
        self.xWt = self._encode_wt()

    def _get_position_internal(self, position:int):
        """
        Description
        -----------
        Returns the "internal position" of an amino acid, e.g. D19V is the desired substitution,
        but the fasta sequence starts from residue 3, i.e. the first two residues are "missing".
        The DCA model will then recognize D19 as D17. In order to avoid wrong assignments,
        it is inevitable to calculate the "internal position" 'i'.

        Parameters
        ----------
        position : int
            Position of interest.
        
        Returns
        -------
        i : int
            "Internal position" that may differ due to different starting residue.
        None
            If the requested position is not an active site.
        """
        offset=self.startingPosition-1
        i=position-offset
        if i in self.offsetMap:
            return i
        else:
            return None

    def Ji(self, i:int, Ai_index:str, sequence:np.ndarray) -> float:
        """
        Description
        -----------
        Caluclates the sum of all site-site interaction terms when site 'i' is occupied with amino acid 'Ai'.

        Parameters
        ----------
        i : int
            Index of the position.
        Ai_index : int
            Index of the introduced amino acid at 'i'.
        sequence: np.ndarray
            Sequence of the variant as numpy array.

        Returns
        -------
        Ji : float
            Sum of all site-site interaction terms acting on position 'i' when occupied with 'Ai'.
        """
        Ji=0.0
        for j,Aj in enumerate(sequence):
            Aj_index = self.alphabet2index[Aj]        
            Ji+=self.Jij[i,j,Ai_index,Aj_index]
        return Ji

    @staticmethod
    def _unpack_substitution(substitution:str) -> tuple:
        """
        Description
        -----------
        Converts string representation of variant into tuple.

        Parameters
        ----------
        substitution : str
            Substitution as string: Integer enclosed by two letters representing
            the wild-type (first) and variant amino acid (last) in one letter code.

        Returns
        -------
        substitution : tuple
            (wild-type amino acid, position, variant amino acid)
        """
        return (substitution[0],int(substitution[1:-1]),substitution[-1])

    def _encode_variant(self, variant:str, separator=',') -> np.ndarray:
        """
        Description
        -----------
        Encodes the variant using its "DCA representation".

        Parameters
        ----------
        variant : str
            Joined string of integers enclosed by two letters representing the wild-type
            and variant amino acid in the single letter code. -> Check separator
        separator : str
            Character to split the variant to obtain the single substitutions (default=',').
        
        Returns
        -------
        X_var : np.ndarray
            Encoded sequence of the variant.
        """
        sequence=self.targetSeq.copy()
        for substitution in get_single_substitutions(variant,separator):
            wild_type_aa,position,Ai=self._unpack_substitution(substitution)
         
            i=self._get_position_internal(position)
            if not i:
                raise ActiveSiteError(position,variant)

            i_mapped=self.position2index[i]
            sequence[i_mapped]=Ai

        X_var=np.zeros(sequence.size,dtype=float)
        for i,Ai in enumerate(sequence):
            Ai_index=self.alphabet2index[Ai]
            X_var[i]=self.hi[i,Ai_index] + 0.5*self.Ji(i,Ai_index,sequence)
        return X_var

    def _encode_wt(self) -> np.ndarray:
        """
        Description
        -----------
        Encodes the wild-type using its "DCA representation".
        
        Returns
        -------
        X_wt : np.ndarray
            Encoded sequence of the wild-type.
        """
        X_wt=np.zeros(self.targetSeq.size,dtype=float)
        for i,Ai in enumerate(self.targetSeq):
            Ai_index=self.alphabet2index[Ai]
            X_wt[i]=self.hi[i,Ai_index] + 0.5*self.Ji(i,Ai_index,self.targetSeq)
        return X_wt

def _get_data(variants:list, fitnesses:list, dcaEncode:object, data:list) -> list:
    """
    Description
    -----------
    Get the variant name, the associated fitness value, and its (DCA-)encoded sequence.

    Parameters
    ----------
    variants : list
        List of strings containing the variants to be encoded.
    fitnesses : list
        List of floats (1d) containing the fitness values associated to the variants.
    dcaEncode : object
        Initialized 'Encode' class object.
    data : manager.list()
        Manager.list() object to store the output of multiple processors. 

    Returns
    -------
    data : manager.list()
        Filled list with variant names, fitnesses, and encoded sequence.
    """
    for variant,fitness in zip(variants,fitnesses):
        try:
            data.append([variant,dcaEncode._encode_variant(variant),fitness])
            
        except ActiveSiteError:
            pass 

def get_data(fitnessKey:str, csvFile:str, dcaEncode:object, nProcesses=6):
    """
    Description
    -----------
    This function allows to generate the encoded sequences based on the variants
    given in 'csvFile' in a parallel manner.
    
    Parameters
    ----------
    fitnessKey : str
        Name of column containing the fitness values.
    csvFile : str
        Name of the csv file containing variant names and associated fitness values.
    dcaEncode : object
        Initialized 'Encode' class object.
    nProcesses : int
        Number of processes to be used for parallel execution (default=6).

    Returns
    -------
    data : np.ndarray
        Filled numpy array including variant names, fitnesses, and encoded sequences.
    """

    df=pd.read_csv(csvFile,sep=';',comment='#')

    fitnesses=df[fitnessKey].to_numpy()
    variants=df['mutant'].to_numpy()

    idxs_nan=np.array([i for i,b in enumerate(np.isnan(fitnesses)) if b]) # find NaNs
    if idxs_nan.size>0: # remove NaNs if presented
        print('NaNs are:', idxs_nan)
        fitnesses=np.delete(fitnesses,idxs_nan)
        variants=np.delete(variants,idxs_nan)

    fitnesses_split=np.array_split(fitnesses,nProcesses)
    variants_split=np.array_split(variants,nProcesses)

    manager=multiprocessing.Manager()
    data=manager.list()

    processes=[]
    for variants,fitnesses in zip(variants_split,fitnesses_split):
        p=multiprocessing.Process(target=_get_data, args=[variants,fitnesses,dcaEncode,data])
        p.start()
        processes.append(p)

    for p in processes:
        p.join()

    return np.array(data,dtype=object)
    
def generate_dataframe(data:np.ndarray, csvFile:str, chunksize=100):
    """
    Description
    -----------
    Takes 'csvFile' and generates new csv file containing the variants, fitnesses, and encoded sequences.
    Parameters
    ----------
    data : np.ndarray
        Filled numpy array including variant names, fitnesses, and encoded sequence.
    csvFile : str
        Name of the output csv-file.
    """
    variants,X,fitnesses=np.array(data,dtype=object).T # Can cause error if data.size==0 ?!
    X=np.stack(X)

    columns = ['variant', 'y']
    [columns.append('X%d'%(i)) for i in range(X.shape[1])]

    append = False
    for i in range(0,variants.size,chunksize):
        variantsChunked = variants[i:i+chunksize]
        fitnessesChunked = fitnesses[i:i+chunksize]
        encodedSequences = X[i:i+chunksize]

        df_data = []
        for variant, fitness, encodedSequence in zip(variantsChunked, fitnessesChunked, encodedSequences):
            df_row = [variant, fitness]
            [df_row.append(xi) for xi in encodedSequence]
            df_data.append(df_row)


        df = pd.DataFrame(data = df_data, columns=columns)

        if append:
            df.to_csv(csvFile, sep=';', index=False, mode='a', header=False)

        else:
            df.to_csv(csvFile, sep=';', index=False)
            append = True
