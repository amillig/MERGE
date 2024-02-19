# version         v0.1.7
# date            30.01.2024
# author          Alexander-Maurice Illig
# affilation      Institute of Biotechnology, RWTH Aachen
# email           a.illig@biotec.rwth-aachen.de

import pickle
import numpy as np
from ._errors import InvalidVariantError

def is_valid_substitution(substitution: str) -> bool:
    """
    Description
    -----------
    A substitution has to follow the scheme:
    First character: (wild-type/substituted) amino acid in one letter code representation
    Last character: (introduced) amino acid in one letter code representation
    In between: position (of substitution)

    If the entered substitution does not follow this scheme (integer enclosed by two one
    letter code representations of amino acids) return False, else return True.

    Parameters
    -----------
    substitution : str
        Substitution as string: Integer enclosed by two letters representing
        the wild-type (first) and variant amino acid (last) in one letter code.

    Returns
    -------
    boolian
    """

    """
    Valid characters for one letter codes of amino acids. 
    """
    aminoAcidsOlc=[
    'A','C','D','E','F',
    'G','H','I','K','L',
    'M','N','P','Q','R',
    'S','T','V','W','Y'
    ]

    if not substitution[0] in aminoAcidsOlc:
        return False

    if not substitution[-1] in aminoAcidsOlc:
        return False

    try:
        int(substitution[1:-1])
    except ValueError:
        return False

    return True

def is_valid_variant(variant:str, separator=',') -> bool:
    """
    Description
    -----------
    Gets the single substitutions of the variant and checks if they follow the required scheme.

    If the entered substitution does not follow this scheme (integer enclosed by two one
    letter code representations of amino acids) return False, else return True.

    Parameters
    ----------
    variant : str
        Joined string of integers enclosed by two letters representing the wild type
        and variant amino acid in the single letter code. -> Check separator
    separator : str
        Character to split the variant to obtain the single substitutions (default=',').
    
    Returns
    -------
    boolian
    """
    for substitution in variant.split(separator):
        if not is_valid_substitution(substitution):
            return False

    return True

def get_single_substitutions(variant:str, separator=',') -> object:
    """
    Description
    -----------
    Generator that extracts and returns the single substitutions of the entered variant.

    Parameters
    ----------
    See 'is_valid_variant' for an explanation.

    Returns
    -------
    Generator object
    """
    if is_valid_variant(variant, separator):
        for substitution in variant.split(separator):
            yield substitution

    else:
        raise InvalidVariantError(variant)

def X_to_deltaE(x:np.ndarray, xWt:np.ndarray) -> float:
    """
    Description
    -----------
    Converts the encoded variant sequence to its deltaE value.

    Parameters
    ----------
    x : np.ndarray
        Encoded variant sequence.

    xWt : np.ndarray
        Encoded wild type sequence.

    Returns
    -------
    DeltaE of between variant and wild type.
    
    """
    return np.sum(np.subtract(x, xWt), axis=1)

def save_pickle(filename:str, data):
    """
    Description
    -----------
    Function to pickle 'data' in 'filename.pkl'.

    Parameters
    ----------
    filename : str
        Name of the output pkl file.
    data :
        Data to pickle. 
    """
    with open(filename, 'wb') as f:
        pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)

def load_pickle(filename:str):
    """
    Description
    -----------
    Function to load pkl file.

    Parameters
    ----------
    filename : str
        Name of the pkl file.
    """
    with open(filename, 'rb') as f:
        return pickle.load(f)