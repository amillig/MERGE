# version         0.1.7
# date            19.02.2024
# author          Alexander-Maurice Illig
# affilation      Chair of Biotechnology, RWTH Aachen
# email           a.illig@biotec.rwth-aachen.de

import argparse
import numpy as np
from Bio import AlignIO

parser = argparse.ArgumentParser()
parser.add_argument('-sto', help="Filename of the alignment in stockholm format to be converted into fasta format.")
parser.add_argument('-interGap', help="Fraction to delete all positions with more than 'interGap' * 100 %% gaps (columnar trimming). | default=0.3",default=0.3, type=float)
parser.add_argument('-intraGap', help="Fraction to delete all sequences with more than 'intraGap' * 100 %% gaps after being columnar trimmed (line trimming). | default=0.5",default=0.5, type=float)
args = parser.parse_args()

def convert_sto2a2m(stoFile:str, interGap:float, intraGap:float):
    # Generate the a2m output filename
    a2mFile="%s.a2m"%(stoFile.split(".sto")[0])

    # Load the stockholm alignment
    stoAlignment=AlignIO.read(stoFile, 'stockholm')

    # Save this 'raw' multiple sequence alignment as numpy array 
    rawMsa=[]
    for record in stoAlignment:
        rawMsa.append(np.array(record.seq))
    rawMsa=np.array(rawMsa)

    # 1st processing step
    # Delete all positions, where WT has a gap to obtain the 'trimmed' MSA
    ungapPos=np.where(rawMsa[0]=="-")
    msaTrimmed=np.array([np.delete(seq, ungapPos) for seq in rawMsa])

    # 2nd processing step
    # Remove ("lower") all positions with more than 'interGap'*100 % gaps (columnar trimming)
    countGaps=np.count_nonzero(msaTrimmed=='-', axis=0)/msaTrimmed.shape[0]
    lower=[idx for idx, count in enumerate(countGaps) if count > interGap]     
    msaTrimmedT=msaTrimmed.T
    for idx in lower:
        msaTrimmedT[idx]=np.char.lower(msaTrimmedT[idx])    
        # replace all columns that are "removed" due to high gap content and have an "-" element by "." 
        msaTrimmedT[idx]=np.where(msaTrimmedT[idx] == '-', '.' , msaTrimmedT[idx])
    msaTrimmedInterGap=msaTrimmedT.T

    # 3rd processing step
    # Remove all sequences with more than 'intraGap'*100 % gaps (line trimming)
    targetLen=len(msaTrimmedInterGap[0])
    gapContent=(np.count_nonzero(msaTrimmedInterGap=="-",axis=1) + np.count_nonzero(msaTrimmedInterGap==".",axis=1))/targetLen
    delete=np.where(gapContent > intraGap)[0]
    msaFinal=np.delete(msaTrimmedInterGap,delete,axis=0)
    seqsCls=[seqCls for idx, seqCls in enumerate(stoAlignment) if not idx in delete]
    chunkSize=60
    with open(a2mFile, 'w') as f:
        for seq, seqCls in zip(msaFinal, seqsCls):
            f.write('>' + seqCls.id + '\n')
            for chunk in [seq[x:x+chunkSize] for x in range(0, len(seq), chunkSize)]:
                f.write("".join(chunk) + '\n')

    # Get number of sequences and active sites in the aligment
    nSeqs=msaFinal.shape[0]
    nSites=sum(1 for char in msaFinal[0] if char.isupper())
    return nSeqs,nSites,targetLen
    

if __name__=="__main__":
    nSeqs,nActiveSites,nSites=convert_sto2a2m(args.sto, args.interGap, args.intraGap)

    print("nSeqs: %d"%(nSeqs))
    print("nActiveSites: %d (out of %d sites)"%(nActiveSites,nSites))
    print("-le: %.1f"%(0.2*(nActiveSites-1)))