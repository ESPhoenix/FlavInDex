
########################################################################################
# import libraries
import os
from os import path as p
import pandas as pd
from subprocess import run
import numpy as np  
from Bio.Blast.Applications import NcbiblastpCommandline
from io import StringIO
from Bio.Blast import NCBIXML
import multiprocessing
########################################################################################
# utilities
def ifnotmkdir(dir):
    if not p.isdir(dir):
        os.mkdir(dir)
    return dir

#################
def fasta2df(fastaFile):
    fastaList = []
    entry = []

    for line in open(fastaFile,"r"):
        line = line.strip()

        if line.startswith(">"):
            if entry:
                fastaList.append(entry)
                entry = []
            entry.append(line[1:])
        else:
            entry.append(line)

    if entry:
        fastaList.append(entry)

    df = pd.DataFrame(fastaList, columns=["ID", "Sequence"])
    return df

########################################################################################
def gen_similarity_matrix(fastaFile,outputName):
    # Parse the .clstr file and extract cluster representatives and members
    fastaDf =   fasta2df(fastaFile)
    numclusters=fastaDf.shape[0]
    # Calculate pairwise similarities between cluster representatives
    similarity_matrix = np.zeros((numclusters, numclusters))  
    pairList=[]
    for i in range(0,numclusters):
        for j in range(0,numclusters):
            if j > i:
                pairList.append([i,j])

    # Create a pool of 15 worker processes
    num_cpus = 15
    pool = multiprocessing.Pool(processes=num_cpus)
    # run blast search for all ssequence pairs
    for pair in pairList:
        print(f"Calculating similarity for proteins {pair[0]} and {pair[1]}")
        similarity=pool.apply_async(calculate_similarity,(pair,fastaDf)).get()
        similarity_matrix[pair[0]][pair[1]] = similarity
        similarity_matrix[pair[1]][pair[0]] = similarity

    # Close the pool and wait for all processes to complete
    pool.close()
    pool.join()


    # write similarity matrix to file
    similarityDf=pd.DataFrame(similarity_matrix)
    similarityDf.to_csv(f"similarity_matrix_{outputName}.csv")

########################################################################################

def calculate_similarity(sequencePair,fastaDf) :
    #extract IDs and sequences from fastaDf for pair
    ID_i=fastaDf.iloc[sequencePair[0]]["ID"]
    ID_j=fastaDf.iloc[sequencePair[1]]["ID"]
    seq_i=fastaDf.iloc[sequencePair[0]]["Sequence"]
    seq_j=fastaDf.iloc[sequencePair[1]]["Sequence"]

    # Create temporary files for input and output
    queryFile = 'query.fasta'
    subjectFile='subject.fasta'
    
    # Write sequences to the input file in FASTA format
    with open(queryFile, 'w') as q, open(subjectFile,'w') as s:
        q.write(f'>{ID_i}\n{seq_i}\n')
        s.write(f'>{ID_j}\n{seq_j}\n')
    
    # Set up the BLAST command
    blastOutput = NcbiblastpCommandline(cmd='blastp', query=queryFile, subject=subjectFile, outfmt=5)()[0]
    blastRead=NCBIXML.read(StringIO(blastOutput))
    # Retrieve the score for the overall alignment
    score = None
    for alignment in blastRead.alignments:
        for hsp in alignment.hsps:
            score = float(hsp.score)  # Retrieve the score
            return score
            break  # Exit the loop after retrieving the first alignment score
        if score is not None:
            return float("nan")
            break  # Exit the loop if a score is retrieved

    # clean up
    os.remove(queryFile)
    os.remove(subjectFile)
    

########################################################################################
def main():
    # input directory and files
    clusterDir="/home/esp/dataset_generation/flavin_dataset/021_post_blast_CD-Hit_results"
    fastaFile=p.join(clusterDir,"output_data_50.fasta")
    os.chdir(clusterDir)
    # run functions
    gen_similarity_matrix(fastaFile=fastaFile,outputName="50percent")
########################################################################################
# run main (if statement prevents running if this script is imported)
if __name__ == '__main__':
    main()