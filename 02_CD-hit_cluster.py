########################################################################################
# import libraries
import os
from os import path as p
import pandas as pd
from subprocess import run
########################################################################################
# utilities
def ifnotmkdir(dir):
    if not p.isdir(dir):
        os.mkdir(dir)
    return dir

# converts sequence and id to fasta
def to_fasta(sequence, identifier):
    return f">{identifier}\n{sequence}\n"
########################################################################################
# manual inputs
def inputs():
    excelDir='/home/esp/dataset_generation/flavin_dataset/00_Excel_and_Fasta/03_post_BLAST'
    clusterDir='/home/esp/dataset_generation/flavin_dataset/03_post_blast_Clustering'
    data=pd.read_csv(p.join(excelDir,'Nr_dataset_inputs_included.csv'))

    return excelDir, clusterDir, data
########################################################################################
# runs cdhit from a given list of fasta sequences, returns a list of centriod sequences 
def run_clustering(fastaFile,tolerance,clusterDir):
    # select correct wordlength
    if 0.4 <= tolerance <= 0.5:
        wordSize = "2"
    elif 0.5 < tolerance <= 0.6:
        wordSize = "3"
    elif 0.6 < tolerance <= 0.7:
        wordSize = "4"
    elif 0.7 < tolerance <= 1.0:
        wordSize = "5"

    # make input/ output names
    tolPercent=str(int(tolerance*100))
    outputFile=p.join(clusterDir,f'output_data_{tolPercent}.fasta')

    # Run CD-HIT
    run(["cd-hit", "-i", fastaFile, "-o", outputFile, "-c", str(tolerance), "-n", wordSize])

    # Read the filtered sequences and return a list of fasta sequences
    centroidFastas = []
    with open(outputFile, "r") as f:
        thisFasta=''
        for line in f:
            if line.startswith(">"):
                thisFasta=line
            else:
                thisFasta+=line
                centroidFastas.append(thisFasta)
    return centroidFastas

########################################################################################
def main():
    # get input directories and dataframe of proteins
    excelDir, clusterDir, dataDf = inputs()
    print(dataDf["FASTA"])
    # write a fasta file from teh FASTA column in dataframe
    inputFastaFile=p.join(clusterDir,f'input_data.fasta')
    fastaSequences=dataDf["FASTA"].to_list()

    with open(inputFastaFile,"w") as outFile:
        for fasta in fastaSequences:
            outFile.write(fasta+'\n')

    # loop through sequence similarity tolerances
    for tol in [0.5,0.6,0.7,0.8,0.9]:
        # extract FASTA from dataframe to a list
        fastaSequences=dataDf["FASTA"].to_list()
        # run clustering
        centroidFastas=run_clustering(fastaFile=inputFastaFile, tolerance=tol,
                                         clusterDir=clusterDir)
        # filter original dataframe using result of clustering
        mask=dataDf["FASTA"].isin(centroidFastas)
        centroidDf=dataDf[mask]
        # write to excel file
        tolPercent=str(int(tol*100))
        outputExcel=p.join(excelDir,f'data_{tolPercent}%_seq_similarity.xlsx')
        centroidDf.to_excel(outputExcel)

########################################################################################
# run main (if statement prevents running if this script is imported)
if __name__ == '__main__':
    main()