# import libraries
import os
from os import path as p
from subprocess import call

## 
# utilities
def ifnotmkdir(dir):
    if not p.isdir(dir):
        os.mkdir(dir)
    return dir
##
def split_fasta_to_dirs(fastaFile,inputDir):
    fastaList=[]
    with open(fastaFile,"r") as f:
        for line in f:
            if line.startswith(">"):
                enzName=line.strip("\n")[1:]
                accessionLine=line
            else:
                sequenceLine=line
                fastaDir=ifnotmkdir(p.join(inputDir,enzName))
                outputFasta=p.join(fastaDir,f"{enzName}.fasta")
                with open(outputFasta,"w") as outFasta:
                    outFasta.write(accessionLine+sequenceLine)
########################################################################################
def inputs():
    inputFasta="/home/eugene/AlphaFold/known_photoenzymes_for_folding.fasta"
    inputDir="/home/eugene/AlphaFold/fasta_inputs"
    predictionsDir="/scratch/photoenzymes_alphafold_predictions"

    return inputFasta, inputDir, predictionsDir
########################################################################################
def main():
    inputFasta, inputDir, predictionsDir = inputs()
    inputDir=ifnotmkdir(inputDir)
    predictionsDir=ifnotmkdir(predictionsDir)
    split_fasta_to_dirs(inputFasta,inputDir)

    fastaDirs=os.listdir(inputDir)
    for enzymeName in fastaDirs:
         # check output dir
        if p.isdir(p.join(predictionsDir,enzymeName)):
            continue

        fastaDir=p.join(inputDir,enzymeName)
        outputDir=ifnotmkdir(p.join(predictionsDir,enzymeName))


        command= ["colabfold_batch",fastaDir,outputDir]
        call(command)



########################################################################################
# run main (if statement prevents running if this script is imported)
if __name__ == '__main__':
    main()

