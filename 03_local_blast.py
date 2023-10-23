################################################################################################################
#   --> The aim of this script is to take a large set of sequences in fasta format (eg.the result of clustering)
#       and to run a BLAST search on each sequence. The script will then collect the results in a csv file
#       for further processing
#
#   --> This script has been written with the use of docker in mind. It should be accompanied by the files:
#       "Dockerfile", "{inputData}.fasta", "get_libs.txt", and the empty directory "/blast_results"
#
#   --> You will need to supply you own {inputData}.fasta and change the variable at the top of the main():
#       function
################################################################################################################
#       For use in a miniconda / conda environment you will need to do the following:
#           --> pip install biopython
#           --> pip install bioservices
#           --> download ncbi-blast-2.14.0+-x64-linux.tar.gz 
#               from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
#           --> transfer tarball to ~/bin
#           --> tar -xvzf ncbi-blast-2.14.0+-x64-linux.tar.gz
#           --> add to .bashrc:     export PATH="/home/eugene/bin/ncbi-blast-2.14.0+/bin:$PATH"
#           --> source .bashrc
#
################################################################################################################
### Import Basic Libraries
import sys
import os
from os import path as p
import time
# Import multiprocessing for paralell CPUs
import multiprocessing
# Import commandline functions from NCBI
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
################################################################################################################
## utilities
# basic dir creation with a check
def ifnotmkdir(dir):
    if not p.isdir(dir):
        os.mkdir(dir)
    return dir
# basic file creation with a check
def ifnotmkfile(file):
    if not p.isfile(file):
        open(file,"w").close()
    return file
# split fasta file into blocks of n - makes a new fasta file for each of these blocks in ./fasta_blocks
def splitFasta(mainDir,fastaFile,blockSize):
    fastaBlockDir=ifnotmkdir(p.join(mainDir,"fasta_blocks"))
    blocks=[]

    tmpFasta=''
    with open(fastaFile,'r') as f:
        for line in f:
            if line.startswith('>') and tmpFasta=='':
                tmpFasta=line
            elif line.startswith('>'):
                blocks.append(tmpFasta)
                tmpFasta=line
            else:
                tmpFasta+=line
        blocks.append(tmpFasta)
    blockNum=0
    for block in blocks:
        blockNum+=1
        with open(p.join(fastaBlockDir,f"fasta_block_{blockNum}.fasta"),"w") as outFile:
            outFile.write(block)
    numBlocks=len(blocks)
    return fastaBlockDir, numBlocks
# converts fasta file (result of cdhit clustering) to list of fasta
def fastaFileToList(mainDir,fastaFile,blockSize):
    fastaList=[]
    fastaTmp=''
    with open(fastaFile,'r') as f:
        for line in f:
            if line.startswith(">"):
                fastaTmp=line
            else:
                fastaTmp+=line
                fastaList.append(fastaTmp)
    blocks=[]       
    numBlocks= len(fastaList) // blockSize
    remainder= len(fastaList) % blockSize
    for i in range(numBlocks):
        block = ''.join(fastaList[i * blockSize : (i + 1) * blockSize])
        blocks.append(block)
    if remainder > 0:
        block = ''.join(fastaList[numBlocks*blockSize :])
        blocks.append(block)

    fastaBlockDir=ifnotmkdir(p.join(mainDir,"fasta_blocks"))
    blockNum=0
    for block in blocks:
        blockNum+=1
        with open(p.join(fastaBlockDir,f"fasta_block_{blockNum}.fasta"),"w") as outFile:
            outFile.write(block)
    numBlocks=len(blocks)
    return fastaBlockDir, numBlocks

# basic stopwatch
def stopwatch(command,timeLogFile=None,comment=None,startTime=None):
    if command=="START":
        startTime=time.time()
        return startTime
    elif command=="STOP":
        stopTime=time.time()
        executionTime=stopTime-startTime
        # format nicely
        minutes,seconds=divmod(executionTime,60)
        hours,minutes=divmod(minutes,60)
        executionTime="%02d:%02d:%02d" % (hours, minutes, seconds)
        #write to file
        with open(timeLogFile,"a") as f:
            f.write('{} {}\n'.format(comment,executionTime))
################################################################################################################
def run_blast_local(outDir,databaseDir,fasta,count):
    # start stopwatch
    startTime=stopwatch("START",None,None,None)
    # run qBLAST search on local nr database, write result to xml file
    outFile = p.join(outDir,f"BLAST_Block_{count}.xml")
    errFile = p.join(outDir,f"Errors_Block_{count}.txt")
    print(f'Running BLAST for batch {count}')
    blastCommand=NcbiblastpCommandline(cmd="blastp",query=fasta, db=databaseDir, out=outFile,  outfmt=5)
    stdout, stderr = blastCommand()

    with open(errFile,"w") as f:
        f.write(stderr)
    print(f'Completed BLAST for batch {count}')

    # stop stopwatch, write to file
    timeLogFile=ifnotmkfile(p.join(outDir,"times_report.log"))
    comment=f'Local BLAST for batch {count}:\t'
    stopwatch("STOP",startTime,timeLogFile,comment)

################################################################################################################
# parses results of run_BLAST
def parse_BLAST_result(xmlFile):
        # initialse empty masterlist
        hitsList=[]
        # open and parse output xml file from NCBIWWW search 
        result_handle=open(xmlFile,'r')
        blast_records=NCBIXML.parse(result_handle,debug=0)     
        # extract desired info from each hit   
        for blast_record in blast_records:
            query_id=blast_record.query_id
            for alignment in blast_record.alignments:
                # new empty temp list
                hitList=[]
                for hsp in alignment.hsps:
                    alignment_length = alignment.length         
                    hsp_evalue = hsp.expect
                    hsp_bitscore = hsp.bits
                    hsp_identity = hsp.identities
                    alignment_accession=alignment.accession
                    hsp_sequence=hsp.spjct
                # collect info into the temp list
                hitList=','.join([query_id,alignment_accession,hsp_sequence,alignment_length,hsp_evalue,hsp_bitscore,hsp_identity])
            # append to master list
            hitsList.append(hitList)
        #r output master list
        return hitsList
################################################################################################################  
# main function to run
def main():
    # database directory
    database="/scratch/non_redundant_protein_database/nr_db"
    #load file
    clusterFile="data50fasta.fasta"
    # convert to list of fasta seqs
    blockSize=5
    fastaBlockDir,numBlocks=fastaFileToList(os.getcwd(),clusterFile,blockSize)
    print(f'{numBlocks} blocks of {blockSize} fasta sequences')
    # new dir for results
    blastDir=ifnotmkdir("./blast_results")
    if p.isdir(blastDir):
        print(f'Output Directory made at:\t{blastDir}')
    else:
        print("No Output Directory made: Exiting!")
        exit

    inputData=[]
    for blockNum in range(1,numBlocks):
        blockFile=p.join(fastaBlockDir,f"fasta_block_{blockNum}.fasta")
        if p.isfile(blockFile):
            inputData.append((blastDir,database,blockFile,blockNum))

    # Create a pool of 16 worker processes
    num_cpus = 15
    pool = multiprocessing.Pool(processes=num_cpus)
    # run blast search for all sequences in list - this will make xml files in blastDir
    for blockNum in range(1,numBlocks):
        blockFile=p.join(fastaBlockDir,f"fasta_block_{blockNum}.fasta")
        pool.apply_async(run_blast_local,(blastDir,database,blockFile,blockNum))
    # Close the pool and wait for all processes to complete
    pool.close()
    pool.join()
###########################S#####################################################################################

# run main (if statement prevents running if this script is imported)
if __name__ == '__main__':
    main()