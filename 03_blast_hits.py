
## Import Basic Libraries
import os
from os import path as p
import pandas as pd
# xml handling
import xml.etree.ElementTree as ET
import csv
################################################################################################################
# # utilities
def ifnotmkdir(dir):
    if not p.isdir(dir):
        os.mkdir(dir)
    return dir

def makeFasta(enzymeId,sequence):
    fasta=f'>{enzymeId}\n{sequence}'
    return fasta

################################################################################################################
def xml2csv(xmlFile,outDir):

    # Parse the XML file
    tree = ET.parse(xmlFile)
    root = tree.getroot()   
    blockNum=p.splitext(p.basename(xmlFile))[0].split("_")[-1]
    outFile=p.join(outDir,f"blast_result_{blockNum}.csv")
    if p.isfile(outFile):
        return
    print(f"-->\t Writing {outFile}")
    # Open the CSV file for writing
    with open(outFile, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)    
        # Write the header row
        writer.writerow(["Query Accession", "Alignment Accession", "Alignment Info","Alignment Sequence", "HSP E-value", "HSP Bit Score", "HSP Identity"])   
        # Iterate over the XML data and extract the required information
        for query in root.findall(".//Iteration"):
            query_accession = query.find("Iteration_query-def").text.split()[0] 
            # for each hit, get accession and sequence
            for hit in query.findall(".//Hit"):
                alignment_accession = hit.find("Hit_accession").text
                alignment_info = hit.find("Hit_def").text.split()[0]    
                # get alignment scores
                for hsp in hit.findall(".//Hsp"):
                    hsp_evalue = hsp.find("Hsp_evalue").text
                    hsp_bitscore = hsp.find("Hsp_bit-score").text
                    hsp_identity = hsp.find("Hsp_identity").text    
                    alignment_sequence=hsp.find("Hsp_hseq").text


                    # Write the row to the CSV file
                    writer.writerow([query_accession, alignment_accession, alignment_info,alignment_sequence, hsp_evalue, hsp_bitscore, hsp_identity])
################################################################################################################
## main function
def main():
    blastResultsDir="/home/esp/dataset_generation/flavin_dataset/02_BLAST/blast_results"
    csvDir=ifnotmkdir("/home/esp/dataset_generation/flavin_dataset/02_BLAST/blast_csv_files")
    excelDir="/home/esp/dataset_generation/flavin_dataset/00_Excel_and_Fasta/03_post_BLAST"
    inputFasta="/home/esp/dataset_generation/flavin_dataset/01_initial_CD-Hit_results/data50fasta.fasta"

    #convert xml to csv
    for file in os.listdir(blastResultsDir):
        if not p.splitext(file)[1]==".xml":
              continue
        xmlFile=p.join(blastResultsDir,file)
        xml2csv(xmlFile,csvDir)
    # make a dataframe from results
    dfs=[]
    for csvFile in os.listdir(csvDir):
        csvFile=p.join(csvDir,csvFile)
        df=pd.read_csv(csvFile)
        dfs.append(df)
    blastDf=pd.concat(dfs,ignore_index=True)

    redundantDataSet=p.join(excelDir,"fully_redundant_post-BLAST.xlsx")
    if not p.isfile(redundantDataSet):
        blastDf.to_excel(redundantDataSet,index=False)
    # write fully redundant dataframe to xlsx file
    blastDf=blastDf.drop_duplicates(subset="Alignment Info")
    blastDf["FASTA"] = blastDf.apply(lambda row:makeFasta(row["Alignment Accession"],row['Alignment Sequence']),axis=1)
    # write a non-redundant dataset 
    nonRedundantDataSet=p.join(excelDir,"non_redundant_post-BLAST.xlsx")
    if not p.isfile(nonRedundantDataSet):
        blastDf.to_excel(nonRedundantDataSet,index=False)
###########################S#####################################################################################
# run main (if statement prevents running if this script is imported)
if __name__ == '__main__':
    main()