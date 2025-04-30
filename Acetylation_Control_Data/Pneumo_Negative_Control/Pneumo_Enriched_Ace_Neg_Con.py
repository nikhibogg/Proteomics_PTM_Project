# Requires xlrd, requests, openpyxl, and pandas packages
# Settings for input excel file
inputFile = "modificationSpecificPeptides.xlsx"
inputSheet = "modificationSpecificPeptides"
idColName = "Proteins"
seqColName = "Sequence"
modColName = "Acetyl (K)"


import pandas
import requests
import sys
# Read in list of peptides
dataframe = pandas.read_excel(inputFile, sheet_name=inputSheet)

negList = set()
posList = set()
seqDict = {}
# Loop through each peptide
for rowNum in range(len(dataframe.index)):

    # Check if the peptide is modified
    ID = dataframe.loc[rowNum, idColName]
    ID = ID.split("|")
    ID = ID[1]
    peptideSeq = dataframe.loc[rowNum, seqColName]
    if dataframe.loc[rowNum, modColName] == 1:
        modded = True
    elif dataframe.loc[rowNum, modColName] == 0:
        modded = False
    else:
        print("Error: Not 0 or 1 in mod column")
        sys.exit(0)
    # Search uniprot for each protein sequence and store in a dictionary
    if ID not in seqDict:
        # ChatGPT was used to make this code for finding the sequence from uniprot id.
        url = f"https://rest.uniprot.org/uniprotkb/{ID}.fasta"
        response = requests.get(url)

        if response.status_code == 200:
            fasta_data = response.text
            # Extract the sequence (ignore the FASTA header)
            sequence = "".join(fasta_data.split("\n")[1:])
            seqDict[ID] = sequence
        else:
            print(f"Error: Unable to fetch data for UniProt ID {ID}")
            sys.exit(0)
    proteinSeq = seqDict[ID]
        
    # Determine where in the protein the peptide is
    peptideStart = proteinSeq.find(peptideSeq)
    if peptideStart == -1:
        print("Didn't find peptide " + peptideSeq + " in " + ID + ", will discard this peptide")
        continue
    

    # For each K in the peptide, determine the location in the full protein and store in either positive or negative list
    proteinPos = peptideStart
    for AA in peptideSeq:
        proteinPos = proteinPos + 1
        if AA != "K":
            continue
        toAppend = ID + "-" + str(proteinPos)
        if modded:
            posList.add(toAppend)
        else:
            negList.add(toAppend)

     
# If something is in positive list, remove from negative list
finalList = negList - posList
finalList = list(finalList)

# Use the dictionary of protein sequences to add each protein and its identified negative sites to a FASTA
finalList.sort()
finalList.append("XXXX-XXXX")
positionList=[]
for index, element in enumerate(finalList):
    if element == "XXXX-XXXX":
        break
    splitElement = element.split("-")
    ID = splitElement[0]
    position = splitElement[1]
    positionList.append(position)
    nextElement = finalList[index+1]
    nextID = nextElement.split("-")[0]
    if ID != nextID:
        sequence = seqDict[ID]
        header = ">" + ID 
        for pos in positionList:
            header = header + "_" + str(pos)
        positionList = []
        with open("AcePneumoEnrichNegControl.fasta", "a") as output:
            output.write(header)
            output.write("\n")
            output.write(sequence)
            output.write("\n")






