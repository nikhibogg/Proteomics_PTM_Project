# Requires xlrd, requests, openpyxl, and pandas packages
# Settings for input excel files
negInputFile = "WT1_msfragger_3_0_pmap_unified_percolator_3_4_0_validated_accepted.xlsx"
negInputSheet = "WT1_msfragger_3_0_pmap_unified_"
negIdColName = "Protein ID"
negSeqColName = "Sequence"
negModColName = "Modifications"
negStartColName = "Sequence Start"
negStopColName = "Sequence Stop"

posInputFile = "emi412907-sup-0001-tables1.xls"
posInputSheet = "Annotation summary"
posIdColName = "Protein"
posPositionColName = "Position"


import pandas
import requests
import sys
# Read in database search output that will be parsed
dataframe = pandas.read_excel(negInputFile, sheet_name=negInputSheet)

posList = set()
negList = set()
# Loop through each entry in the DB search
for rowNum in range(len(dataframe.index)):
    modCell = dataframe.loc[rowNum, negModColName]
    peptideSeq = dataframe.loc[rowNum, negSeqColName]
    ID = dataframe.loc[rowNum, negIdColName]
    peptideStart = dataframe.loc[rowNum, negStartColName]
    acetylPresent = modCell.find("Acetyl")
    # If the peptide is not acetylated add each lysine position to the negative set
    if acetylPresent == -1:
        proteinPos = peptideStart - 1
        for AA in peptideSeq:
            proteinPos = proteinPos + 1
            if AA != "K":
                continue
            toAppend = ID + "-" + str(proteinPos)
            negList.add(toAppend)
    # If the peptide is acetylated add the acetylated positions to the positive set
    else:
        modList = modCell.split(";")
        for modification in modList:
            if modification[:6] == "Acetyl":
                proteinPos = peptideStart - 1 + int(modification.split(":")[1])
                toAppend = ID + "-" + str(proteinPos)
                posList.add(toAppend)
# Add all the acetylated positions from the positive control data to the positive list
dataframe = pandas.read_excel(posInputFile, sheet_name=posInputSheet)
for rowNum in range(len(dataframe.index)):
    ID = dataframe.loc[rowNum, posIdColName]
    position = dataframe.loc[rowNum, posPositionColName]
    toAppend = ID + "-" + str(position)
    posList.add(toAppend)

# Remove everything in the positive set from the negative set to get the final list of negative sites
finalList = negList - posList
finalList = list(finalList)
print("Finished parsing, now creating FASTA")
# Sort the list so every position of each protein is together
finalList.sort()
# Placeholder for when the final list element looks for the next element
finalList.append("XXXX-XXXX")
# Loop through all the positions, checking if each position is the final one listed for a protein
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
    # After collecting all the positions for a protein, search uniprot for the sequence
    if ID != nextID:

        url = f"https://rest.uniprot.org/uniprotkb/{ID}.fasta"
        response = requests.get(url)
    
        if response.status_code == 200:
            fasta_data = response.text
            # Extract the sequence (ignore the FASTA header)
            sequence = "".join(fasta_data.split("\n")[1:])
        else:
            print(f"Error: Unable to fetch data for UniProt ID {ID}")
            sys.exit(0)
    # Add the protein with all the listed negative positions to a fasta
        header = ">" + ID 
        for pos in positionList:
            header = header + "_" + str(pos)
        positionList = []
        with open("AceMutansUnenrichNegControl.fasta", "a") as output:
            output.write(header)
            output.write("\n")
            output.write(sequence)
            output.write("\n")

        

