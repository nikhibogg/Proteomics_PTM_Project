# This program is specifically for Appendix Table S1 of Lei et al 2020
# "Quantitative acetylome analysis reveals involvement of glucosyltransferase acetylation in Streptococcus mutans biofilm formation"
# Requires xlrd, requests and pandas packages
inputFile = "emi412907-sup-0001-tables1.xls"
inputSheet = "Annotation summary"
idColName = "Protein"
posColName = "Position"

import pandas
import requests

dataframe = pandas.read_excel(inputFile, sheet_name=inputSheet)
# Add a placeholder at the end for when the last entry looks at the next row
dataframe.loc[len(dataframe), idColName] = "XXXXX"

positionList = []
for rowNum in range(len(dataframe.index)):

    ID = dataframe.loc[rowNum, idColName]
    if ID == "XXXXX":
        break
    position = dataframe.loc[rowNum, posColName]
    positionList.append(position)
    
    nextRow = rowNum + 1
    if (dataframe.loc[nextRow, idColName] != ID):

        header = ">" + ID 
        for pos in positionList:
            header = header + "_" + str(pos)

        # ChatGPT was used to make this code for finding the sequence from uniprot id.
        url = f"https://rest.uniprot.org/uniprotkb/{ID}.fasta"
        response = requests.get(url)
    
        if response.status_code == 200:
            fasta_data = response.text
            # Extract the sequence (ignore the FASTA header)
            sequence = "".join(fasta_data.split("\n")[1:])
        else:
            print(f"Error: Unable to fetch data for UniProt ID {uniprot_id}")
            break


        positionList = []
        with open("AcPosCon.fasta", "a") as output:
            output.write(header)
            output.write("\n")
            output.write(sequence)
            output.write("\n")

    
    # If next ID not the same, print the header with ID and each position, clear position list, and search + print sequence from ID



