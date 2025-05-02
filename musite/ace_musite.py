# Usage: 
# python ace_musite.py <Name of input fasta>

import requests
import json
import sys
from Bio import SeqIO
import os

if len(sys.argv) == 0:
    print("No input provided")
    sys.exit(0)

if not os.path.exists(sys.argv[1]):
    print("Filename doesn't exist in this folder")
    sys.exit(0)


outputFile = "musite_" + sys.argv[1] + ".csv"
# For each protein in the input FASTA, check if the sequence is too long for musite
n = 0
for record in SeqIO.parse(sys.argv[1], "fasta"):
    header = record.name
   # print(header)
    sequence = record.seq
    if len(record) > 1000:
        n = n + 1
        print("Too long: skipped seq #" + str(n))
        continue
    # Send the sequence to musite to search for lysine acetylation sites
    url = "https://api.musite.net/musitedeep/N6-acetyllysine/"+sequence
    myResponse = requests.get(url)
    if(myResponse.ok):
        jData = json.loads(myResponse.content.decode('utf-8'))
        if "Error" in jData.keys(): 
            print(jData["Error"])

    else:
        myResponse.raise_for_status()
        sys.exit(0)
    # Parse through the musite results and add each position with a score above 0.5 to a list
    positionList = []
    for result in jData['Results']:
        position = result['Position']
        ptmScore = result['PTMscores']
        ptmScore = float(ptmScore[16:])
        if ptmScore >= 0.5:
            positionList.append(position)
    # Add the protein header and identified sites to a CSV in the correct format for Parse_Ace_Tool_Output.py
    with open(outputFile, "a") as output:
        output.write(header)
        output.write("\n")
        for position in positionList:
            output.write(position)
            output.write("\n")
    n = n + 1
    print("Finished with seq #" + str(n))

