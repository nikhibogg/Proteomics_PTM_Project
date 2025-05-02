# Usage: 
# python sitetack_data_transform.py <Name of input csv>

import pandas
import sys

if len(sys.argv) < 2:
    print("No input file given")
    sys.exit(0)
# Open sitetack CSV output and add a placeholder at the end
inputFilename = sys.argv[1]
dataframe = pandas.read_csv(inputFilename)
outputFilename = inputFilename
dataframe.loc[len(dataframe), "sequence name"] = "XXXXX"

# Overwrite current file and add in a header line
with open(inputFilename + ".txt", "w") as output:
    output.write("ColumnName")
    output.write("\n")
# Loop through each predicted site above 0.5 confidence and add to a list
positionList = []
for rowNum in range(len(dataframe.index)):
  #  position = dataframe.loc[rowNum, posColName]
    if dataframe.loc[rowNum, "probability"] >= 0.5:
        positionList.append(int(dataframe.loc[rowNum, "site"]))

    ID = dataframe.loc[rowNum, "sequence name"]
    if ID == "XXXXX":
        break
    # If this row is the last for this protein, add the protein and predicted positions to a CSV in the correct format for Parse_Ace_Tool_Output.py
    nextRow = rowNum + 1
    if (dataframe.loc[nextRow, "sequence name"] != ID):
        with open(inputFilename + ".txt", "a") as output:
            output.write(ID)
            output.write("\n")
            for position in positionList:
                output.write(str(position))
                output.write("\n")
        positionList = []
