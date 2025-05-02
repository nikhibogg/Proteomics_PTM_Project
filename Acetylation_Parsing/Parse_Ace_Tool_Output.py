# Requires xlrd, requests, openpyxl, and pandas packages
# Usage:
# python Parse_Ace_Tool_Output.py <input excel file>

# Settings for file containing COG classes of each protein
proteinClassFile = "AceControlClassification.xlsx"
proteinClassSheet = "AllTogether"
classIdColumn = "query"
cogCatColumn = "COG_category"

# Settings for the input excel file
positiveSheet = "Positive"
mutansNegativeSheet = "Negative_(mutans+unenrich)"
pneumoNegativeSheet = "Negative_(pneumo+enrich)"
columnsName = "ColumnName"

inputFile = "" # Keep this empty


import sys
import os
import pandas

def main():

    global cogList
    global cogDict
    cogDict = {}

    cogDF = pandas.read_excel(proteinClassFile, sheet_name=proteinClassSheet)
    # Make dictionary of protein IDs and matching COG classes
    for rowNum in range(len(cogDF.index)):
        ID = cogDF.loc[rowNum, classIdColumn]
        classification = cogDF.loc[rowNum, cogCatColumn]
        cogDict[ID] = classification
    # Make a list of COG classes in a consistent order
    cogList = cogDF[cogCatColumn].tolist()
    cogList = list(set(cogList))
    cogList.sort()

    # Read in prediction results for all three datasets and call the function to parse them, indicating whether it is positive or negative control data

    df = pandas.read_excel(inputFile, sheet_name=positiveSheet)
    infoList = df[columnsName].tolist()
    print("Parsing positive control data")
    parse(True, infoList)

    df = pandas.read_excel(inputFile, sheet_name=mutansNegativeSheet)
    infoList = df[columnsName].tolist()
    print("Parsing mutans negative control data")
    parse(False, infoList)

    df = pandas.read_excel(inputFile, sheet_name=pneumoNegativeSheet)
    infoList = df[columnsName].tolist()
    print("Parsing pneumo negative control data")
    parse(False, infoList)



def parse(positiveControl, infoList):

    # Make dataframe to store correct and incorrect info for each COG class (and total)
    dataframe = pandas.DataFrame(columns=["correct", "incorrect"], index= ["Total"] + cogList)
    with pandas.option_context('future.no_silent_downcasting', True):
        dataframe = dataframe.fillna(0)
    # Make a set to store predicted positions
    toolSet = set()
    # Add placeholder to input list
    infoList.append("Placeholder string")
    # Start with the first protein in the input list
    currentHeader = infoList[0]
    # Loop through the input list
    for entry in infoList:
        # If the next entry is a protein name, check the list of predicted positions against correct positions
        # Start by getting the COG class of the protein from the dictionary and splitting the correct answers in the protein name
        if type(entry) == str:
            if currentHeader in cogDict:
                cogType = cogDict[currentHeader]
            else:
                cogType = "-"
            correctAnswers = currentHeader.split("_")
            correctAnswers = correctAnswers[1:]
        # Loop through the correct sites from the protein name
            for correctAnswer in correctAnswers:
                # Turn the correct answer into an integer so it can be compared to the predicted sites
                correctAnswer = int(correctAnswer.split(".")[0])
                inToolSet = correctAnswer in toolSet
                # Add a tally for either correct or incorrect (for both total and cog rows) depending on whether the control site was predicted correctly as acetylated or not
                if (inToolSet and positiveControl) or (not inToolSet and not positiveControl):
                    dataframe.loc["Total", "correct"] += 1
                    dataframe.loc[cogType, "correct"] += 1
                 
                if (inToolSet and not positiveControl) or (not inToolSet and positiveControl):
                    dataframe.loc["Total", "incorrect"] += 1
                    dataframe.loc[cogType, "incorrect"] += 1
                
            # Reset the list of predicted sites and set the current protein now that the last protein is finished
            toolSet = set()
            currentHeader = entry

        # If the entry is not a header, add it to the list of predicted positions
        elif type(entry == int):
            toolSet.add(entry)


    # Print the results for this control dataset
    print(dataframe)







if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("No input provided")
        sys.exit(0)

    if not os.path.exists(sys.argv[1]):
        print("Filename doesn't exist in this folder")
        sys.exit(0)
    inputFile = sys.argv[1]
    main()
