READ ME

This file contains information for all the items needed to run the various programs, results from running the program and any data analysis done.

Acetylation:

	Positive control data:
	- Provided is the supplementary table 1 from Lei et al (2020) listing the MS determined KACE sites, a script, and the output
	- Install xlrd, requests, and pandas then type "python AcPosCon.py"
	- The output will be a fasta file containing each protein with KACE sites, and the header will list what residues those sites are

	Negative control data (mutans unenriched):
	- Included is the generate_target_decoy_from_fasta.py script unchanged from the one provided in the homework assignment, along with the
	theoretical proteome needed to run the script and the crapDatabase.fasta output used in the database search (also unchanged from the HW)
	- Install ursgal and run "get_file_from_ftp_server.py" to obtain the needed .raw file
	- Use convertRaw.py (unchanged from the script provided in HW) to convert to .mzml
	- Place the .mzml file in a folder within a folder containing the DB_Search script and crapDatabase.fasta
	- Type "python Mutans_Unenriched_Ace_DB_Search.py <mzml folder> trypsin crapDatabase.fasta" to perform a database search that 
	includes lysine acetylation
	- Provided is a database search output edited to remove entries as noted in the report, along with changing the protein ID column to
	only have 6 character IDs and adding "none" to blank cells in the modifications column
	- Place the edited database search output, the "emi412907-sup-0001-tables1.xls" file (same as used for positive control), and the
	Mutans_Unenrich_Ace_Neg_Con.py script in the same folder
	- Install pandas and requests then run the script to generate a fasta containing proteins with true negative sites, with the
	headers listing the true negative residues.

	Negative control data (pneumo enriched):
	- Provided is a MS database search file from PXD007687 (Liu et al 2018) with certain entries removed as noted in the report
	- Install xlrd, requests, openpyxl, and pandas, then type "python Pneumo_Enriched_Ace_Neg_Con.py"
	- The output will be a fasta file containg each protein with true negative sites, with the header listing the true negative residues

	MusiteDeep:
	- Only non-standard required librarys are Biopython and requests
	- Run with "python ace_musite.py <input fasta>
	- The program will use the MusiteDeep web server to predict lysine acetylation in each protein of the input file that has 1000 AA or less
	with a confidence level of 0.5
	- The results will be output as a csv file with one column that lists the input header followed by a list of predicted residues.
	- Provided in the folder are the three Ace control datasets and their output when run with the script
	
	Sitetack:
	- After removing rows with low confidence levels, sitetack_data_transform.py will turn the CSV files output
	by the Sitetack webserver into a format usable by the "Parse_Ace_Tool_Output.py" script (in the Acetylation_Parsing folder)
	- Provided are the three files of Sitetack output ready to be run through the script
	- Simply type "python sitetack_data_transform.py <csv file>" to reformat the data
	
	to run GPS-PAIL, just make sure to install the program, that is all that is needed to get acetylation ptm sites
	GPS PAIL - Nikhi
	- Contains program
	
	For PTM GPT 2, make sure to install: torch, evaluate, numpy, pandas, scikit-learn, transformers
	PTM_GPT_2_Acetylation - Nikhi
	- Contains script required to run it - scripts mainly modified Dakota, modified by Nikhi a little bit to make it work for Acteylation
	- Output
	- Filtered output to reformat into excel
	- Pull_Info script to format the output from PTM GPT 2
	- excel file with results from PTM GPT 2
	- not in folder but needed to run
		- fasta files

	Counting true/false positives/negatives for acetylation data (Acetylation_Parsing folder):
	- Data from the tools must be in the exact format as the excel files in the folder, with sheet and column names
	as defined in the script and the control positive/negative sites in the protein names.
	- Install xlrd, requests, pandas, and openpyxl
	- Place the script and AceControlClassification file in the same folder
	- Type "python Parse_Ace_Tool_Output.py <input excel file>"
	- On the command line the program will print results for each of the three datasets for each COG category (and total)
	- For the positive dataset, "correct" is true positive and "incorrect" is false negative
	- For the negative datasets, "correct" is true negative and "incorrect" is false positive
