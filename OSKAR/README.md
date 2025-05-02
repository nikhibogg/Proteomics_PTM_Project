
%%THIS README CORRESPONDS TO THE SCRIPTS USED TO FORMAT THE NETPHOS 3.1 AND SITETACK PREDICTION DATA FOR S PNEUMONIAE AND S MUTANS.
## THE SITETACK RELATED DATA IS IN THE SITETACK FOLDER IN THE SAME STRUCTURE IT WAS DOWNLOADED FROM GITHUB LINK https://github.com/clair-gutierrez/sitetack
## IT IS PREFERRED TO RUN SITETACK THROUGH A PYENV (CHECK REQUIREMENTS,TXT)

##DUE TO THE NATURE OF SOME OF THE NOTEBOOKS AND SCRIPTS MOST OF THE FILEPATHS ARE PROBABLY BROKEN. HOWEVER, AS I AM UNSURE OF HOW EXACTLY TO FIX THIS
## GIVEN SOME FUNCITONS REQUIRE YOU TO SPECIFY YOUR HOME DIRECTORY, I AM SIMPLY LEAVING THE FILEPATHS UNCHANGED,

#FastaSplit.py
Takes an input .txt file corresponding to a multifasta, and splits it into chuncks based on each protein ID.

#Fix_phospho_PLEASE.py
Converts an input .csv into a table usable for downstream quantification of NetPhos prediction. Takes data about residue location, phosphorylation type, and 
what protein ID it corresponds to. The it outputs them as protein ID and residuepositions columns. This data was utilized by Dakota for quantification of 
predictions.

#Combine_ptms_y.py
Converts an input .txt file into an xslx file for use downstream for training the Sitetack model.
Was used for both S,T, and Y predictions despite the name. 

#generate_target_decoy_from_fasta.py
This code was used to generate the decoy fasta for S pneumoniae seen in:/home/kleptoswirlig/CODE_PROTEOMICS/Negative Control Data/Decoy/Pneumongus_decoy.fasta

#protein_database_search.py
This code was utilized in an Ursgal (Ver 0,6.9) python environment, using the above generated decoy fasta, MsFragger (ver3.0.0) search tool, and percolator (ver 3.4.0) validation tool. It runs the search for S/T/and Y phosphorylation sites. The output combined validated.csv files were
renamed QE_S_pneumoniae_neg_with_pos_hits.csv and Sp1_pneumo_neg1_With_pos_hits.csv and placed in /home/kleptoswirlig/CODE_PROTEOMICS/Negative Control Data/SPneumo_1_2_Neg/.


#Net_csv_convert.py
Converts the output Prediction data from Sitetack into a format that can be used for downstream analysis. Filters out any predictions with a quality score < 0.5.

##Split_fastas_<Bacteria>
These folders contain batch files made using SplitFasta.py, and their corresponding .csv files, as formatted using Fix_phospho_PLEASE.py from the NetPhos Prediction. THe N_chunk_out_x.txt are the output NetPhos files, before running the Fix_phospho_PLEASE.pu function.



####Jupyter Notebooks->
##jupyter notebook /sitetack_app/train/Construct_datasets_train_models.ipynb
This function was utilized to train a dataset using filtered reads predicted by Musitedeep. This prediction data was formatted as a xslx sheet with Sheet1 
containing columns "Protein ID" and "Position" for each protein ID and the integer corresponding to the sites of Phosphorylation.
Each cell in the notebook has annotations for their respective function, aside from the cells containing exclusivly funtions ironically.
The output of this file generates a folder named 'PTM_name' into a folder in the same directory as the .ipynib folder. 
It is currently set to generate an output 10 test datasets in to a folder called "Phosphorylation (Y) Pneumo". It is currently set to only uptake data from the Y phosphorylation.

##jupyter notebook /sitetack_app/train/Predict_on_sequence_pneumo.ipynb   
Predicts the phosphorylation site based on the training model data outlined in Construct_datasets_train_models.ipynb
In order to avoid timing out the jupyter notebook, there are 13-14 cells corresponding to the 13-14 batches of S mutans and S pneumonia fastas that the prediction was run on.
The batch files for S mutans and S pneumonia proteomes is in ____ and ____ respectively.



REFERENCES:

MusiteDeep
 Wang, D., et al. (2020) MusiteDeep: a deep-learning based webserver for protein post-translational modification site prediction and visualization, Nucleic Acids Research,Volume 48, Issue W1, 02 July 2020, Pages W140–W146

Sitetack
Gutierrez CS, Kassim AA, Gutierrez BD, Raines RT. 2024. Sitetack: A deep learning model that improves PTM prediction by using known PTMs. Bioinformatics. 40(11). doi:https://doi.org/10.1093/bioinformatics/btae602.

NetPhos3.1
Sequence- and structure-based prediction of eukaryotic protein phosphorylation sites.Blom, N., Gammeltoft, S., and Brunak, S.Journal of Molecular Biology: 294(5): 1351-1362, 1999. doi:https://pubmed.ncbi.nlm.nih.gov/10600390/