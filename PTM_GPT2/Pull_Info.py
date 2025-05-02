# Needed to reformat output files from PTM_GPT2
def reformat_predictions_file(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        lines = [line.strip() for line in infile if line.strip()]
        
        for i in range(0, len(lines), 2):
            seq_id = lines[i].replace("Sequence ID: ", "")
            predictions = lines[i+1].replace("Prediction sites: ", "").split(';')
            # seperarted each ptm into it's own line
            prediction_numbers = [p[1:] for p in predictions if p.startswith('K')]
                # output from ptm_gpt_2 hasa k before each site that needed to be removed

            outfile.write(seq_id + '\n')
            for number in prediction_numbers:
                outfile.write(number + '\n')

reformat_predictions_file('C:/Users/nikhi/Downloads/AceMutansUnenrichNegControl_Results', 'C:/Users/nikhi/Downloads/AceMutansUnenrichNegControl_PTM2_Results_Filtered.txt')
# writes to a new files
