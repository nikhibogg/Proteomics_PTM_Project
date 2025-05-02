# Import necessary libraries for the model
import os
import concurrent.futures
from Bio import SeqIO
from itertools import islice
import torch
import evaluate
import numpy as np
import pandas as pd
from torch.utils.data import DataLoader,Dataset
from transformers import DataCollatorForSeq2Seq
from transformers import AutoTokenizer, GPT2LMHeadModel,TrainingArguments, Trainer,GPT2Config
from sklearn.metrics import average_precision_score,matthews_corrcoef,f1_score, precision_score, recall_score, balanced_accuracy_score

def find_subsequences(sequence:str, chars:list, left=10, right=10):
    subsequences = []
    length = len(sequence)
    # Iterate through the sequence to find the character
    for i, c in enumerate(sequence):
        if c in chars:
            # Calculate the start and end indices for the subsequence
            start = max(0, i - left)  # Ensure start is not less than 0
            end = min(length, i + right + 1)  # Ensure end does not exceed the sequence length
            
            # Append the subsequence to the list
            subsequences.append({'Seq':sequence[start:end],
                                 'Pos':i+1,
                                 'text':f'<startoftext>SEQUENCE:{sequence[start:end]}\nLABEL:'
                                 })
    return subsequences

def read_fasta(file_path):
    """
    Reads a FASTA file and returns a dictionary with sequence identifiers as keys
    and sequences as values.

    :param file_path: str, path to the FASTA file
    :return: dict, dictionary with sequence IDs as keys and sequences as values
    """
    sequences = {}
    sequence_id = None
    sequence_data = []

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if sequence_id is not None:
                    sequences[sequence_id] = ''.join(sequence_data)
                sequence_id = line[1:]
                sequence_data = []
            else:
                sequence_data.append(line)
        
        # Add the last sequence
        if sequence_id is not None:
            sequences[sequence_id] = ''.join(sequence_data)

    return sequences

def load_model(mdl_pth):
    """
    Loads a pre-trained GPT-2 model from the specified path.

    :param mdl_pth: str, path to the model directory.
    :return: GPT2LMHeadModel, the loaded GPT-2 model in evaluation mode on the CPU.
    """
    model_config = GPT2Config.from_pretrained(mdl_pth)
    model = GPT2LMHeadModel.from_pretrained(mdl_pth,config=model_config,ignore_mismatched_sizes=True)
    return model.cpu().eval()

def tokenize(sub_sequences,tokenizer):
    """
    Tokenizes the given subsequences using the specified tokenizer.

    :param sub_sequences: list of dicts, each containing a 'text' field with the subsequence to tokenize.
    :param tokenizer: AutoTokenizer, the tokenizer to use for tokenizing the subsequences.
    :return: dict, the tokenized subsequences with padding applied.
    """
    sub_sequences=[x['text'] for x in sub_sequences]
    encoded=tokenizer(sub_sequences,return_tensors='pt',padding='longest')
    return encoded

def inference(input_seq,tokenizer_pth,model_pth,chars:list):
    """
    Performs inference on the input sequence using a specified tokenizer and model, and extracts labels.

    :param input_seq: str, the input sequence to process.
    :param tokenizer_pth: str, path to the tokenizer directory.
    :param model_pth: str, path to the model directory.
    :param chars: list of str, characters to find subsequences for.
    :return: dict, a JSON-like dictionary containing the input sequence, model type, and labeled results.
    """
    tokenizer = AutoTokenizer.from_pretrained(tokenizer_pth,padding_side='left')
    model: GPT2LMHeadModel =load_model(model_pth)
    sub_sequences=find_subsequences(input_seq,chars=chars)
    inputs_encode=tokenize(sub_sequences=sub_sequences, tokenizer=tokenizer)

    predicted = model.generate(input_ids=inputs_encode['input_ids'], attention_mask=inputs_encode['attention_mask'], pad_token_id=50259, do_sample=True, temperature=0.1, top_k=50, max_new_tokens=2, top_p=0.15)
    
    predicted_text=tokenizer.batch_decode(predicted,skip_special_tokens=True)
    predicted_labels=[x.split('LABEL:')[-1] for x in predicted_text]
    json_results={'Sequence':input_seq,
                'Type':model_pth,
                'Results':[]
                }
    for label,sub_seq in zip(predicted_labels,sub_sequences):
        json_results['Results'].append({sub_seq['Pos']:label})
    return json_results

# Function to predict phosphorylation sites for each sequence
def predict_phosphorylation(fasta, i, num_sequences, tokenizer_path):
    id = fasta.id
    seq = fasta.seq

    print(f'Predicting Phosphorylation sites for {id} ({i + 1}/{num_sequences})')

    result = None

    if 'S' in seq or 'T' in seq:
        result_st = inference(seq, tokenizer_path, 'PTMGPT2/Phosphorylation (S,T)/', ['S', 'T'])
        result = result_st

    if 'Y' in seq:
        result_y = inference(seq, tokenizer_path, 'PTMGPT2/Phosphorylation (Y)/', ['Y'])

        if result != None:
            result['Results'] = result['Results'] + result_y['Results']
            result['Type'] = 'Phosphorylation (S,T,Y)'

    pred_ptm_sites = list()

    if result != None:
        for entry in result['Results']:
            residue = next(iter(entry))
            if entry[residue] == 'POSITIVE':
                pred_ptm_sites.append(f'{seq[residue - 1]}{residue}')

    if len(pred_ptm_sites) == 0:
        pred_ptm_sites = ''
    else:
        pred_ptm_sites = ';'.join(pred_ptm_sites)
    
    return (id, pred_ptm_sites)


def main(org_name, fasta_file, output_path):
    print(f'Predicting {org_name} Phosphorylation sites...')

    tokenizer_path = 'PTMGPT2/Tokenizer/'

    fasta_sequences = SeqIO.parse(open(fasta_file), 'fasta')

    # Count the number of sequences in the FASTA file
    fh = open(fasta_file)
    num_sequences = 0
    for line in fh:
        if line.startswith(">"):
            num_sequences += 1
    fh.close()

    # Open the output file
    file = open(output_path, 'w')
    file.write('ID,Phosphorylated residue\n')
    file.flush()
    os.fsync(file.fileno())

    # Use ThreadPoolExecutor to parallelize the sequence processing
    with concurrent.futures.ProcessPoolExecutor(max_workers=5) as executor:
        # Submit all the tasks to the thread pool
        futures = []
        for i, fasta in enumerate(fasta_sequences):
            futures.append(executor.submit(predict_phosphorylation, fasta, i, num_sequences, tokenizer_path))

        # Collect and write results as they complete
        for future in concurrent.futures.as_completed(futures):
            id, pred_ptm_sites = future.result()
            file.write(f'{id},{pred_ptm_sites}\n')
            file.flush()
            os.fsync(file.fileno())

    file.close()

if __name__ == "__main__":
    main('SM', './proteome_sm.fasta', './predictions/ptmgpt2_sm.csv')
    main('SP', './proteome_sp.fasta', './predictions/ptmgpt2_sp.csv')