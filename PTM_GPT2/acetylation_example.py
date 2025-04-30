# Import necessary libraries for the model
import os
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

tokenizer_path = 'C:/Users/nikhi/Downloads/PTMGPT2-main/PTMGPT2-main/Tokenizer'

#result = None

fasta_path = 'C:/Users/nikhi/Downloads/AceMutansUnenrichNegControl.fasta'
# Path was changed according to fasta file
sequences = read_fasta(fasta_path)

# Loop through all sequences in the FASTA
for seq_id, seq in sequences.items():
    if 'K' in seq:
        result = inference(seq, tokenizer_path, 'C:/Users/nikhi/Downloads/PTMGPT2-models-Part1/Acetylation (K)/', ['K'])

        pred_ptm_sites = []
        for entry in result['Results']:
            residue = next(iter(entry))
            if entry[residue] == 'POSITIVE':
                pred_ptm_sites.append(f'{seq[residue - 1]}{residue}')

        pred_ptm_sites = ';'.join(pred_ptm_sites)
        print(f'Sequence ID: {seq_id}')
        print(f'Prediction sites: {pred_ptm_sites}')


pred_ptm_sites = list()

# If a ptm is detected, then it is added to the list pred_ptm_sites

for entry in result['Results']:
    residue = next(iter(entry))
    if entry[residue] == 'POSITIVE':
        pred_ptm_sites.append(f'{seq[residue - 1]}{residue}')

pred_ptm_sites = ';'.join(pred_ptm_sites)

print(f'Prediction sites: {pred_ptm_sites}')