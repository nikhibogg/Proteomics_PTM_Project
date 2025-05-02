from Bio import SeqIO
import os
import math
##Used to split input fastas into batches for easier input into NetPhos 3.1 webserver. (I.e. what inputtting data into a website with a depreciated API does to someone)
# Defines a function with easily modifiable inputs, and should be splitting the files based on the last sequence seperated by length of input fasta.

def split_fasta_by_count(input_fasta, output_folder, chunk_size=150):
    # Make output directory if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)
    
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    total_sequences = len(sequences)
    num_chunks = math.ceil(total_sequences / chunk_size)

    for i in range(num_chunks):
        start = i * chunk_size
        end = start + chunk_size
        chunk_seqs = sequences[start:end]
        output_file = os.path.join(output_folder, f"chunk_{i+1}.fasta")
        SeqIO.write(chunk_seqs, output_file, "fasta")
        print(f"✅ Wrote {len(chunk_seqs)} sequences to {output_file}")

# Example usage
# Example usage
split_fasta_by_count("/home/kleptoswirlig/Downloads/S mutans Pred/Streptococcus_mutans_UP000002512_2025_01_16.fasta", "split_fastas_mutans", chunk_size=150)