import pandas as pd
import re

##Run on the html output of Netphos copied into a text file format.
##Parses the output of the .txt into something way more managable.
##NetPhos specificed parameter: classical output, only highest confidence score per residue, S,T, and Y predictions, not kinase specific.
# Loada the file
with open("/home/kleptoswirlig/split_fastas_mutans/chunk_14_Net_out.txt", "r") as file:
    ptm_raw_text = file.read()

# Split entries at each protein section Gives us individual fastas.
entries = ptm_raw_text.strip().split('>')[1:]

data_rows = []

for entry in entries:
    #splits each fasta into seperate lines.
    lines = entry.strip().splitlines()
    if not lines:
        continue
    # pulls out protein id from the lines ripped from each protein fasta.
    protein_id = lines[0].strip().split()[0]

    # Find header line using the worst syntax ever concieved for a table output
    header_line_index = None
    for i, line in enumerate(lines):
        if re.match(r"^Sequence\s+x\s+Context\s+Score\s+Kinase\s+Answer", line):
            header_line_index = i
            break

    if header_line_index is None:
        continue  # Skip if no header found

    # Start 1 line after header to dodge the atrocity that is the Netphos row splitter
    for line in lines[header_line_index + 1:]:
        if not line.strip():
            continue

        # Use re to split on 2 or more whitespace characters
        parts = re.split(r'\s{2,}', line.strip())

        if len(parts) >= 6:
            protein_id = parts[0]
            position_residue = parts[1].split()
            if len(position_residue) != 2:
                continue  # Skip malformed lines (avoids the ACII PTM residue visualization)
            #Segregates the input NetPhos format into seperate parts, based on space delineation, this should hopefully align with the format for their respective values.
            position, residue = position_residue
            context = parts[2]
            score = parts[3]
            kinase = parts[4]
            answer = parts[5]
            # attempts to append a functional lsit from the input fasta.
            try:
                data_rows.append([
                    protein_id,
                    int(position),
                    residue,
                    context,
                    float(score),
                    kinase,
                    answer
                ])
            except Exception as e:
                print("Skipping row due to error:", parts, e)


# Create and append DataFrame using the data appended by the try function.
columns = ["ProteinID", "Position", "Residue", "Context", "Score", "Kinase", "Answer"]
df = pd.DataFrame(data_rows, columns=columns)

# Save to CSV
df.to_csv("/home/kleptoswirlig/split_fastas_mutans/Chunk_14_Net_out_predictions_combined.csv", index=False)
print(f"Parsed {len(df)} entries.")
