import pandas as pd

#Used to convert MusitePrediction output for teaching Sitetack.
# Input and output filenames
input_file = '/home/kleptoswirlig/Downloads/S mutans pred Y/prediction_results.txt'  # Replace with your actual .txt file name
output_file = '/home/kleptoswirlig/Downloads/S mutans pred Y/prediction_results_Y.xlsx'

# list of final data
data = []

with open(input_file, 'r') as f:
    for line in f:
        line = line.strip()
        # Skip header and description lines
        if not line or line.startswith('ID') or line.startswith('>'):
            continue
        # Split the line by tabs
        parts = line.split('\t')
        if len(parts) >= 3:
            full_id = parts[0]
            position = parts[1]
            residue = parts[2]

            # Extract protein ID from the full_id (e.g., tr|A0A0H2ZLH7|...)
            protein_id = full_id.split('|')[1] if '|' in full_id else full_id

            data.append([protein_id, residue, int(position)])

# Create DataFrame
df = pd.DataFrame(data, columns=['Protein ID', 'Residue', 'Position'])

# Save to Excel
df.to_excel(output_file, index=False)

print(f"Saved to {output_file}")
