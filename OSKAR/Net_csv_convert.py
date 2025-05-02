import pandas as pd


## Converts phosphorylation predictions into a format easily readable by later functions.
# Load the CSV file
df = pd.read_csv("/home/kleptoswirlig/split_fastas_pneumo/Chunk_13_Net_out_predictions_combined.csv")

# Filter rows with Score > 0.5
df_filtered = df[df['Score'] > 0.5].copy()

# Create a new column combining Residue and Position
df_filtered['ResiduePosition'] = df_filtered['Residue'] + df_filtered['Position'].astype(str)

# Group by ProteinID and join the residue-position strings with commas
result = df_filtered.groupby("ProteinID")["ResiduePosition"].apply(lambda x: ', '.join(x)).reset_index()

# Rename columns
result.columns = ["ProteinID", "ResiduePositions"]

# Save to CSV
result.to_csv("/home/kleptoswirlig/split_fastas_pneumo/Residue_Ranges_output_chunk_5_13.csv", index=False)

print("Conversion complete!")
