import pandas as pd
from Bio import Entrez, SeqIO

# From https://pubs.acs.org/doi/full/10.1021/pr900612v#
# For Streptococcus pneumoniae
df = pd.read_csv('pr900612v_si_001.csv', delimiter=',', index_col=False)
unique_ncbi_acc = df['NCBI_Acc#'].unique()

phosphorylation_sites = {}

for i in unique_ncbi_acc:
    pr = df.loc[df['NCBI_Acc#'] == i, 'Phosphorylated residue']
    pr = sum([r.replace(',', ';').split(';') for r in pr], [])
    phosphorylation_sites[i] = pr

with open('sp_pos_con.csv', 'w') as file:
    file.write('NCBI_Acc#,Phosphorylated residue\n')

    for i, acc in enumerate(phosphorylation_sites.keys()):
        file.write(f'{acc},')
        file.write(f'{";".join(phosphorylation_sites[acc])}')

        if i < len(phosphorylation_sites.keys()) - 1:
            file.write('\n')

    file.close()

Entrez.email = 'dmk2050@rit.edu'

sp_seqs = {}

print(f'Looking for {len(unique_ncbi_acc)} sequences...')
for i, acc in enumerate(unique_ncbi_acc):
    handle = Entrez.efetch(db='protein', id=acc, rettype='gb', retmode='text')
    record = SeqIO.read(handle, 'genbank')
    handle.close()

    #print(f"ID: {record.id}")
    #print(f"Description: {record.description}")
    #print(f"Sequence: {record.seq}")

    print(f'({i+1}/{len(unique_ncbi_acc)}) Found {acc}...')
    sp_seqs[f'{acc},{record.id}'] = record.seq

with open('sp_pos_con_seq.fasta', 'w') as file:
    for i, id in enumerate(sp_seqs.keys()):
        file.write(f'>{id}\n')
        file.write(f'{sp_seqs[id]}')

        if i < len(sp_seqs.keys()) - 1:
            file.write('\n')