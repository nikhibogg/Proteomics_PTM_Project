import pandas as pd
import requests
import os
import numpy as np

# MusiteDeep threshold of 0.5

# for all 4 prediction files
#   get how many times those residues appear in the proteome
# make graphs

def read_proteome_fasta(path):
    prot = {}
    id = ''
    with open(path, 'r') as file:
        for i, line in enumerate(file.readlines()):
            line = line.strip()
            if line[0] == '>':
                id = line.split('|')[1]
                prot[id] = ''
            else:
                prot[id] += line
    return prot

def musite_site_breakdown(df, cutoff, proteome_path):
    raw_residue_list = df['Phosphorylated residue'].values
    residues = {'S': 0, 'T': 0, 'Y': 0}
    sites_per_protein = {}

    prot = read_proteome_fasta(proteome_path)

    for i, x in enumerate(raw_residue_list):
        protein_id = df.iloc[i, 0].split('|')[1].strip()

        if x == 'SEQ>1000':
            sites_per_protein[protein_id] = (-1, -1)
        else:
            sites = 0
            for y in x.split(';'):
                z = y.split(':')
                r = z[0][0]

                score = float(z[1])

                if score >= cutoff:
                    residues[r] += 1
                    sites += 1
            sites_per_protein[protein_id] = (sites, sites / len(prot[protein_id]))

    return [residues, sites_per_protein]

def ptmgpt2_site_breakdown(df, proteome_path):
    prot = read_proteome_fasta(proteome_path)

    residues = {'S': 0, 'T': 0, 'Y': 0}
    sites_per_protein = {}

    for i, x in df.iterrows():
        protein_id = x['ID'].split('|')[1]
        raw_residues = x['Phosphorylated residue'].split(';')

        if raw_residues == ['NA']:
            sites_per_protein[protein_id] = (0, 0)
        else:
            sites = 0
            for y in raw_residues:
                z = y[0]
                residues[z] += 1
                sites += 1
            sites_per_protein[protein_id] = (sites, sites / len(prot[protein_id]))

    return [residues, sites_per_protein]

def sitetack_site_breakdown(df, cutoff, proteome_path):
    raw_residue_list = df['Phosphorylated residue'].values
    residues = {'S': 0, 'T': 0, 'Y': 0}
    sites_per_protein = {}

    prot = read_proteome_fasta(proteome_path)

    for i, x in enumerate(raw_residue_list):
        protein_id = df.iloc[i, 0].split('|')[1].strip()

        sites = 0
        for y in x.split(';'):
            z = y.split(':')
            if z[0] == '':
                break
            r = z[0][0]

            score = float(z[1])

            if score >= cutoff:
                residues[r] += 1
                sites += 1
        sites_per_protein[protein_id] = (sites, sites / len(prot[protein_id]))

    return [residues, sites_per_protein]

def map_gi(df):
    results = {}
    
    for i in df:
        gi = i
        genbank_accession = df[i][0]
        # have ABJ53652.1, need to get associated ID in musite_pneumoniae
        url = f"https://rest.uniprot.org/uniprotkb/search?query={genbank_accession}&format=json&fields=accession,id,protein_name"

        response = requests.get(url)
        data = response.json()
        result = data.get("results", [])[0]

        uniprot_accession = result.get("primaryAccession")
        uniprot_entry_id = result.get("uniProtkbId")
        protein_name = result.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value")

        print("UniProt Accession:", uniprot_accession)
        print("UniProt Entry ID:", uniprot_entry_id)
        print("Protein name:", protein_name)

        results[uniprot_accession] = list()
        results[uniprot_accession].append(gi)
        results[uniprot_accession].append(genbank_accession)
        results[uniprot_accession].append(uniprot_accession)
        results[uniprot_accession].append(uniprot_entry_id)
    return results

def map_uniprot_accession(ids):
    results = {}
    for uniprot_accession in ids:
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_accession}.json"

        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            uniprot_entry_id = data.get('uniProtkbId')

            # Default to empty strings
            genbank_accession = ''

            # Search cross-references for GenBank/EMBL info
            xrefs = data.get("uniProtKBCrossReferences", [])
            for ref in xrefs:
                genbank_accession = ''
                for mref in ref.get('properties'):
                    if mref['key'] == 'ProteinId':
                        genbank_accession = mref['value']
                        break
                if genbank_accession != '':
                    break

            print("UniProt Accession:", uniprot_accession)
            print("UniProt Entry ID:", uniprot_entry_id)
            print("GenBank Accession:", genbank_accession)

            results[uniprot_accession] = list()
            results[uniprot_accession].append('') # unknown gi identifier
            results[uniprot_accession].append(genbank_accession)
            results[uniprot_accession].append(uniprot_accession)
            results[uniprot_accession].append(uniprot_entry_id)

    return results

def calculate_statistics(map, tool, proteome, pred_df, pos_df):
    site_counts = {'S': 0, 'T': 0, 'Y': 0} # The total number of phosphorylatable sites in s. pneumoniae
    true_positive = {'S': 0, 'T': 0, 'Y': 0} # Resides that were correctly predicted to be phosphorylated
    true_negative = {'S': 0, 'T': 0, 'Y': 0} # Resides that were correctly predicted to NOT be phosphorylated
    false_positive = {'S': 0, 'T': 0, 'Y': 0} # Resides that were incorrectly predicted to be phosphorylated
    false_negative = {'S': 0, 'T': 0, 'Y': 0} # Resides that were incorrectly predicted to NOT be phosphorylated

    for gene in proteome:
        seq = proteome[gene]
        # if map is an empty array, the protein has no phosphorylation sites
        gene_map = map[map['UniProt Accession'] == gene]
        
        # get total number of phosphorylatable sites
        site_counts['S'] += seq.count('S')
        site_counts['T'] += seq.count('T')
        site_counts['Y'] += seq.count('Y')

        possible_phosphorylation_sites = list()
        for i, c in enumerate(seq):
            if c in ['S', 'T', 'Y']:
                possible_phosphorylation_sites.append(''.join([c, str(i)]))

        # use UniProt accession to get predicted residues in UniProt Accession
        pred_result = pred_df[pred_df['ID'].str.split(pat='|').apply(lambda x: x[1]) == gene_map['UniProt Accession'].values[0]]
        pred_residues = list()
        if len(pred_result) != 0:
            pred_residues = pred_result['Phosphorylated residue'].values[0].split(';')

        if tool == 'musite':
            temp = list()
            for x in pred_residues:
                if x == 'SEQ>1000':
                    break
                r = x.split(':')[0]
                score = float(x.split(':')[1])
                if score > 0.5:
                    temp.append(r)
            pred_residues = temp

        if tool == 'sitetack':
            temp = list()
            for x in pred_residues:
                if x == 'NA' or len(x) == 0:
                    break
                y = x.split(':')[0]
                r = y[0]
                loc = int(y[1:]) - 1
                score = float(x.split(':')[1])
                if score > 0.5:
                    temp.append(''.join([r, str(loc)]))
            pred_residues = temp

        if len(pred_residues) == 1:
            if pred_residues[0] == 'NA':
                pred_residues = list()

        if tool == 'ptmgpt2' or tool == 'netphos':
            temp = list()
            for x in pred_residues:
                r = x[0]
                loc = int(x[1:]) - 1
                temp.append(''.join([r, str(loc)]))
            pred_residues = temp

        # get TP, TN, FP, FN
        if gene_map['GI'].values[0] == '': # there are no real phosphorylation sites 
            # any predictions are false positives
            for x in pred_residues:
                r = x[0]
                false_positive[r] += 1
                if x in possible_phosphorylation_sites:
                    possible_phosphorylation_sites.remove(x)
                else:
                    print(x)

            # the S,T, and Y sites where there were no predictions are true negatives
            for x in possible_phosphorylation_sites:
                r = x[0]
                true_negative[r] += 1
            
        else:
            # use Genbank GI from map to get positive residues in s_pneumoniae_pos
            pos_result = pos_df[pos_df['NCBI_Acc#'] == gene_map['GI'].values[0]]
            pos_residues = pos_result['Phosphorylated residue'].values[0].split(';')
            pos_residues = [''.join([x[0], str(int(x[1:]) - 1)]) for x in pos_residues] # musite indexed at 1 and not 0

            # see what predictions are real (true postive) or not predicted correctly (false negative)
            modified_pred_residues = pred_residues
            for x in pos_residues:
                r = x[0]
                #loc = x[1:]

                if x in pred_residues: # real phos site that was predicted (true positive)
                    true_positive[r] += 1
                    modified_pred_residues.remove(x)
                else: # was a real phos site but was not predicted (false negative)
                    false_negative[r] += 1

                if x in possible_phosphorylation_sites:
                    possible_phosphorylation_sites.remove(x)
                    
            pred_residues = modified_pred_residues

            # these are the residues that were predicted to be phos sites but really are not (false positive)
            for x in pred_residues: 
                r = x[0]
                #loc = x[1:]
                false_positive[r] += 1

                if x in possible_phosphorylation_sites:
                    possible_phosphorylation_sites.remove(x)

            # these are the sites that did not have a real phos ptm and were not predicted to have it (true negative)
            for x in possible_phosphorylation_sites: 
                r = x[0]
                #loc = x[1:]
                true_negative[r] += 1
    
    print('TP', true_positive)
    print('TN', true_negative)
    print('FP', false_positive)
    print('FN', false_negative)

def safe_concat_arrays(a, b):
    # Convert NaNs or non-arrays to empty arrays
    a = a if isinstance(a, np.ndarray) else np.array([])
    b = b if isinstance(b, np.ndarray) else np.array([])
    return np.concatenate((a, b))

def format_residue_scores(r, locations, scores):
    if not isinstance(locations, np.ndarray) or not isinstance(scores, np.ndarray):
        return ''
    return ';'.join([f'{r}{int(loc)}:{round(score, 8)}' for loc, score in zip(locations, scores)])

def process_sitetack(st_path, y_path):
    sitetack_st = pd.read_csv(st_path, index_col=False)
    sitetack_y = pd.read_csv(y_path, index_col=False)

    sitetack_st['No labels model'] = sitetack_st['No labels model'].apply(lambda x: np.fromstring(x.strip('[]'), sep=' '))
    sitetack_st['With PTM labels model'] = sitetack_st['With PTM labels model'].apply(lambda x: np.fromstring(x.strip('[]'), sep=' '))
    sitetack_st['S'] = sitetack_st['S'].apply(lambda x: np.fromstring(x.strip('[]'), sep=',', dtype=np.int64))
    sitetack_st['T'] = sitetack_st['T'].apply(lambda x: np.fromstring(x.strip('[]'), sep=',', dtype=np.int64))
    
    sitetack_y['No labels model'] = sitetack_y['No labels model'].apply(lambda x: np.fromstring(x.strip('[]'), sep=' '))
    sitetack_y['With PTM labels model'] = sitetack_y['With PTM labels model'].apply(lambda x: np.fromstring(x.strip('[]'), sep=' '))
    sitetack_y['Y'] = sitetack_y['Y'].apply(lambda x: np.fromstring(x.strip('[]'), sep=',', dtype=np.int64))

    sitetack_st = sitetack_st.rename(columns={
        'No labels model': 'No labels model_df1',
        'With PTM labels model': 'With PTM labels model_df1'
    })
    sitetack_y = sitetack_y.rename(columns={
        'No labels model': 'No labels model_df2',
        'With PTM labels model': 'With PTM labels model_df2'
    })

    sitetack = pd.merge(sitetack_st, sitetack_y, on='ID', how='outer')

    sitetack['No labels model'] = sitetack.apply(
        lambda row: safe_concat_arrays(row['No labels model_df1'], row['No labels model_df2']), axis=1
    )
    sitetack['With PTM labels model'] = sitetack.apply(
        lambda row: safe_concat_arrays(row['With PTM labels model_df1'], row['With PTM labels model_df2']), axis=1
    )

    S = sitetack.apply(
        lambda row: format_residue_scores('S', row['S'], row['With PTM labels model']), axis=1
    )
    T = sitetack.apply(
        lambda row: format_residue_scores('T', row['T'], row['With PTM labels model']), axis=1
    )
    Y = sitetack.apply(
        lambda row: format_residue_scores('Y', row['Y'], row['With PTM labels model']), axis=1
    )

    sitetack['Phosphorylated residue'] = S + ';' + T + ';' + Y

    sitetack = sitetack.drop(columns=['No labels model_df1', 'No labels model_df2',
                                    'With PTM labels model_df1', 'With PTM labels model_df2'])
    sitetack = sitetack[['ID', 'Phosphorylated residue']].fillna('NA')

    return sitetack

def process_netphos(path):
    netphos = pd.read_csv(path, index_col=False).fillna('NA')
    netphos['ProteinID'] = netphos['ProteinID'].apply(lambda id: id.replace('_', '|'))
    netphos['ResiduePositions'] = netphos['ResiduePositions'].apply(lambda res: ';'.join([element.strip() for element in res.split(',')]))
    
    netphos = netphos.rename(columns={
        'ProteinID': 'ID',
        'ResiduePositions': 'Phosphorylated residue'
    })

    return netphos

def main():
    netphos_mutans = process_netphos('./predictions/Net_mutans_merge.csv')
    netphos_pneumoniae = process_netphos('./predictions/Net_Pneumo_merged__5.csv')

    sitetack_mutans = process_sitetack('./predictions/Merged_Sitetack_Mutans_Run_1_ST.csv', './predictions/Mutans_Y_Sitetack_pred_merged.csv')
    sitetack_pneumoniae = process_sitetack('./predictions/Merged_pneumo_Sitetac_ST.csv', './predictions/Pneumo_Y_sitetack_pred_merged.csv')

    # IDs in these prediction file are in uniprot format
    musite_mutans = pd.read_csv('./predictions/musite_sm.csv', index_col=False).fillna('NA')
    musite_pneumoniae = pd.read_csv('./predictions/musite_sp.csv', index_col=False).fillna('NA')
    ptmgpt2_mutans = pd.read_csv('./predictions/ptmgpt2_sm.csv', index_col=False).fillna('NA')
    ptmgpt2_pneumoniae = pd.read_csv('./predictions/ptmgpt2_sp.csv', index_col=False).fillna('NA')

    musite_mutans_breakdown = musite_site_breakdown(musite_mutans, 0.5, 'proteome_sm.fasta')
    musite_pneumoniae_breakdown = musite_site_breakdown(musite_pneumoniae, 0.5, 'proteome_sp.fasta')
    ptmgpt2_mutans_breakdown = ptmgpt2_site_breakdown(ptmgpt2_mutans, 'proteome_sm.fasta')
    ptmgpt2_pneumoniae_breakdown = ptmgpt2_site_breakdown(ptmgpt2_pneumoniae, 'proteome_sp.fasta')

    sitetack_mutans_breakdown = sitetack_site_breakdown(sitetack_mutans, 0.5, 'proteome_sm.fasta')
    netphos_mutans_breakdown = ptmgpt2_site_breakdown(netphos_mutans, 'proteome_sm.fasta')

    print('MusiteDeep', musite_mutans_breakdown[0])
    print('PTM GPT2', ptmgpt2_mutans_breakdown[0])
    print('Sitetack', sitetack_mutans_breakdown[0])
    print('NetPhos 3.1', netphos_mutans_breakdown[0])

    s_pneumoniae_proteome = read_proteome_fasta('proteome_sp.fasta')
    temp = list()
    x = 0

    for i in list(musite_pneumoniae_breakdown[1].keys()):
        if musite_pneumoniae_breakdown[1][i][0] == -1:
            x += s_pneumoniae_proteome[i].count('S')
            x += s_pneumoniae_proteome[i].count('T')
            x += s_pneumoniae_proteome[i].count('Y')

    print(x)

    # ID in this file uses deprecated gi identifiers
    s_pneumoniae_pos = pd.read_csv('sp_pos_con.csv', index_col=False)

    s_pneumoniae_pos_seq = {}
    with open('sp_pos_con_seq.fasta', 'r') as file:
        id = ''
        for i, line in enumerate(file.readlines()):
            line = line.strip()
            if line[0] == '>':
                id = line.split(',')[0][1:]
                accession = line.split(',')[1]
                s_pneumoniae_pos_seq[id] = [accession, '']
            else:
                s_pneumoniae_pos_seq[id][1] += line

    #s_pneumoniae_accessions = read_proteome_fasta('proteome_sp.fasta').keys()
    #_pneumoniae_whole_proteome_map_data = map_uniprot_accession(s_pneumoniae_accessions)
    #s_pneumoniae_phos_pos_map_data = map_gi(s_pneumoniae_pos_seq)

    #file = open('sp_wp_map.csv', 'w')
    #file.write('GI,GenBank Accession,UniProt Accession,UniProt Entry ID\n')
    #for id in s_pneumoniae_whole_proteome_map_data:
    #    entry = s_pneumoniae_whole_proteome_map_data[id]
    #    file.write(f'{entry[0]},{entry[1]},{entry[2]},{entry[3]}\n')
    #file.close()

    #file = open('sp_phos_pos_map.csv', 'w')
    #file.write('GI,Protein ID,UniProt Accession,UniProt Entry ID\n')
    #for id in s_pneumoniae_phos_pos_map_data:
    #    entry = s_pneumoniae_phos_pos_map_data[id]
    #    file.write(f'{entry[0]},{entry[1]},{entry[2]},{entry[3]}\n')
    #file.close()
    
    s_pneumoniae_whole_proteome_map = pd.read_csv('sp_wp_map.csv', index_col=False).fillna('')
    s_pneumoniae_whole_proteome_map.set_index('UniProt Accession', inplace=True)

    s_pneumoniae_pos_map = pd.read_csv('sp_phos_pos_map.csv', index_col=False).fillna('')
    s_pneumoniae_pos_map.set_index('UniProt Accession', inplace=True)

    # get the gi information from s_pneumoniae_pos_map for the real phosphorylated proteins into the main map df
    s_pneumoniae_whole_proteome_map.update(s_pneumoniae_pos_map)
    s_pneumoniae_whole_proteome_map.reset_index(inplace=True)

    s_pneumoniae_proteome = read_proteome_fasta('proteome_sp.fasta')

    print('\nSP NetPhos')
    calculate_statistics(s_pneumoniae_whole_proteome_map, 'netphos', s_pneumoniae_proteome, netphos_pneumoniae, s_pneumoniae_pos)

    print('\nSP Sitetack')
    calculate_statistics(s_pneumoniae_whole_proteome_map, 'sitetack', s_pneumoniae_proteome, sitetack_pneumoniae, s_pneumoniae_pos)

    print('\nSP MusiteDeep')
    calculate_statistics(s_pneumoniae_whole_proteome_map, 'musite', s_pneumoniae_proteome, musite_pneumoniae, s_pneumoniae_pos)
    
    print('\nSP PTMGPT2')
    calculate_statistics(s_pneumoniae_whole_proteome_map, 'ptmgpt2', s_pneumoniae_proteome, ptmgpt2_pneumoniae, s_pneumoniae_pos)



if __name__ == '__main__':
    main()