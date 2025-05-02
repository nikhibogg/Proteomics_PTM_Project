import requests
import json
import os
import concurrent.futures
from Bio import SeqIO

def predict_phosphorylation(fasta, i, num_sequences, model):
    id = fasta.id
    seq = fasta.seq

    url = f'https://api.musite.net/musitedeep/{model}/{seq}'
    myResponse = requests.get(url)

    print(f'Predicting Phosphorylation sites for {id} ({i + 1}/{num_sequences})')

    ptm_sites = []

    if(myResponse.ok):
        # In this Example, jData are prediction results from MusiteDeep predictor
        jData = json.loads(myResponse.content.decode('utf-8'))
        if "Error" in jData.keys(): 
            return (id, "SEQ>1000")
        else:
            results = jData['Results']
            for r in results:
                ps = int(r['Position']) - 1
                score = float(r['PTMscores'].split(':')[1])
                residue = r['Residue']

                ptm_sites.append(f'{residue}{ps}:{score}')

            return (id, ';'.join(ptm_sites))
    return (id, 'RESPONSE_ERROR')

def main():
    #users can only select from the following PTM models:
    modeloptions = ["Phosphoserine_Phosphothreonine",
                    "Phosphotyrosine",
                    "N-linked_glycosylation",
                    "O-linked_glycosylation",
                    "Ubiquitination",
                    "SUMOylation",
                    "N6-acetyllysine",
                    "Methylarginine",
                    "Methyllysine",
                    "Pyrrolidone_carboxylic_acid",
                    "S-palmitoyl_cysteine",
                    "Hydroxyproline",
                    "Hydroxylysine"]
    
    model = f'{modeloptions[0]};{modeloptions[1]}'

    print('Predicting SP Phosphorylation sites...')

    fasta_file = 'proteome_sp.fasta'
    fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')

    fh = open(fasta_file)
    num_sequences = 0
    for line in fh:
        if line.startswith(">"):
            num_sequences += 1
    fh.close()

    file = open('./predictions/musite_sp.csv', 'w')
    file.write('ID,Phosphorylated residue\n')
    file.flush()
    os.fsync(file.fileno())

    with concurrent.futures.ProcessPoolExecutor(max_workers=10) as executor:
        futures = []
        for i, fasta in enumerate(fasta_sequences):
            futures.append(executor.submit(predict_phosphorylation, fasta, i, num_sequences, model))

        # Collect and write results as they complete
        for future in concurrent.futures.as_completed(futures):
            id, pred_ptm_sites = future.result()
            file.write(f'{id},{pred_ptm_sites}\n')
            file.flush()
            os.fsync(file.fileno())

    file.close()

if __name__ == '__main__':
    main()