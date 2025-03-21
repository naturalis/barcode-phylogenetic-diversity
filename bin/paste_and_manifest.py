# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 13:57:12 2024

@author: chuis
"""
## This is a script to connect with the galaxy server and upload sequence data.

##Dependencies
from bioblend import galaxy
import os
import logging
import csv

# Configure logging
logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG,
                    format='%(asctime)s:%(levelname)s:%(message)s')

# Endpoints for the Naturalis production instance. Domain can also be an IP address.
domain = 'https://galaxy.naturalis.nl'
logging.debug(f"Accessing galaxy instance: {domain}")

# The API key can be obtained from https://galaxy.naturalis.nl/user/api_key
api_key = os.environ.get('GALAXY_API_KEY')
logging.debug(f"API Key provided: {bool(api_key)}")
gi = galaxy.GalaxyInstance(domain, key=api_key)

# Get or create a new history. The civilized thing would be to delete this when done.
history_name = 'barcode-phylogenetic-diversity'
histories = gi.histories.get_histories(name=history_name)
history = gi.histories.create_history(history_name) if not histories else histories[0]
logging.info(f"Using history: {history_name}")

# Loop for pasting each sequences file
# Samples is an empty list where the manifest data will be stored
samples = []
data_dir = snakemake.input[0]

for file in os.listdir(data_dir):
    logging.info(f"Processing dataset: {file}")
# Dataset is pasted and id retrieved
    if file.endswith('.fastq.gz'):
        seq_path = f'{data_dir}/{file}'

        try:
            seq_id = gi.tools.upload_file(seq_path, history['id'])['outputs'][0]['id']
            logging.info(f"Dataset {file} pasted successfully in Galaxy with ID: {seq_id}")

        except Exception as e:
            logging.error(f"Failed to paste dataset {file}: {e}")
            continue

# Now sample_id, absolute filepath and direction are obtained for manifest file
        logging.info(f"Extracting sample_id, absolute_filepath, and direction from {file} ")
        seq_uuid = gi.datasets.show_dataset(seq_id)['uuid']
        absolute_filepath = f"/data/files/{seq_uuid[0]}/{seq_uuid[1]}/{seq_uuid[2]}/dataset_{seq_uuid}.dat"
        base = file.replace('.fastq.gz', '')

        if 'R1' in base:
            direction = 'forward'
            remove = 'R1'

        elif 'R2' in base:
            direction = 'reverse'
            remove = 'R2'

        sample_id = base.replace(f'_{remove}', '')
        manifest_data = {
            'sample-id': f'{sample_id}',
            'absolute-filepath': f'{absolute_filepath}',
            'direction': f'{direction}'
        }
        samples.append(manifest_data)
        continue

    else:
        continue

# Saving the manifest data as a csv file
man_csv_out = snakemake.output[0]
field_names = ['sample-id', 'absolute-filepath', 'direction']

with open(f'{man_csv_out}', 'w') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=field_names)
    writer.writeheader()
    writer.writerows(samples)

# Pasting manifest in Galaxy
try:
    man_id = gi.tools.upload_file(man_csv_out, history['id'])['outputs'][0]['id']
    man_uuid = gi.datasets.show_dataset(man_id)['uuid']
    man_uuid_path = f"/data/files/{man_uuid[0]}/{man_uuid[1]}/{man_uuid[2]}/dataset_{man_uuid}.dat"
    logging.info(f"Manifest file pasted successfully in Galaxy with UUID: {man_uuid}")

except Exception as e:
    logging.error(f"Failed to paste manifest file: {e}")

# Saving id of manifest file
man_uuid_out = snakemake.output[1]

with open(f'{man_uuid_out}', 'w') as manout_f:
    manout_f.write(man_uuid_path)


