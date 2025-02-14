##Script to extract data from galaxy server
#dependencies
from bioblend import galaxy
import os
import logging

# configure logging
logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG,
                    format='%(asctime)s:%(levelname)s:%(message)s')

#connect to galaxy
logging.debug('Accessing galaxy instance with API key.')
domain = 'https://galaxy.naturalis.nl'
api_key = os.environ.get('GALAXY_API_KEY')
gi = galaxy.GalaxyInstance(domain, key=api_key)

#get id(s) and download data
logging.debug('Getting the id of the exported qiime2 object.')
extract_ids = []
id_path = snakemake.input[0]
with open(f'{id_path}', 'r') as id_file:
    extract_ids.append(id_file.read())

if not snakemake.input[0] == snakemake.input[-1]:
    id_path2 = snakemake.input[1]
    with open(f'{id_path2}', 'r') as id_file2:
        extract_ids.append(id_file2.read())

logging.debug('Downloading exported qiime2 object...')
for metric in range(len(extract_ids)):
    output_path = snakemake.output[metric]
    gi.datasets.download_dataset(dataset_id=extract_ids[metric], file_path=output_path, use_default_filename=False)