##Script to extract data from galaxy server
#dependencies
from bioblend import galaxy
import os

#connect to galaxy
domain = 'galaxy.naturalis.nl'
api_key = os.environ.get('GALAXY_API_KEY')
gi = galaxy.GalaxyInstance(domain, key=api_key)

#get id and download data
id_path = (snakemake.input[0])
with open(f'{id_path}', 'r') as id_file:
    extract_id = id_file.read()

output_path = (snakemake.output[0])
gi.datasets.download_dataset(dataset_id=extract_id, file_path=output_path, use_default_filename=False)