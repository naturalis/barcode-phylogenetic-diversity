##Script to extract data from galaxy again
#dependencies
from bioblend import galaxy
import os

#define id and output paths, later defined as babysnake input
id_path = 'results/galaxy-ids/paste_id.txt'
output_path = 'results/div-metrics/'

#connect to galaxy, set up history
domain = 'galaxy.naturalis.nl'
api_key = os.environ.get('GALAXY_API_KEY')
gi = galaxy.GalaxyInstance(domain, key=api_key)

history_name = 'barcode-phylogenetic-diversity'
histories = gi.histories.get_histories(name=history_name)
history = gi.histories.create_history(history_name) if not histories else histories[0]

#get id and download data
with open(f'{id_path}', 'r') as id_file:
    extract_id = id_file.read()

extract_data = gi.datasets.download_dataset(dataset_id=extract_id, file_path=output_path)


