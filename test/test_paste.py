# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 13:57:12 2024

@author: chuis
"""
## This is a script adjusted from bioblend to connect with the galaxy server and upload sequence data.
## Based on: https://github.com/naturalis/bioblend-test/blob/main/client.py

##Dependencies
from bioblend import galaxy
#from bioblend.galaxy.tools.inputs import inputs, dataset
import os
#import requests
#from dotenv import load_dotenv, find_dotenv

# Endpoints for the Naturalis production instance. Domain can also be an IP address.
domain = 'galaxy.naturalis.nl'
dlbase = f'https://{domain}/api/datasets'

# The API key can be obtained from https://galaxy.naturalis.nl/user/api_key
#load_dotenv(find_dotenv(usecwd=True))
api_key = os.environ.get('GALAXY_API_KEY')
gi = galaxy.GalaxyInstance(domain, key=api_key)

# Get or create a new history. The civilized thing would be to delete this when done.
history_name = 'barcode-phylogenetic-diversity'
histories = gi.histories.get_histories(name=history_name)
history = gi.histories.create_history(history_name) if not histories else histories[0]

## Load sequence data into galaxy
with open('../data/sequences.fasta', 'r') as file:
    seqs_contents = file.read()
seqs_id = gi.tools.paste_content(seqs_contents, history['id'], file_type='fasta')['outputs'][0]['id']

##Write id to textfile next
with open('../results/galaxy-ids/paste_id.txt', 'w') as idfile:
    idfile.write(f'{seqs_id}')

##gi.histories.delete_history(history['id'])