### Testing things manually using Moving Pictures data, mainly params
## Dependencies
from bioblend import galaxy
from bioblend.galaxy.tools.inputs import inputs, dataset
import os
import requests
from dotenv import load_dotenv, find_dotenv

# Endpoints for the Naturalis production instance. Domain can also be an IP address.
domain = 'galaxy.naturalis.nl'
dlbase = f'https://{domain}/api/datasets'

# The API key can be obtained from https://galaxy.naturalis.nl/user/api_key
load_dotenv(find_dotenv(usecwd=True))
api_key = os.environ.get('GALAXY_API_KEY')
gi = galaxy.GalaxyInstance(domain, key=api_key)

# Get or create a new history. The civilized thing would be to delete this when done.
histories = gi.histories.get_histories(name='qiime-dentity')
if len(histories) == 0:
    history = gi.histories.create_history('qiime-dentity')
else:
    history = histories[0]

## Load data into galaxy
with open('data/dna-sequences.fasta', 'r') as file:
    seqs_contents = file.read()
seqs_id = gi.tools.paste_content(seqs_contents, history['id'], file_type='fasta')['outputs'][0]['id']

## Get import tool & import data into qiime2
q2_import = gi.tools.get_tools(name='qiime2 tools import')[0]

params_imp = inputs().set('import_root|type','FeatureData[Sequence]')\
    .set('__q2galaxy__GUI__cond__format__|format', 'DNAFASTAFormat')\
    .set('__q2galaxy__GUI__cond__format__|data', dataset(seqs_id))

## Run import tool and get response
try:
    imp_results = gi.tools.run_tool(history['id'], q2_import['id'], params_imp)

    job = imp_results['jobs'][0]
    jc = galaxy.jobs.JobsClient(galaxy_instance=gi)
    jc.wait_for_job(job_id=job['id'])

    imp_id = imp_results['outputs'][0]['id']
    url = f'{dlbase}/{imp_id}/display?to_ext=dat'
    imp_response = requests.get(url)

    with open("imported_data.qza", "w") as file:
        file.write(imp_response.text)

except ConnectionError as e:
    print(e.body)

print("Response to request imported data is:", imp_response)

## Just try to export it again for now
q2_export = gi.tools.get_tools(name='qiime2 tools export')[0]

params_exp = inputs().set('input', dataset(imp_id))\
    .set('input|type_peek', 'FeatureData__ob__Sequence__cb__')\
    .set('input|fmt_peek', 'FeatureData__ob__Sequence__cb__')\
    .set('output_format', 'DNASequencesDirectoryFormat')

## Run export tool and get response
try:
    exp_results = gi.tools.run_tool(history['id'], q2_export['id'], params_exp)

    job2 = exp_results['jobs'][0]
    jc = galaxy.jobs.JobsClient(galaxy_instance=gi)
    jc.wait_for_job(job_id=job2['id'])

    exp_id = exp_results['outputs'][0]['id']
    url2 = f'{dlbase}/{exp_id}/display?to_ext=dat'
    exp_response = requests.get(url2)

    with open("exported_data.fa", "w") as file:
        file.write(exp_response.text)

except ConnectionError as e:
    print(e.body)

print("Response to request exported data is:", exp_response)


