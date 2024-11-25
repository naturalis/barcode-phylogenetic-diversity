
from bioblend import galaxy
from bioblend.galaxy.tools.inputs import inputs, dataset
import os
import requests
#from dotenv import load_dotenv, find_dotenv

# Load environment variables
#load_dotenv(find_dotenv(usecwd=True))
api_key = os.environ.get('GALAXY_API_KEY')

# Galaxy instance
domain = 'galaxy.naturalis.nl'
try:
    gi = galaxy.GalaxyInstance(domain, key=api_key)
    print("Access obtained from Galaxy")
except requests.RequestException as req_err:
    print(f"Request error Galaxy access: {req_err}")
    raise
except Exception as e:
    print(f"Unexpected error Galaxy access: {e}")
    raise

# Create or get history
history_name = 'qiime-dentity'
histories = gi.histories.get_histories(name=history_name)
history = gi.histories.create_history(history_name) if not histories else histories[0]

# Load data into Galaxy
fasta_path = 'data/sequences.fasta'
with open(fasta_path, 'r') as file:
    seqs_contents = file.read()

try:
    seqs_id = gi.tools.paste_content(seqs_contents, history['id'], file_type='fasta')['outputs'][0]['id']
    print(f"Data pasted into history: {seqs_id}")
except Exception as e:
    print(f"Error pasting data: {e}")
    raise


# Get import tool and import data into qiime2
try:
    q2_import_tool = gi.tools.get_tools(name='qiime2 tools import')[0]
    params_imp = inputs().set('import_root|type', 'FeatureData__ob__Sequence__cb__')\
                         .set('import_root|__q2galaxy__GUI__cond__format__|format', 'DNAFASTAFormat')\
                         .set('import_root|__q2galaxy__GUI__cond__format__|data', dataset(seqs_id))

    imp_results = gi.tools.run_tool(history['id'], q2_import_tool['id'], params_imp)
    print(f"Import tool run initiated: {imp_results}")

    job = imp_results['jobs'][0]
    jc = galaxy.jobs.JobsClient(galaxy_instance=gi)
    jc.wait_for_job(job_id=job['id'])

    imp_id = imp_results['outputs'][0]['id']

    # Print detailed dataset permissions info
    dataset_info = gi.datasets.show_dataset(imp_id)
    print(f"Dataset detailed info: {dataset_info}")

    url = f'https://{domain}/api/datasets/{imp_id}/display?to_ext=data'
    imp_response = requests.get(url)
    imp_response.raise_for_status()

    with open("ignore/imported_data.qza", "w") as file:
        file.write(imp_response.text)

    print("Import and retrieval successful")

except requests.RequestException as req_err:
    print(f"Request error importing data: {req_err}")
    raise
except Exception as e:
    print(f"Unexpected error importing data: {e}")
    raise

print("Response to request imported data is:", imp_response)

#gi.histories.delete_history(histories[0]['id'])



