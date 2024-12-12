## Running the right tools with the right parameters
#Dependencies
from bioblend import galaxy
from bioblend.galaxy.tools.inputs import inputs, dataset
import os
import ast
#import dotenv

# Endpoints for the Naturalis production instance. Domain can also be an IP address.
domain = 'galaxy.naturalis.nl'
dlbase = f'https://{domain}/api/datasets'

# The API key can be obtained from https://galaxy.naturalis.nl/user/api_key
#dotenv.load_dotenv(dotenv.find_dotenv(usecwd=True))
api_key = os.environ.get('GALAXY_API_KEY')
gi = galaxy.GalaxyInstance(domain, key=api_key)

# Get or create a new history. The civilized thing would be to delete this when done.
history_name = 'barcode-phylogenetic-diversity'
histories = gi.histories.get_histories(name=history_name)
history = gi.histories.create_history(history_name) if not histories else histories[0]


# Saving the snakemake input in variables to use for parameter setting
tool_path = snakemake.input[0]
with open(f'{tool_path}', 'r') as name_f:
    tool_name = name_f.read()

data_path = snakemake.input[1]
with open(f'{data_path}', 'r') as data_f:
    data_id = data_f.read()


## Get the right parameters, data and tool name goes in here
if tool_name == 'qiime2 tools import':
    params = inputs().set('import_root|type', 'FeatureData__ob__Sequence__cb__') \
        .set('__q2galaxy__GUI__cond__format__|format', 'DNAFASTAFormat') \
        .set('__q2galaxy__GUI__cond__format__|data', dataset(data_id))

elif tool_name == 'qiime2 tools export':
    params = inputs().set('input', dataset(data_id)) \
        .set('input|type_peek', 'FeatureData__ob__Sequence__cb__') \
        .set('input|fmt_peek', 'FeatureData__ob__Sequence__cb__') \
        .set('output_format', 'DNAFASTAFormat')
else:
    print('No parameters for this tool. Check your tool list or definition of the tool.')


## Get the id of the tool from previously made dictionary (converted back from text file)
dict_path = snakemake.input[2]

with open(f'{dict_path}', 'r') as tools_f:
    tools_str = tools_f.read()

tools_dict = ast.literal_eval(tools_str)
tool_id = tools_dict[f'{tool_name}']


## Running the tool in question using dictionary with tool name and ids
try:
    results = gi.tools.run_tool(history['id'], tool_id, params)

    job = results['jobs'][0]
    jc = galaxy.jobs.JobsClient(galaxy_instance=gi)
    jc.wait_for_job(job_id=job['id'])

except ConnectionError as e:
    print(e.body)


##Writing output into id files
result_id = results['outputs'][0]['id']
output_path = snakemake.output[0]
with open(f'{output_path}', 'w') as output_f:
    output_f.write(result_id)

#Clause for dada2 having two results (concept, confusing which is seq and which is tab)?
#if len(results) == 2:
#    result_id2 = results['outputs'][1]['id']
#    with open(f'../galaxy-ids/{tool_name}_out2.txt', 'w') as file:
#        file.write(result_id2.text)
