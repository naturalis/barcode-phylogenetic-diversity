## Running the right tools with the right parameters

#Major dependencies
from bioblend import galaxy
import os
#import dotenv
import logging
import parameters

# configure logging
logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG,
                    format='%(asctime)s:%(levelname)s:%(message)s')

# Endpoints for the Naturalis production instance. Domain can also be an IP address.
# The API key can be obtained from https://galaxy.naturalis.nl/user/api_key
logging.debug('Accessing galaxy instance with API key.')
domain = 'galaxy.naturalis.nl'
dlbase = f'https://{domain}/api/datasets'

#dotenv.load_dotenv(dotenv.find_dotenv(usecwd=True))
api_key = os.environ.get('GALAXY_API_KEY')
gi = galaxy.GalaxyInstance(domain, key=api_key)

# Get or create a new history. The civilized thing would be to delete this when done.
logging.debug('Getting the right history.')
history_name = 'barcode-phylogenetic-diversity'
histories = gi.histories.get_histories(name=history_name)
history = histories[0]


# Storing the snakemake tool name to use for parameter setting and getting the tool
logging.debug('Preparing snakemake input for parameter setting.')
tool_path = snakemake.input[0]
with open(f'{tool_path}', 'r') as name_f:
    tool_name = name_f.read()

# Making a list of input ids to use for parameter setting
data_ids =[]
index = 0
while snakemake.input[index] != snakemake.input[-1]:
    data_path = snakemake.input[index+1]
    with open(f'{data_path}', 'r') as data_f:
        data_ids.append(data_f.read())
    index += 1

## Get the right parameters and tool according to tool_name and data_ids
# out_file variable is set for dealing with different export parameters (i.e. alpha/beta)
logging.debug('Setting parameters and getting tool.')

out_file = os.path.basename(snakemake.output[0])
params = parameters.set(tool_name, data_ids, out_file)

tool = gi.tools.get_tools(name=tool_name)
tool_id = tool[0]['id']

## Running the tool in question
logging.debug('Running the appropriate tool...')
try:
    results = gi.tools.run_tool(history['id'], tool_id, params)

    job = results['jobs'][0]
    jc = galaxy.jobs.JobsClient(galaxy_instance=gi)
    jc.wait_for_job(job_id=job['id'])

except ConnectionError as e:
    print(e.body)


##Writing output into id files
logging.debug('Saving output id(s) to temporary file(s) for further use.')
for res in range(len(results['outputs'])):
    output_path = snakemake.output[res]
    result_id = results['outputs'][res]['id']
    with open(f'{output_path}', 'w') as output_f:
        output_f.write(result_id)

