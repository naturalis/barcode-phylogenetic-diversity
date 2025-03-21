## Script that runs the right tools with the right parameters

# Dependencies
from bioblend import galaxy
import os
import dotenv
import logging
import parameters
import requests

# Configure logging
logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG,
                    format='%(asctime)s:%(levelname)s:%(message)s')

# Endpoints for the Naturalis production instance. Domain can also be an IP address.
# The API key can be obtained from https://galaxy.naturalis.nl/user/api_key
logging.info('Accessing galaxy instance with API key.')
domain = 'https://galaxy.naturalis.nl'
response = requests.get(domain)
logging.debug(f'Response from {domain}: {response}')

dotenv.load_dotenv(dotenv.find_dotenv(usecwd=True))
api_key = os.environ.get('GALAXY_API_KEY')
logging.debug(f"API Key provided: {bool(api_key)}")
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
data_ids = []
for data_path in snakemake.input[1:]:
    with open(f'{data_path}', 'r') as data_f:
        data_id = data_f.read().split('\n')
    data_ids.extend(data_id)

# Get the right tool according to tool_name
# The out_file variable is set for dealing with different export parameters (i.e. alpha/beta)
logging.info('Getting the tool and setting out_file variable.')

tool = gi.tools.get_tools(name=tool_name)
if not tool:
    raise ValueError(f"Tool '{tool_name}' not found in Galaxy instance. Check definition.")

tool_id = tool[0]['id']
out_file = os.path.basename(snakemake.output[0])
params = parameters.set(tool_name, data_ids, out_file)

# Parameter setting and running the tool
logging.info(f"Running the appropriate tool: {tool_name}")
try:
    results = gi.tools.run_tool(history['id'], tool_id, params)
    job = results['jobs'][0]
    jc = galaxy.jobs.JobsClient(galaxy_instance=gi)
    jc.wait_for_job(job_id=job['id'])

except ConnectionError as e:
    logging.error(f"Failed to run {tool_name}: {e}")

# Saving output ids
logging.info(f"Saving output id(s) of {tool_name} to temporary file(s) for further use.")
for res in range(len(results['outputs'])):
    output_path = snakemake.output[res]
    result_id = results['outputs'][res]['id']
    with open(f'{output_path}', 'w') as output_f:
        output_f.write(result_id)
