######## snakemake preamble start (automatically inserted, do not edit) ########
import sys;sys.path.extend(['/home/chunix/mambaforge2/envs/barcode-phylogenetic-diversity/lib/python3.12/site-packages', '/mnt/c/users/chuis/documents/biology_biosus/masterstage_2/barcode-phylogenetic-diversity/workflow', '/home/chunix/mambaforge2/envs/barcode-phylogenetic-diversity/bin', '/home/chunix/mambaforge2/envs/barcode-phylogenetic-diversity/lib/python3.12', '/home/chunix/mambaforge2/envs/barcode-phylogenetic-diversity/lib/python3.12/lib-dynload', '/home/chunix/mambaforge2/envs/barcode-phylogenetic-diversity/lib/python3.12/site-packages', '/home/chunix/.cache/snakemake/snakemake/source-cache/runtime-cache/tmp009dhr5_/file/mnt/c/users/chuis/documents/biology_biosus/masterstage_2/barcode-phylogenetic-diversity/workflow/../bin', '/mnt/c/users/chuis/documents/biology_biosus/masterstage_2/barcode-phylogenetic-diversity/workflow/../bin']);import pickle;from snakemake import script;script.snakemake = pickle.loads(b'\x80\x04\x95*\x04\x00\x00\x00\x00\x00\x00\x8c\x10snakemake.script\x94\x8c\tSnakemake\x94\x93\x94)\x81\x94}\x94(\x8c\x05input\x94\x8c\x0csnakemake.io\x94\x8c\nInputFiles\x94\x93\x94)\x81\x94(\x8c$../results/galaxy-ids/dada2_tool.txt\x94\x8c$../results/galaxy-ids/import_out.txt\x94e}\x94(\x8c\x06_names\x94}\x94\x8c\x12_allowed_overrides\x94]\x94(\x8c\x05index\x94\x8c\x04sort\x94eh\x11h\x06\x8c\x0eAttributeGuard\x94\x93\x94)\x81\x94}\x94\x8c\x04name\x94h\x11sbh\x12h\x14)\x81\x94}\x94h\x17h\x12sbub\x8c\x06output\x94h\x06\x8c\x0bOutputFiles\x94\x93\x94)\x81\x94(\x8c#../results/galaxy-ids/dada2_out.txt\x94\x8c$../results/galaxy-ids/dada2_out2.txt\x94\x8c$../results/galaxy-ids/dada2_out3.txt\x94e}\x94(h\r}\x94h\x0f]\x94(h\x11h\x12eh\x11h\x14)\x81\x94}\x94h\x17h\x11sbh\x12h\x14)\x81\x94}\x94h\x17h\x12sbub\x8c\r_params_store\x94h\x06\x8c\x06Params\x94\x93\x94)\x81\x94}\x94(h\r}\x94h\x0f]\x94(h\x11h\x12eh\x11h\x14)\x81\x94}\x94h\x17h\x11sbh\x12h\x14)\x81\x94}\x94h\x17h\x12sbub\x8c\r_params_types\x94}\x94\x8c\twildcards\x94h\x06\x8c\tWildcards\x94\x93\x94)\x81\x94}\x94(h\r}\x94h\x0f]\x94(h\x11h\x12eh\x11h\x14)\x81\x94}\x94h\x17h\x11sbh\x12h\x14)\x81\x94}\x94h\x17h\x12sbub\x8c\x07threads\x94K\x01\x8c\tresources\x94h\x06\x8c\tResources\x94\x93\x94)\x81\x94(K\x01K\x01\x8c\x04/tmp\x94e}\x94(h\r}\x94(\x8c\x06_cores\x94K\x00N\x86\x94\x8c\x06_nodes\x94K\x01N\x86\x94\x8c\x06tmpdir\x94K\x02N\x86\x94uh\x0f]\x94(h\x11h\x12eh\x11h\x14)\x81\x94}\x94h\x17h\x11sbh\x12h\x14)\x81\x94}\x94h\x17h\x12sbhHK\x01hJK\x01hLhEub\x8c\x03log\x94h\x06\x8c\x03Log\x94\x93\x94)\x81\x94\x8c\x19../results/logs/dada2.log\x94a}\x94(h\r}\x94h\x0f]\x94(h\x11h\x12eh\x11h\x14)\x81\x94}\x94h\x17h\x11sbh\x12h\x14)\x81\x94}\x94h\x17h\x12sbub\x8c\x06config\x94}\x94\x8c\x04rule\x94\x8c\x05dada2\x94\x8c\x0fbench_iteration\x94N\x8c\tscriptdir\x94\x8cc/mnt/c/users/chuis/documents/biology_biosus/masterstage_2/barcode-phylogenetic-diversity/workflow/../bin\x94ub.');del script;from snakemake.logging import logger;from snakemake.script import snakemake; logger.printshellcmds = False;__real_file__ = __file__; __file__ = '/mnt/c/users/chuis/documents/biology_biosus/masterstage_2/barcode-phylogenetic-diversity/bin/run_tools.py';
######## snakemake preamble end #########
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
data_ids = []
for data_path in snakemake.input[1:]:
    with open(f'{data_path}', 'r') as data_f:
        data_id = data_f.read().split('\n')
    data_ids.extend(data_id)

## Get the right tool according to tool_name
# out_file variable is set for dealing with different export parameters (i.e. alpha/beta)
logging.info('Getting the tool and setting out_file variable.')

tool = gi.tools.get_tools(name=tool_name)
if not tool:
    raise ValueError(f"Tool '{tool_name}' not found in Galaxy instance. Check definition.")

tool_id = tool[0]['id']
out_file = os.path.basename(snakemake.output[0])
params = parameters.set(tool_name, data_ids, out_file)

## Parameter setting and running the tool
logging.info(f"Running the appropriate tool: {tool_name}")
try:
    results = gi.tools.run_tool(history['id'], tool_id, params)
    job = results['jobs'][0]
    jc = galaxy.jobs.JobsClient(galaxy_instance=gi)
    jc.wait_for_job(job_id=job['id'])

except ConnectionError as e:
    logging.error(f"Failed to run {tool_name}: {e}")

# Saving output ids. Tools with multiple outputs
logging.info(f"Saving output id(s) of {tool_name} to temporary file(s) for further use.")
for res in range(len(results['outputs'])):
    output_path = snakemake.output[res]
    result_id = results['outputs'][res]['id']
    with open(f'{output_path}', 'w') as output_f:
        output_f.write(result_id)
