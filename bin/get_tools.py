# Script for importing tools
# Might get rid of the whole tool dictionary structure.

# Dependencies
from bioblend import galaxy
import os
#import dotenv #test in python console
import logging


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


# Putting the tool names inside a list for the creation of a dictionary with tool ids
# Splitting the tool list text file into multiple files for run_tools.py input/parameters.py.
logging.debug('Creating temporary files with tool names for further use.')
input_path = snakemake.input[0]
with open(f'{input_path}', 'r') as list_f:
    tool_list = list_f.read().split('\n')
    for name in tool_list:
        index = tool_list.index(name)
        tools_output = (snakemake.output[index])
        with open(f'{tools_output}', 'w') as name_f:
            name_f.write(name)


# Make a tool id list to zip together with the tool names list
logging.debug('Getting tool ids. Creating temporary file with tool names and ids (dictionary).')
tool_ids = []
for tool in tool_list:
    added_tool = gi.tools.get_tools(name=tool)
    tool_id = added_tool[0]['id']
    tool_ids.append(tool_id)

# zip them together and write to file
tools_dict = dict(zip(tool_list, tool_ids))
ids_output = snakemake.output[-1]
with open(f'{ids_output}', 'w') as ids_f:
    print(tools_dict, file=ids_f)