## Imports all tools at once and puts them in a list used in params_set and tool_run
#Maybe make this into a dictionairy with tool_ids
from bioblend import galaxy
import os

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

#getting the tool names, creating textfiles for each snakefile_concept rule
tools = open('src/tools.txt', 'r')
tool_list = tools.read().split('\n')

#==#


#Import tools
imported_tools = []
for tool in range(0, len(tool_list)):
    added_tool = gi.tools.get_tools(name=tool_list[tool])[0]
    imported_tools.append(added_tool)

#Create dictionary
tools_dict = dict(zip(tool_list, imported_tools))
#==#