# Script for obtaining tool names from tools.txt for use by run_tools.py and parameters.py.

# configure logging (for consistency)
import logging
logging.basicConfig(filename=snakemake.log[0], level=logging.DEBUG,
                    format='%(asctime)s:%(levelname)s:%(message)s')

# Splitting the tool list text file into multiple files.
logging.debug('Creating temporary files with tool names for further use.')
input_path = snakemake.input[0]
with open(f'{input_path}', 'r') as list_f:
    tool_list = list_f.read().split('\n')
    for name in tool_list:
        index = tool_list.index(name)
        tools_output = (snakemake.output[index])
        with open(f'{tools_output}', 'w') as name_f:
            name_f.write(name)
