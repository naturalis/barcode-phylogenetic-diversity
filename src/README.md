 ## Workflow details

### Overview

This workflow is based on five bin scripts and the three files inside this directory to operate correctly. Most rules 
deal with the running of tools within the Galaxy environment and rely on _bin/run_tools.py_ script, while the other 
rules and bin scripts facilitate this process and the visualization of the results.

### Framework: _Snakefile_

The snakefile describes the whole workflow through rules of input, shell/script, and output. Most rules 
have two types of input: a temporary text file with a tool name and another temporary text file with an id 
of the dataset to be used by the tool. The output consists of another text file with an id of the produced 
dataset, which is often used by another rule as input. These temporary files are stored in the 'results/galaxy-ids' 
directory, until no longer required.

### Tool name input: _tools.txt_

The tool name files are produced by the _bin/get_tool_names.py_ script, which uses the tool list in the _tools.txt_ 
file to do so. 

### Parameter setting: _parameters.py_

This script sets the parameters needed for the running of the tool. It is imported by the _bin/run_tools.py_
script and takes as arguments the aforementioned tool name as well as the id of the dataset, to produce 
the right parameter settings. The out_file variable is only needed when the same tool is used for different 
rules, which is in this case the _qiime2_export_alpha_ and _qiime2_export_alpha_ rules that both employ the 
'qiime2 tools export' tool.