# barcode-phylogenetic-diversity
## About the workflow
Repository for experiments with QIIME2-based alpha/beta phylogenetic diversity (PD) metrics in the Galaxy environment.
The goal of the described pipeline is to provide automated visualizations of PD metrics when provided with 
ITS metabarcoding data, in the context of the eDentity project at Naturalis Biodiversity Center. It is written 
in Python and R and implemented in Snakemake. It also makes use of conda for the management of the virtual environment
needed for the workflow and of the API 'bioblend' to access the Galaxy environment. The workflow is part of 
a master internship project of Ciel Huisman.

All rules of the workflow and their dependencies can be represented like this:
![DAG of entire workflow.](https://github.com/naturalis/barcode-phylogenetic-diversity/blob/main/full_dag.svg)

## Usage
### Installation and setup
1. Operating system: make sure you can access a Linux or WSL terminal, such as the Ubuntu (22.04.3) WSL that was used 
in making this pipeline and that the rest of the setup is based on.
2. Software: install conda and Python via your terminal, such as by installing 'miniforge' or 'miniconda'. At 
the time of building the workflow, conda was installed via 'mambaforge', but this is currently no longer supported 
(R need not be installed, as conda supports R packages). For miniforge, it should look something like this:
```
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash ./ Miniforge3-Linux-x86_64.sh
```
3. Virtual environment: use conda inside your terminal to create the virtual environment 
'barcode-phylogenetic-diversity' according the _environment.yml_ file in the root of this repository. Then use conda 
to activate this environment:
```
conda create barcode-phylogenetic-diversity -filename environment.yml
conda activate barcode-phylogenetic-diversity
```
4. API key: log into Galaxy and navigate to 'Users > Preferences > Manage API Key' and copy
your API key. Now, in your terminal, export it:
```
export GALAXY_API_KEY='pasteyourapikey'
```

### Running of the workflow 
1. Prepare the data directory: make sure the ITS sequences are present and are in the appropriate format, 
as described in the README.md in the data directory
2. Prepare the src directory: check/modify the tool list (src/tools.txt), parameters (src/parameters.py), 
and snakefile (src/snakefile), as described in the README.md in the src directory 
3. Navigate to the src directory in your terminal
4. Run snakemake. You must at least specify the number of cores to be used:

```
snakemake -c 1
```