# barcode-phylogenetic-diversity
This is a master internship project of Ciel Huisman at Naturalis Biodiversity Center. 
Repository for experiments with QIIME-based alpha/beta phylogenetic diversity metrics in the Galaxy environment.
The aim of this workflow is to provide visualizations of such PD metrics when provided with ITS metabarcoding data.
Current example data comes from the qiime2 tutorial 'Moving Pictures'.

**Preparation**

- Make sure the ITS sequences are FASTQ files in Casava 1.8 demultiplexed format and present in the 'data' folder
- Check/modify the tool list (src/tools.txt) and parameters (src/parameters.py)
- Check/modify the snakefile (src/snakefile), specifically the target rule for the correct output

**Executing the workflow**

1. Open your terminal and export your API key for access to Galaxy
2. Create/conda activate the virtual environment 'barcode-phylogenetic-diversity' (will export yaml at some point)
3. Navigate to 'src' folder
4. Run snakemake

**Still under construction**

- dada2 tool is not implemented yet
- might get rid of the tool dictionary structure
- bin/galaxy_paste.py needs to handle multiple inputs
- bin/visualize.R needs tinkering and is not integrated yet
