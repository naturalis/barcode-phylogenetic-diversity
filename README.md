# barcode-phylogenetic-diversity
This is a master internship project of Ciel Huisman at Naturalis Biodiversity Center. 
Repository for experiments with QIIME-based alpha/beta phylogenetic diversity metrics in the Galaxy environment.
The aim of this workflow is to provide visualizations of such PD metrics when provided with ITS metabarcoding data.

**Preparation**

- Make sure the ITS sequences are FASTQ files in demultiplexed Phred33 format and present in the 'data' folder
- Check/modify the tool list (src/tools.txt) and parameters (src/parameters.py)
- Check/modify the snakefile (src/snakefile), specifically the target rule for the correct output

**Executing the workflow**

1. Open your terminal and export your API key for access to Galaxy
2. Conda create and activate the virtual environment 'barcode-phylogenetic-diversity' with environment.yml
3. Navigate to 'src' directory
4. Run snakemake