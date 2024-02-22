# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 13:57:12 2024

@author: chuis
"""
## This is a script adjusted from bioblend, meant to run a qiime2 pipeline with Galaxy API. Probably needs to split.
## Based on: https://github.com/naturalis/bioblend-test/blob/main/client.py
## Notes with one hashtag are from bioblend by rvosa, two are mine

##Dependencies
from bioblend import galaxy
from bioblend.galaxy.tools.inputs import inputs, dataset
import os
import requests
from dotenv import load_dotenv, find_dotenv

# Endpoints for the Naturalis production instance. Domain can also be an IP address.
domain = 'galaxy.naturalis.nl'
dlbase = f'https://{domain}/api/datasets'

# The API key can be obtained from https://galaxy.naturalis.nl/user/api_key
load_dotenv(find_dotenv(usecwd=True))
api_key = os.environ.get('GALAXY_API_KEY')
gi = galaxy.GalaxyInstance(domain, key=api_key)

# Get or create a new history. The civilized thing would be to delete this when done.
histories = gi.histories.get_histories(name='qiime-dentity')
if len(histories) == 0:
    history = gi.histories.create_history('qiime-dentity')
else:
    history = histories[0]

## Upload fasta file into history (can be many formats, consult qiime documentation):
## https://docs.qiime2.org/2024.2/tutorials/importing/
## Trees (.tre) can also be imported or produced in the process.
with open('qiime_ready.fa', 'r') as file:
    fasta_contents = file.read()
fasta_id = gi.tools.paste_content(fasta_contents, history['id'], file_type='fa')['outputs'][0]['id']

## Files must be converted to .qza format using the import tool from qiime.
## The fasta file will become a file with semantics according to input data: SampleData[...].
q2_import = gi.tools.get_tools(name='qiime2 tools import')[0]
params_fa = inputs().set('infile', dataset(fasta_id))
gi.tools.run_tool(history['id'], qiime2_import['id'], params_fa)

## These are demultiplexed (demux) and denoised (dada2), resulting in FeatureTable[Frequency] and FeatureData[Sequence] files.
## Maybe this only works when they are added like last time
demux = gi.tools.get_tools(name='qiime2 demux emp-single')[0] #or emp-paired, depends on data
dada2 = gi.tools.get_tools(name='qiime2 dada2 denoise-single')[0] #or paired, or pyro, again

## Then they can be dereplicated, clustered and chimera filtered (vsearch) without change in semantics.
vsearch_derep = gi.tools.get_tools(name='qiime2 vsearch dereplicate-sequences')[0]
vsearch_cluster = gi.tools.get_tools(name='qiime2 vsearch cluster-features-open-reference')[0]
vsearch_chimfil = gi.tools.get_tools(name='qiime2 vsearch uchime-denovo')[0] #denovo seems standard?

## The next tools are for more chimera filtering (?) and abundance filtering of both (feature_table filter)
## Also not changing the semantics.
freq_filter = gi.tools.get_tools(name='qiime2 feature-table filter-features')[0]
seq_filter = gi.tools.get_tools(name='qiime2 feature-table filter-seqs')[0]

## Classification of sequences by vsearch or blast (or scikitlearn package from Python?) creating a FeatureTable[Taxonomy] file.
## The latter two tools are used to filter sequences and feature table alike (taxa filter seqs/table)
vsearch_classify = gi.tools.get_tools(name='qiime2 feature-classifier classify-consensus-vsearch')[0] #or *-consensus-blast
taxa_tabfil = gi.tools.get_tools(name='qiime2 taxa filter-table')[0]
taxa_seqfil = gi.tools.get_tools(name='qiime2 taxa filter-seqs')[0]

## Alignment of the FeatureData[Sequence] and phylogeny building can be done in one step by the alignment tool mafft
## combined with the mask tool to remove ambiguously aligned sequences, and a phylogeny builder such as fasttree or raxml).
## The alignment and masking result in FeatureData[AlignedSequence] files. The phylogeny building
## results in a Phylogeny[Unrooted] file that is converted it into Phylogeny[Rooted] by the midpoint-root tool.
## The FeatureTable[Frequency] can then be phylogenetically filtered for further analysis (phylogeny filter-table).
mafft_fasttree = gi.tools.get_tools(name='qiime2 phylogeny align-to-tree-mafft-fasttree')[0] #or *-mafft-raxml
phyl_midpoint = gi.tools.get_tools(name='qiime2 phylogeny midpoint-root')[0]
phyl_tabfil = gi.tools.get_tools(name='qiime2 phylogeny filter-table')[0]

## Some of which are alpha and beta diversity, which first need a rarefied FeatureTable[Frequency] (feature table-rarify).
## and a Phylogeny[Rooted]. Then the qiime tools for alpha and beta diversity can be fetched, the first resulting in
## a SampleData[AlphaDiversity] object and the latter in a DistanceMatrix file.
rarefy_tab = gi.tools.get_tools(name='qiime2 feature table-rarefy')[0]
alpha_div = gi.tools.get_tools(name='qiime2 diversity alpha-phylogenetic')[0]
beta_div = gi.tools.get_tools(name='qiime2 diversity beta-phylogenetic')[0]

## Alpha and beta diversity can be visualized with the following.
## Could not find emporer plots or other pcoa based visualizations on galaxy.
alpha_vis = gi.tools.get_tools(name='qiime2 diversity alpha-correlation')[0] #or alpha-group-significance
beta_vis = gi.tools.get_tools(name='qiime2 diversity beta-correlation')[0]


## Leftovers
params_beta = inputs().set('table', dataset(FeatureTable[Frequency]))\
    .set('phylogeny', dataset(Phylogeny[Rooted]))\
    .set('metric', 'unweighted_unifrac')\ ##Multiple options here
    .set('__q2galaxy__GUI__section__extra_opts__|threads', 'auto')\ ##different identation in xml
    .set('__q2galaxy__GUI__section__extra_opts__|variance_adjusted', False)\
    .set('__q2galaxy__GUI__section__extra_opts__|bypass_tips', False)\

try:

    # This output is written to outputs.json to figure out the structure.
    some_results = gi.tools.run_tool(history['id'], tool['id'], params)

    # We have to poll the job until it is done. For this we need the job ID.
    job = some_results['jobs'][0]
    jc = galaxy.jobs.JobsClient(galaxy_instance=gi)
    jc.wait_for_job(job_id=job['id'])

    # In this case there is a list of multiple outputs. The order seems to be unpredictable, so we have
    # to probe what we're looking for.
    for output in some_results['outputs']:
        if output['output_name'] == 'result':

            # The download URL should be composed as follows:
            # https://galaxy.naturalis.nl/api/datasets/8baaa5f7fc118c6e/display?to_ext=nhx
            url = f'{dlbase}/{output["id"]}/display?to_ext=nhx'
            response = requests.get(url)
            with open("qiime_output.nhx", "w") as file:
                file.write(response.text)

except ConnectionError as e:
    print(e.body)



