## Set parameters (
## Depends on tool_import
#qiime2 tools import, example for moving pictures data
params_imp = inputs().set('import_root|type','FeatureData__ob__Sequence__cb__')\
    .set('__q2galaxy__GUI__cond__format__|format', 'DNAFASTAFormat')\
    .set('__q2galaxy__GUI__cond__format__|data', dataset(seqs_id))

##qiime2 dada2 denoise-paired
params_dada = input().set('', '')

##qiime2 phylogeny align-to-tree-mafft-raxml
params_rax = input().set('', '')

##qiime2 feature-table rarefy
params_ftr = input().set('', '')

##qiime2 diversity alpha-phylogenetic
params_alpha = input().set('', '')

##qiime2 diversity beta-phylogenetic, example
params_beta = inputs().set('table', dataset())\
    .set('phylogeny', dataset())\
    .set('metric', 'unweighted_unifrac')\
    .set('__q2galaxy__GUI__section__extra_opts__|threads', 'auto')\
    .set('__q2galaxy__GUI__section__extra_opts__|variance_adjusted', False)\
    .set('__q2galaxy__GUI__section__extra_opts__|bypass_tips', False)

##qiime2 tools export, example from moving pictures data
params_exp = inputs().set('input', dataset(imp_id))\
    .set('input|type_peek', 'FeatureData__ob__Sequence__cb__')\
    .set('input|fmt_peek', 'FeatureData__ob__Sequence__cb__')\
    .set('output_format', 'DNAFASTAFormat')

## dummy parameters
params_imp = 0
params_dada = 1
params_rax = 2
params_ftr = 3
params_alpha = 4
params_beta = 5
params_exp = 6

params_list = [params_imp,
    params_dada,
    params_rax,
    params_ftr,
    params_alpha,
    params_beta,
    params_exp
]

params_dict = dict(zip(tool_list, params_list))