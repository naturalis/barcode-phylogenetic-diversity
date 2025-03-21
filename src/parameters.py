## Parameter setting script, imported by run_tools.py

# Dependencies
from bioblend.galaxy.tools.inputs import inputs, dataset

# Put the parameters you want to use in here
def set(tool_name, data_ids, out_file):
    if tool_name == 'qiime2 tools import':
        params = inputs().set('import_root|__q2galaxy__GUI__cond__format__|data', dataset(data_ids[0])) \
            .set('import_root|__q2galaxy__GUI__cond__format__|format', 'PairedEndFastqManifestPhred33') \
            .set('import_root|type', 'SampleData__ob__PairedEndSequencesWithQuality__cb__')

    elif tool_name == 'qiime2 dada2 denoise-paired':
        params = inputs().set('demultiplexed_seqs', dataset(data_ids[0])) \
            .set('trunc_len_f', 0) \
            .set('trunc_len_r', 0)

    elif tool_name == 'qiime2 phylogeny align-to-tree-mafft-raxml':
        params = inputs().set('sequences', dataset(data_ids[0])) \
            .set('__q2galaxy__GUI__section__extra_opts__|substitution_model', 'GTRGAMMA') \
            .set('__q2galaxy__GUI__section__extra_opts__|raxml_version', 'Standard')

    elif tool_name == 'qiime2 feature-table rarefy':
        params = inputs().set('table', dataset(data_ids[0])) \
            .set('sampling_depth', 51992)

    elif tool_name == 'qiime2 diversity alpha-phylogenetic':
        params = inputs().set('table', dataset(data_ids[0])) \
            .set('phylogeny', dataset(data_ids[1])) \
            .set('metric', 'faith_pd')

    elif tool_name == 'qiime2 diversity beta-phylogenetic':
        params = inputs().set('table', dataset(data_ids[0])) \
            .set('phylogeny', dataset(data_ids[1])) \
            .set('metric', 'weighted_unifrac')

    elif tool_name == 'qiime2 tools export' and 'alpha' in out_file:
        params = inputs().set('input', dataset(data_ids[0])) \
            .set('input|type_peek', 'SampleData__ob__AlphaDiversity__cb__') \
            .set('input|fmt_peek', 'AlphaDiversityDirectoryFormat') \

    elif tool_name == 'qiime2 tools export' and 'beta' in out_file:
        params = inputs().set('input', dataset(data_ids[0])) \
            .set('input|type_peek', 'DistanceMatrix') \
            .set('input|fmt_peek', 'DistanceMatrixDirectoryFormat')

    else:
        print('No parameters for this tool. Check tool list, or definition of tool_name or out_file.')

    return params
