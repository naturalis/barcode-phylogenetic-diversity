#Parameter setting script.
# Needs dada2 parameters
# If you want to use different tools and/or parameters, copy the existing if/elif structure
# for the tool and/or parameters of your choosing from XML files in galaxy toolshed.


# Dependencies
from bioblend.galaxy.tools.inputs import inputs, dataset

# Put the parameters you want to use in here
def set(tool_name, data_ids, in_file):
    if tool_name == 'qiime2 tools import':
        params = inputs().set('import_root|type', 'FeatureData__ob__Sequence__cb__') \
            .set('__q2galaxy__GUI__cond__format__|format', 'DNAFASTAFormat') \
            .set('__q2galaxy__GUI__cond__format__|data', dataset(data_ids[0]))

    elif tool_name == 'qiime2 phylogeny align-to-tree-mafft-raxml':
        params = inputs().set('sequences', dataset(data_ids[0])) \
            .set('__q2galaxy__GUI__section__extra_opts__|substitution_model', 'GTRGAMMA') \
            .set('__q2galaxy__GUI__section__extra_opts__|raxml_version', 'Standard')

    elif tool_name == 'qiime2 feature-table rarefy':
        params = inputs().set('table', dataset(data_ids[0])) \
            .set('sampling_depth', 1103)

    elif tool_name == 'qiime2 diversity alpha-phylogenetic':
        params = inputs().set('table', dataset(data_ids[0])) \
            .set('phylogeny', dataset(data_ids[1])) \
            .set('metric', 'faith_pd')

    elif tool_name == 'qiime2 diversity beta-phylogenetic' :
        params = inputs().set('table', dataset(data_ids[0])) \
            .set('phylogeny', dataset(data_ids[1])) \
            .set('metric', 'weighted_unifrac')

    elif tool_name == 'qiime2 tools export' and 'alpha' in in_file:
        params = inputs().set('input', dataset(data_ids[0])) \
            .set('input|type_peek', 'SampleData__ob__AlphaDiversity__cb__') \
            .set('input|fmt_peek', 'AlphaDiversityDirectoryFormat') \

    elif tool_name == 'qiime2 tools export' and 'beta' in in_file:
        params = inputs().set('input', dataset(data_ids[0])) \
            .set('input|type_peek', 'DistanceMatrix') \
            .set('input|fmt_peek', 'DistanceMatrixDirectoryFormat')

#    elif tool_name == 'qiime2 tools export' and 'import' in in_file:
#            params = inputs().set('input', dataset(data_id)) \
#                .set('input|type_peek', 'FeatureData__ob__Sequence__cb__') \
#                .set('input|fmt_peek', 'FeatureData__ob__Sequence__cb__') \
#                .set('output_format', 'DNAFASTAFormat')

#    elif tool_name == 'qiime2 tools export' and 'mafft-raxml' in in_file:
#        params = inputs().set('input', dataset(data_ids[0])) \
#            .set('input|type_peek', 'Phylogeny__ob__Rooted__cb__') \
#            .set('input|fmt_peek', 'Phylogeny__ob__Rooted__cb__') \
#            .set('output_format', 'NewickDirectoryFormat')
    else:
        print('No parameters for this tool. Check tool list, or definition of tool_name or in_file.')

    return params
