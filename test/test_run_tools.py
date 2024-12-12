## Running the right tools with the right parameters

##Need more info in babysnake inputs to understand this
##input 1 is supposed to be tool_output of previous run
##input 2 is supposed to be list(tool_name, tool_id)
data_id = input(1)[0]
tool_name = input(2)[0]
tool_id = input(2)[1]

## Get the right parameters, data goes in here
if tool_name == 'qiime2 tools import':
    params = inputs().set('import_root|type', 'FeatureData__ob__Sequence__cb__') \
        .set('__q2galaxy__GUI__cond__format__|format', 'DNAFASTAFormat') \
        .set('__q2galaxy__GUI__cond__format__|data', dataset(data_id))

elif tool_name == 'qiime2 dada2 denoise-paired':
    params = input().set('', '')

elif tool_name == 'qiime2 phylogeny align-to-tree-mafft-raxml':
    params = input().set('', '')

elif tool_name == 'qiime2 feature-table rarefy':
    params = input().set('', '')

elif tool_name == 'qiime2 diversity alpha-phylogenetic':
    params = input().set('', '')

elif tool_name == 'qiime2 diversity beta-phylogenetic':
    params_beta = inputs().set('table', dataset()) \
        .set('phylogeny', dataset()) \
        .set('metric', 'unweighted_unifrac') \
        .set('__q2galaxy__GUI__section__extra_opts__|threads', 'auto') \
        .set('__q2galaxy__GUI__section__extra_opts__|variance_adjusted', False) \
        .set('__q2galaxy__GUI__section__extra_opts__|bypass_tips', False)

elif tool_name == 'qiime2 tools export':
    params = inputs().set('input', dataset(data_id)) \
        .set('input|type_peek', 'FeatureData__ob__Sequence__cb__') \
        .set('input|fmt_peek', 'FeatureData__ob__Sequence__cb__') \
        .set('output_format', 'DNAFASTAFormat')
else:
    print('No parameters for this tool. Check your tool list or definition of the tool.')


## Running the tool in question
try:
    results = gi.tools.run_tool(history['id'], tool_id, params)

    job = results['jobs'][0]
    jc = galaxy.jobs.JobsClient(galaxy_instance=gi)
    jc.wait_for_job(job_id=job['id'])

except ConnectionError as e:
    print(e.body)


##Writing output into id files
result_id = results['id'][0]
with open(f"../galaxy-ids/{tool_name}_out.txt", "w") as file:
    file.write(result_id.text)

#Clause for dada2 having two results (concept, confusing which is seq and which is tab)
if len(results) == 2:
    result_id2 = results['id'][1]
    with open(f"../galaxy-ids/{tool_name}_out2.txt", "w") as file:
        file.write(result_id2.text)


