## Running the right tools with the right parameters
## Depends on tool_import and params_set
## tool_name should be defined differently
tool_name = 'qiime2 tools import'

## Get the right parameters
params = params_dict[tool_name]

## Extract the right tool id
for entry in range(0, len(imported_tools)):
    if imported_tools[entry]['name'] == tool_name:
        tool_id = imported_tools[entry]['id']

## Running the tool in question
try:
    results = gi.tools.run_tool(history['id'], tool_id, params)

    job = results['jobs'][0]
    jc = galaxy.jobs.JobsClient(galaxy_instance=gi)
    jc.wait_for_job(job_id=job['id'])

except ConnectionError as e:
    print(e.body)

result_id =