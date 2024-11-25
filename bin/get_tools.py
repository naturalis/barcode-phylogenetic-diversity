## Imports all tools at once and puts them in a list used in params_set and tool_run
#Maybe make this into a dictionairy with tool_ids
tools = open('bin/tools.txt', 'r')
tool_list = tools.read().split('\n')

imported_tools = []
for tool in range(0, len(tool_list)):
    added_tool = gi.tools.get_tools(name=tool_list[tool])[0]
    imported_tools.append(added_tool)

tools_dict = dict(zip(tool_list, imported_tools))