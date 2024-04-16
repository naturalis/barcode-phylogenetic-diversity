## Imports all tools at once and puts them in a list used in params_set and tool_run
tools = open('src/tools.txt', 'r')
tool_list = tools.read().split('\n')

imported_tools = []
for tool in range(0, len(tool_list)):
    added_tool = gi.tools.get_tools(name=tool_list[tool])[0]
    imported_tools.append(added_tool)
