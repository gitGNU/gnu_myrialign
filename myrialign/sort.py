
#
#    Copyright 2008 Paul Harrison
#
#    This file is part of Myrialign.
#    
#    Myrialign is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    
#    Myrialign is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with Myrialign.  If not, see <http://www.gnu.org/licenses/>.
#

import heapq

def strongly_connected_components(graph):
    result = [ ]
    stack = [ ]
    low = { }
        
    def visit(node):
        if node in low: return
	
	num = len(low)
        low[node] = num
        stack_pos = len(stack)
        stack.append(node)
	
        for successor in graph[node]:
            visit(successor)
            low[node] = min(low[node], low[successor])
        
        if num == low[node]:
	    component = tuple(stack[stack_pos:])
            del stack[stack_pos:]
            result.append(component)
	    for item in component:
	        low[item] = len(graph)
    
    for node in graph:
        visit(node)
    
    return result


def topological_sort(graph, priority=lambda x:x):
    count = { }
    for node in graph:
        count[node] = 0
    for node in graph:
        for successor in graph[node]:
            count[successor] += 1

    ready = [ (priority(node), node) for node in graph if count[node] == 0 ]
    heapq.heapify(ready)
    
    result = [ ]
    while ready:
        node = heapq.heappop(ready)[1]
        result.append(node)
        
        for successor in graph[node]:
            count[successor] -= 1
            if count[successor] == 0:
                heapq.heappush(ready, (priority(successor), successor))
    
    return result


def robust_topological_sort(graph, priority=lambda x:x):
    components = strongly_connected_components(graph)

    node_component = { }
    for component in components:
        for node in component:
            node_component[node] = component

    component_graph = { }
    for component in components:
        component_graph[component] = [ ]
    
    for node in graph:
        node_c = node_component[node]
        for successor in graph[node]:
            successor_c = node_component[successor]
            if node_c != successor_c:
                component_graph[node_c].append(successor_c) 

    return topological_sort(component_graph, priority)


if __name__ == '__main__':
    print robust_topological_sort({
        0 : [1],
        1 : [2],
        2 : [1],
    })


