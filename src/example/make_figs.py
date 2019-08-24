# example_impl.py
# An example implementation of the Danseiji map projections.
# Uses data from http://www.naturalearthdata.com/

# This is free and unencumbered software released into the public domain.

# Anyone is free to copy, modify, publish, use, compile, sell, or
# distribute this software, either in source code form or as a compiled
# binary, for any purpose, commercial or non-commercial, and by any
# means.

# In jurisdictions that recognize copyright laws, the author or authors
# of this software dedicate any and all copyright interest in the
# software to the public domain. We make this dedication for the benefit
# of the public at large and to the detriment of our heirs and
# successors. We intend this dedication to be an overt act of
# relinquishment in perpetuity of all present and future rights to this
# software under copyright law.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
# OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
# ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
# OTHER DEALINGS IN THE SOFTWARE.

# For more information, please refer to <http://unlicense.org>

import csv
import math
import matplotlib.pyplot as plt
import shapefile


CSV_DIR = '../../output/'
CSV_NAME = 'danseijiO3.csv'
SHP_DIR = '../../data/'
SHP_NAME = ['ne_110m_admin_0_countries', 'ne_110m_graticules_15']

SHOW_MESH = False

nodes = []    # a list of nodes: (ɸ, θ, x, y)
elements = [] # a list of elements: (n0, n1, n2) where n is the index of a node in the nodes list
cells = []    # a table of cells: (eE, eN, eW, eS) where e is the index of an element in the elements list
edge = []     # a list of node indices
table = []    # a table of interpolation points: (ɸ, θ)

with open(CSV_DIR+CSV_NAME, 'r') as f:
	data = csv.reader(f)
	N, m_1, m_2, l, o, p, w, h = [float(x) for x in next(data)]

	for i in range(int(N)): # load node locations
		x, y = [float(x) for x in next(data)]
		nodes.append((None, None, x, y)) 

	for i in range(int(m_1)): # load cell vertexen
		cells.append([])
		for j in range(int(m_2)):
			kind, *vertexen = [int(x) for x in next(data)]
			if kind == 0: # polar cells
				node_positions = [(0,1), (0,0), (1,0), (1,1)]
				if vertexen[0] == vertexen[1]: # north polar
					elements.append([vertexen[0], vertexen[2], vertexen[3]])
				else: # south polar
					elements.append([vertexen[0], vertexen[1], vertexen[2]])
				element_idx = len(elements) - 1
				cells[i].append([element_idx, element_idx, element_idx, element_idx])
			elif kind == 1: # NW-SE cells
				node_positions = [(0,1), (0,1), (0,0), (1,0), (1,0), (1,1)]
				elements.append([vertexen[1], vertexen[2], vertexen[3]])
				elements.append([vertexen[0], vertexen[4], vertexen[5]])
				northwest_idx = len(elements) - 2
				southeast_idx = len(elements) - 1
				cells[i].append([southeast_idx, northwest_idx, northwest_idx, southeast_idx])
			elif kind == -1: # NE-SW cells
				node_positions = [(0,1), (0,0), (0,0), (1,0), (1,1), (1,1)]
				elements.append([vertexen[0], vertexen[1], vertexen[5]])
				elements.append([vertexen[2], vertexen[3], vertexen[4]])
				northeast_idx = len(elements) - 2
				southwest_idx = len(elements) - 1
				cells[i].append([northeast_idx, northeast_idx, southwest_idx, southwest_idx])

			for node_idx, position in zip(vertexen, node_positions): # infer node spherical coordinates
				ɸ = 90 - (i+position[0])*180/m_1
				θ = (j+position[1])*360/m_2 - 180
				if nodes[node_idx][0] is not None:
					assert nodes[node_idx][0] == ɸ, "I thought for sure node {} was at {:.3f},{:.3f}".format(node_idx, ɸ, θ)
				nodes[node_idx] = (ɸ, θ, *nodes[node_idx][2:])

	for i in range(int(l)): # load edge
		node = [int(x) for x in next(data)]
		edge.append(node)

	for i in range(int(o)): # load table
		table.append([])
		for j in range(int(p)):
			φ, θ = [float(x) for x in next(data)]
			table[i].append((φ, θ))

	try: # assert that the entire file has been read
		next(data)
	except StopIteration:
		print("CSV loaded successfully.")
	else:
		raise ValueError("CSV file did not terminate at the expected time.")

plt.figure() # now start plotting

if SHOW_MESH:
	for element in elements:
		xs = [nodes[node_idx][2] for node_idx in element]
		ys = [nodes[node_idx][3] for node_idx in element]
		plt.fill(xs, ys, edgecolor='k', linewidth=1, fill=False) # plot the edges of the elements if desired

for shapefilename in SHP_NAME:
	sf = shapefile.Reader(SHP_DIR+shapefilename) # map actual coordinates onto the mesh
	for shape in sf.shapes():
		for k, part in enumerate(shape.parts):
			start = shape.parts[k]
			stop = shape.parts[k+1] if k+1 < len(shape.parts) else len(shape.points)
			xs, ys = [], []
			for θ, φ in shape.points[start:stop]:
				θ = min(θ, 179.999) # (these lines are kind of janky, but I don't like having to deal with if statements later)
				φ = max(φ, -89.999)
				i = (90-φ)*m_1/180 # use the coordinates to select the correct cell
				j = (θ+180)*m_2/360
				cell = cells[int(i)][int(j)]
				if i%1 <= j%1: # and the correct element
					if 1 - i%1 <= j%1:
						element = elements[cell[0]] # (E)
					else:
						element = elements[cell[1]] # (N)
				else:
					if i%1 <= 1 - j%1:
						element = elements[cell[2]] # (W)
					else:
						element = elements[cell[3]] # (S)

				φ_1, θ_1, x_1, y_1 = nodes[element[0]]
				φ_2, θ_2, x_2, y_2 = nodes[element[1]]
				φ_3, θ_3, x_3, y_3 = nodes[element[2]]
				det = (φ_2 - φ_3)*(θ_1 - θ_3) - (θ_2 - θ_3)*(φ_1 - φ_3)
				weight_1 = ((φ_2 - φ_3)*(θ - θ_3) - (θ_2 - θ_3)*(φ - φ_3))/det
				weight_2 = ((φ_3 - φ_1)*(θ - θ_3) - (θ_3 - θ_1)*(φ - φ_3))/det
				weight_3 = 1 - weight_1 - weight_2
				xs.append(weight_1*x_1 + weight_2*x_2 + weight_3*x_3) # finally, interpolate
				ys.append(weight_1*y_1 + weight_2*y_2 + weight_3*y_3)
			plt.plot(xs, ys, color='k', linewidth=1) # plot the shape

plt.axis('equal')
plt.axis('off')
plt.show()
# plt.savefig("test.svg", bbox_inches=0)