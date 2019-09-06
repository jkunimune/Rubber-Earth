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
CSV_NAME = 'danseijiV40.csv'
SHP_DIR = '../../data/'
SHP_NAME = ['ne_110m_graticules_30', 'ne_110m_land']
OUT_DIR = '../../output/'
OUT_NAME = 'danseijiV.svg'

SHOW_MESH = False

LENGTH_THRESHOLD = .2

nodes = []    # a list of nodes: (x, y)
elements = [] # a list of elements: (nNE, nNW, nSW, nSE) where n is the index of a node in the nodes list
cells = []    # a table of cells: (eE, eN, eW, eS) where e is the index of an element in the elements list
boundary = []     # a list of node indices
table = []    # a table of interpolation points: (ɸ, θ)

with open(CSV_DIR+CSV_NAME, 'r') as f:
	data = csv.reader(f)
	N, m_1, m_2, l, o, p, w, h = [float(x) for x in next(data)]

	for i in range(int(N)): # load node locations
		x, y = [float(x) for x in next(data)]
		nodes.append((x, y)) 

	for i in range(int(m_1)): # load cell vertexen
		cells.append([])
		for j in range(int(m_2)):
			kind, *vertexen = [int(x) for x in next(data)]
			if kind == 0: # polar cells
				assert len(vertexen) == 4
				elements.append([vertexen[0], vertexen[1], vertexen[2], vertexen[3]])
				element_idx = len(elements) - 1
				cells[i].append([element_idx, element_idx, element_idx, element_idx])
			elif kind == 1: # NW-SE cells
				assert len(vertexen) == 6
				elements.append([vertexen[1], vertexen[2], vertexen[3],        None])
				elements.append([vertexen[0],        None, vertexen[4], vertexen[5]])
				northwest_idx = len(elements) - 2
				southeast_idx = len(elements) - 1
				cells[i].append([southeast_idx, northwest_idx, northwest_idx, southeast_idx])
			elif kind == -1: # NE-SW cells
				assert len(vertexen) == 6
				elements.append([vertexen[0], vertexen[1],        None, vertexen[5]])
				elements.append([       None, vertexen[2], vertexen[3], vertexen[4]])
				northeast_idx = len(elements) - 2
				southwest_idx = len(elements) - 1
				cells[i].append([northeast_idx, northeast_idx, southwest_idx, southwest_idx])

	for i in range(int(l)): # load boundary
		node, = [int(x) for x in next(data)]
		boundary.append(node)

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
		xs = [nodes[node_idx][0] for node_idx in element if node_idx is not None]
		ys = [nodes[node_idx][1] for node_idx in element if node_idx is not None]
		plt.fill(xs, ys, edgecolor='k', linewidth=.5, fill=False) # plot the edges of the elements if desired

xs = [nodes[node_idx][0] for node_idx in boundary]
ys = [nodes[node_idx][1] for node_idx in boundary]
plt.fill(xs, ys, edgecolor='k', linewidth=1, fill=False) # plot the boundary of the map

for element in elements: # extrapolate virtual nodes
	for i in range(4):
		if element[i] is None:
			node_1 = nodes[element[(i+1)%4]]
			node_2 = nodes[element[(i+2)%4]]
			node_3 = nodes[element[(i+3)%4]]
			nodes.append((node_1[0] - node_2[0] + node_3[0], node_1[1] - node_2[1] + node_3[1]))
			element[i] = len(nodes) - 1

for shapefilename in SHP_NAME:
	sf = shapefile.Reader(SHP_DIR+shapefilename) # map actual coordinates onto the mesh
	for shape in sf.shapes():
		for k, part in enumerate(shape.parts):
			start = shape.parts[k]
			stop = shape.parts[k+1] if k+1 < len(shape.parts) else len(shape.points)
			xs, ys = [], []
			for θ, φ in shape.points[start:stop] + [shape.points[start]]:
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

				x_NE, y_NE = nodes[element[0]]
				x_NW, y_NW = nodes[element[1]]
				x_SW, y_SW = nodes[element[2]]
				x_SE, y_SE = nodes[element[3]]
				x = i%1*(j%1*x_SE + (1-j%1)*x_SW) + (1-i%1)*(j%1*x_NE + (1-j%1)*x_NW)
				y = i%1*(j%1*y_SE + (1-j%1)*y_SW) + (1-i%1)*(j%1*y_NE + (1-j%1)*y_NW)

				if len(xs) > 0 and math.hypot(x - xs[-1], y - ys[-1]) < LENGTH_THRESHOLD: # if this line is short
					xs.append(x)
					ys.append(y) # add it on
				else: # if it is very long,
					plt.plot(xs, ys, color='k', linewidth=.5) # plot what we have and reset
					xs = [x]
					ys = [y]
			if len(xs) > 0:
				plt.plot(xs, ys, color='k', linewidth=.5)

plt.axis('equal')
plt.axis('off')
# plt.show()
plt.show()
