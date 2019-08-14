# example_impl.py
# An example implementation of the Danseiji map projections.

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


CSV_DIR = 'C:/Users/jkunimune.522MT32/Documents/GitHub/Rubber-Earth/output/'
FILENAME = 'danseijiO30.csv'

with open(CSV_DIR+FILENAME, 'r') as f:
	data = csv.reader(f)
	N, m_1, m_2, l, o, p, w, h = next(data)

	for i in range(N):
		pass 

	for i in range(m_1):
		for j in range(m_2):
			pass

	for i in range(l):
		pass

	for i in range(o):
		for j in range(p):
			pass

	try:
		next(data)
	except StopIteration:
		print("CSV loaded successfully.")
	else:
		raise ValueError("CSV file did not terminate at the expected time.")