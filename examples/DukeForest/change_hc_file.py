import numpy as np

row_range = np.arange(4, 102)
file_name = 'hc'

with open(file_name, 'r') as f:
    content = f.readlines()

for row in row_range:
    if content[row] == 'YES\n':
        content[row] = 'NO\n'

with open(file_name, 'w') as f:
    f.write(''.join(content))
