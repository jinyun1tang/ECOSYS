import numpy as np

year_range = np.arange(1997, 2006)
file_names = []

for year in year_range:
    file_names.append('opt'+str(year))


for file_name in file_names:
    with open(file_name, 'r') as f:
        content = f.readlines()

    if content[4] == 'YES\n':
        content[4] = 'NO\n'
    if content[5] == 'YES\n':
        content[5] = 'NO\n'

    tmp = content[10].split(',')
    tmp[4] = '-1'

    content[10] = ','.join(tmp)

    with open(file_name, 'w') as f:
        f.write(''.join(content))
