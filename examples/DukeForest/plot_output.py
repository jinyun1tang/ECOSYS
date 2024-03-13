import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors


folders = ['/clusterfs/lawrencium/yaningl/Research/DK3_duke_forest_1996_2005_041417/sim_'
           + str(i) for i in range(1, 11)]
year = 2005
col_num = 4

colors = ['red', 'green', 'blue', 'm', 'purple', 'brown', 'teal',
          'cyan', 'maroon', 'gold']

data = []

for folder in folders:
    data_tmp = np.loadtxt(folder+'/01010'+str(year)+'hc', skiprows=1,
                          usecols=col_num)
    data.append(data_tmp)

data_length = np.min([len(d) for d in data])

plt.figure()
for i in range(len(data)):
    plt.subplot(3,4,i+1)
    h = plt.plot(range(1, data_length+1), data[i][:data_length], label='case '+str(i+1),
                 color=colors[i])
    plt.legend(handles=h)

plt.savefig('plot.png')
