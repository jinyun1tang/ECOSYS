from mpi4py import MPI
import glob
import simulator as sm
import os
import numpy as np
import sys


# Note rank 0 is not used for simulations. It causes problems with the subprocess module.
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nproc = comm.Get_size()
# print('nproc = ', nproc)

sim_name = '/global/home/users/yaningl/repos/ecosys_duke_run/ecosys.x'
sim_folder = '/global/scratch/yaningl/Research/Ecosys_DukeForest/ecosys_duke_forest_base'
base_file = 'ecosys_parameters.dat'
ecosys_file = 'rundk3_restart_1996_1997_input'
nparam = 4
variable = ['ECO_CO2_FLUX']
year = ['1997']


data_file = glob.glob(os.path.join(sim_folder, '*'))
data_file = [os.path.basename(df) for df in data_file]
data_file.remove(base_file)

var_list = []
for i in range(0, nparam+1):
    var_list.append('#parameter0' + str(i+1) + '#')

sim_ecosys_DF = sm.sim_ecosys_DukeForest(sim_name, base_file,
                                         var_list, ecosys_file, data_file,
                                         sim_folder, variable, year)

if rank == 0:
    nsample_per_param = 3
    nsample = nsample_per_param**nparam
    lower_end = 0.1
    upper_end = 10.0
    # last parameter value
    lpv = 1.0
    sz = (upper_end-lower_end)/(nsample_per_param-1)
    
    param_idx = 0
    parameters = np.zeros((nsample, nparam+1))

    for i1 in range(nsample_per_param):
        for i2 in range(nsample_per_param):
            for i3 in range(nsample_per_param):
                for i4 in range(nsample_per_param):
                    parameters[param_idx,:] = [lower_end+i1*sz, lower_end+i2*sz,
                                               lower_end+i3*sz, lower_end+i4*sz, lpv]
                    param_idx = param_idx+1

    nsim = param_idx
    
if rank == 0:    
    for i in range(1,nproc):
        comm.send(nsim, dest=i)
if rank != 0:
    nsim = comm.recv(source=0)

comm.Barrier()

# Note the -1's below is to skip rank 0 (no simulation on rank 0)
index_local = list(range(rank-1, nsim, nproc-1))

if rank == 0:
    # print('rank {0}: Sending data to the other ranks'.format(rank))
    # sys.stdout.flush()
    # parameters_local = parameters[index_local, :]
    
    # So if there are too many processors, we do not need to use the excessive ones
    for j in range(1,np.min([nproc, nsim+1])):
        comm.send(parameters[list(range(j-1,nsim,nproc-1)),:], dest=j)

if rank != 0 and rank < nsim+1:
    parameters_local = comm.recv(source=0)
    
comm.Barrier()
# if rank != 0:
#     print('rank {0}: Data have been received'.format(rank))

# print('rank {0}: index_local is '.format(rank), index_local)

if rank != 0 and rank < nsim+1:
    postfix = [str(i) for i in index_local]
    sim_ecosys_DF.create_files(parameters_local, postfix)
    sim_ecosys_DF.run_serial()

    # now get objective functions on each cpu
    sim_out_local, sim_date_local = sim_ecosys_DF.output()

    obs_file = 'observations.dat'
    output_file = 'output.dat'
    std_dev_file = 'std_dev.dat'

    obs = np.loadtxt(obs_file)
    output = np.loadtxt(output_file)
    lik_std = np.loadtxt(std_dev_file)
    # Prior values for x
    x_prior = 2.0*np.ones(nparam+1, )
    # Prior standard deviations for x
    x_std = np.ones(nparam+1, )

    objfunc_local = sim_ecosys_DF.objfunc_nlposterior(parameters_local,
                                                      sim_out_local, x_prior,
                                                      x_std, obs, lik_std)

    comm.send(objfunc_local, dest=0)

comm.Barrier()

if rank == 0:
    objfunc_all = np.zeros(nsim, )
    for i in range(1, np.min([nproc, nsim+1])):
        recved = comm.recv(source=i)
        objfunc_all[list(range(i-1,nsim,nproc-1))] = recved
    idx_max = np.argmax(objfunc_all)
    objfunc_max = objfunc_all[idx_max]
    print('maximum objective function value is: ', objfunc_max)
    print('maximum objective function value is obtained in the simulation folder : ',
          idx_max)
    with open(os.getcwd()+'/sim_'+str(idx_max)+'/'+base_file, 'r') as fh:
        param_opt = [float(s) for s in fh.readline().split()]

    print('optimal parameters after grid search are: ', param_opt)
