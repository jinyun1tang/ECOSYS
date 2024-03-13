import numpy as np
import simulator as sm

year = ['2004', '2005']
nsim = 16
nparam = 11
files = []
obs_file = '/clusterfs/lawrencium/yaningl/Research/Ecosys_DukeForest/ecosys_duke_forest_weekly_data/observations.dat'
output_file = '/clusterfs/lawrencium/yaningl/Research/Ecosys_DukeForest/ecosys_duke_forest_weekly_data/output.dat'
for i in range(1, nsim+1):
    files.append('sim_'+str(i))

param = np.zeros((nsim, nparam))
nobs = sm.sim_ecosys_DukeForest.get_output(files[0], year).size
sim_output = np.zeros((nobs, nsim))
param_prior = 2.0*np.ones(nparam, )
param_mean = np.ones(nparam, )
prior_std = np.ones(nparam, )
obs = np.loadtxt(obs_file)
output = np.loadtxt(output_file)
lik_std = np.absolute(0.2*output)

cost_function = np.zeros(nsim, )
for i in range(nsim):
    print(i)
    param[i, :] = np.loadtxt(files[i]+'/ecosys_parameters.dat')    
    sim_output[:, i] = sm.sim_ecosys_DukeForest.get_output(files[i], year)

cost_function = sm.sim_ecosys_DukeForest.objfunc_nlposterior(param, sim_output, param_prior, 
                                                             prior_std, obs, lik_std)
print(cost_function)

cost_function_at_mean = 0
cost_function_at_mean = 0.5*np.square((obs-output)/lik_std).sum() + 0.5*np.square((param_prior-param_mean)/prior_std).sum()
print('cost_function_at_mean = ', cost_function_at_mean)
    
