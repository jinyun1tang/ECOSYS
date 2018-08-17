# -*- coding: utf-8 -*-
"""
goerge
simulator class

@author: Yaning Liu
"""
from __future__ import print_function
from abc import ABCMeta
import numpy as np
import os
import sys
import shutil
import re
import subprocess
import time
import glob
import csv
import pandas as pd


class simulator(object):
    """ A simulator superclass: all simulators are derived from this class

    Attributes:
        sim_name;
        prefix;
        id;
        nsim;
        base_file;
        dir_list;
        var_list;
        batch_file;
        qsub_command;
        batch_check;
        batch_size = 0;
        is_success;
        verbose = 0;
    """

    __metaclass__ = ABCMeta

    batch_size = 0
    verbose = 0
    id = []
    batch_file = ''

    def __init__(self, sim_name_in, base_file_in, var_list_in,
                 prefix_in='sim'):
        """initializer for the abstract base class

        :param sim_name_in: string, executable name, e.g., '../myit'
        :param base_file_in: string, input file, e.g., 'samROMffi'
        :param var_list_in: a list of string, the strings to be replaced in \
        input file
        :param prefix_in: string, the prefix to be added in the names of \
        simulation folders, e.g., 'sim'
        :returns: None
        :rtype: None

        """
        self.sim_name = sim_name_in  # a str
        self.base_file = base_file_in  # a str
        self.var_list = var_list_in  # a list of str
        self.prefix = prefix_in  # a str

    def create_files(self, param_in, postfix=[]):
        """the function does the following things:

        1) create all folders in the form of:
        if postfix=[]:
        current_dir/prefix_1, current_dir/prefix_2,...
        if postfix is not empty
        current_dir/prefix_postfix[0], current_dir/prefix_postfix[1], ...
        e.g., ~/sim_1, ~/sim_2,...

        2) copy base file to all folders. Example of the base file is samROMffi

        3) In each folder, replace the parameters in the base file
        using the corresponding parameters from param_in

        4) set self.dir_list to be the list of folder names

        :param param_in: 2d array (nsample, nparam)
        :param postfix: a list of string, the postfix to be appended to \
        simulation folders
        :returns: the list of folder names as class variables
        :rtype: a list of string

        """

        cur_dir = os.path.dirname(os.path.realpath(self.base_file))
        fpath = "{}/{}".format(cur_dir, self.base_file)
        status = os.path.isfile(fpath)

        if (not status):
            # sys.exit(self.base_file+' does not exist.')
            raise Exception(self.base_file+' does not exist.')

        # Check if the number of elements of postfix is the same as that of
        # if param_in is 1d array:
        if param_in.ndim == 1:
            param_in = param_in[np.newaxis, :]

        # param_in
        if postfix and len(postfix) != param_in.shape[0]:
            raise Exception('The number of postfixes does not equal that of '
                            'parameters.')

        fdir = []
        for i in range(0, param_in.shape[0]):
            if len(self.id) == 0:
                if not postfix:
                    fdir.append("{}/{}_{}".format(cur_dir, self.prefix,
                                                  str(i+1)))
                else:
                    fdir.append("{}/{}_{}".format(cur_dir, self.prefix,
                                                  postfix[i]))
            else:
                if not postfix:
                    fdir.append("{}/{}_{}_{}".format(cur_dir, self.prefix,
                                                     str(self.id), str(i+1)))
                else:
                    fdir.append("{}/{}_{}_{}".format(cur_dir, self.prefix,
                                                     str(self.id), postfix[i]))

            # status = os.path.isdir(fdir[i])
            # if status == False: os.mkdir(fdir[i])
            # else:
            if os.path.isdir(fdir[i]):
                shutil.rmtree(fdir[i])

            os.mkdir(fdir[i])

            flist = fdir[i] + '/' + self.base_file
            try:
                shutil.copyfile(self.base_file, flist)
            except IOError as e:
                # print(e.errno, e.strerror)
                # sys.exit('File copy failed: check ' + flist + '. ')
                raise Exception('File copy failed: check ' + flist + '. ' +
                                str(e))
            self.update_file(flist, self.var_list, param_in[i, :])

        self.dir_list = fdir

    def setup_dir_list(self, param_in, postfix=[]):
        """ Set up the names of the simulation folders for the parameters
        param_in

        :param param_in: numpy array of size (nsample, nparam), each row of \
        the array is a set of parameters
        :param postfix: a list of strings, containig the postfix (numbers) to \
        be appended to the created folders
        :returns: dir_list as class variable members
        :rtype: a list of strings

        """

        cur_dir = os.path.dirname(os.path.realpath(self.base_file))

        # Check if the number of elements of postfix is the same as that of
        # param_in
        if param_in.ndim == 1:
            param_in = param_in[np.newaxis, :]

        if postfix and len(postfix) != param_in.shape[0]:
            raise Exception('The number of postfixes does not equal that of '
                            'parameters.')

        fdir = []
        for i in range(0, param_in.shape[0]):
            if len(self.id) == 0:
                if not postfix:
                    fdir.append("{}/{}_{}".format(cur_dir, self.prefix,
                                                  str(i+1)))
                else:
                    fdir.append("{}/{}_{}".format(cur_dir, self.prefix,
                                                  postfix[i]))
            else:
                if not postfix:
                    fdir.append("{}/{}_{}_{}".format(cur_dir, self.prefix,
                                                     str(self.id), str(i+1)))
                else:
                    fdir.append("{}/{}_{}_{}".format(cur_dir, self.prefix,
                                                     str(self.id), postfix[i]))

        self.dir_list = fdir

    def run(self):
        """Copy files to simulation folders and run simulations

        :returns: None
        :rtype: None

        """
        command3 = 'cd ..'
        if self.batch_file:
            command2 = 'qsub ' + self.batch_file
            for i in range(0, len(self.dir_list)):
                try:
                    shutil.copyfile(self.batch_file, self.dir_list[i])
                except IOError as e:
                    print(e.errno, e.strerror)
                    sys.exit('Batch file copy failed. ')
                command1 = 'cd ' + self.dir_list[i]
                try:
                    os.chdir(self.dir_list[i])
                except OSError as e:
                    sys.exit(command1 + ' failed!')
                status = subprocess.call(command2.split())
                if status:
                    sys.exit(command2 + ' failed!')
                try:
                    os.chdir('..')
                except OSError as e:
                    sys.exit(command3 + ' failed!')
        else:
            command2 = self.sim_name + ' ' + self.base_file
            if not self.batch_size:
                for i in range(0, len(self.dir_list)):
                    command1 = 'cd ' + self.dir_list[i]
                    status = subprocess.call(command1.split())
                    if status:
                        sys.exit(command1 + ' failed!')
                    status = subprocess.call(command2.split())
                    if status:
                        sys.exit(command2 + ' failed!')
                    status = subprocess.call(command3.split())
                    if status:
                        sys.exit(command3 + ' failed!')
            else:
                idx = 0
                while idx <= len(self.dir_list)-1:
                    for i in range(idx, min(idx+self.batch_size,
                                            len(self.dir_list))):
                        command1 = 'cd ' + self.dir_list[i]
                        status = subprocess.call(command1.split())
                        if status:
                            sys.exit(command1 + ' failed!')
                        status = subprocess.call(command2.split())
                        if status:
                            sys.exit(command2 + ' failed!')
                        status = subprocess.call(command3.split())
                        if status:
                            sys.exit(command3 + ' failed!')
                    idx = i+1

    def obj_func_value(self, param_in, caldata, postfix=[]):
        """Compute objective function values. Note that here the objective
        function is just the sum((sim-obs)**2)

        :param param_in: 2d array of size (nsample, nparam), input parameters
        :param caldata: 1d array, the calibration/observation data
        :param postfix: a list of strings, containig the postfix (numbers) to \
        be appended to the created folders
        :returns: obj_out, the objective function values
        :rtype: 1d array of size (nsample, )

        """

        self.create_files(param_in, postfix)
        self.run()
        sim_out = self.output()
        obj_out = np.zeros(param_in.shape[0],)
        for i in range(0, param_in.shape[0]):
            obj_out[i] = np.square(caldata-sim_out[i]).sum()
        return obj_out

    def output(self):
        """Get the output

        :returns: sim_out
        :rtype: numpy array

        """
        if self.batch_file:
            cmd = 'check'
            while cmd:
                time.sleep(5)
                try:
                    cmd = subprocess.check_output(self.batch_check)
                except subprocess.CalledProcessError as e:
                    sys.exit(e.output)
        sim_out = []
        return sim_out

    @staticmethod
    def update_file(filename, var_list, param_in):
        """Replace the strings in var_list by parameter values given by param_in
        in the input file

        :param filename: string, the input file
        :param var_list: a list of string, the strings to be replaced in \
        input file
        :param param_in: 1d array of size (nparam, ), input parameters
        :returns: None
        :rtype: None

        """

        with open(filename, 'r') as fin:
            input = fin.read()

        for i in range(0, len(var_list)):
            numstr = '{:.13E}'.format(param_in[i])
            input = re.sub(var_list[i], numstr, input)

        with open(filename, 'w') as fout:
            fout.write(input)

class sim_ecosys_DukeForest(simulator):
    """The ecosys simulator for the Duke Forest problem
    """

    def __init__(self, sim_name, base_file,
                 var_list, ecosys_file, data_file, sim_folder,
                 variable, year, prefix='sim'):
        """Initializer for ecosys simulator class

        :param sim_name: string, executable name, including directory part
        e.g., '/global/home/users/yaningl/repos/ecosys_120316_for_dukeforest/ecosys.x'
        :param base_file: string, input file with paramters to change, \
        e.g., ecosys_parameters.dat
        :param var_list: a list of string, the strings to be replaced in \
        input file
        :param ecosys_file: string, ecosys input file, e.g., 'runba'
        :param data_file: a list of strings, the data file, \
        e.g., ['cond0.dat'] in itough2. For ecosys, they are all the files
        needed and that are to be copied to each simulation folder
        :param sim_folder: string, the folder with simulation files, in which the \
        input file, data file, ... are located,
        :param variable: list of str, the output variable, for example ['ECO_CO2_FLUX']
        :param year: a list of string, the year of which the observation data \
        will be used, for example ['2000'], ['2000', '2001']
        :param prefix: string, the prefix to be added in the names of \
        simulation folders, e.g., 'sim'
        :returns: None
        :rtype: None

        """
        super(sim_ecosys_DukeForest, self).__init__(sim_name, base_file,
                                                    var_list, prefix)
        self.sim_folder = sim_folder
        # print("Work Directory is '{}'".format(sim_folder_in))
        self.variable = variable  # e.g., 'ECO_CO2_FLUX'
        self.year = year
        self.var_list = var_list  # a list of str

        self.execlist = dict()
        self.execlist['data_file'] = data_file
        self.execlist['sim_name'] = sim_name  # a str
        self.execlist['base_file'] = base_file  # a str
        self.execlist['ecosys_file'] = ecosys_file  # a str
        self.execlist['prefix'] = prefix  # a str
        self.execlist['sim_folder'] = self.sim_folder

    def output(self):
        """Go to all simulation folders and get the output from
        ecosys simulations

        :returns: sim_out: a dictionary, the output of simulations.
                  sim_out['DOY']: 1d array of float, the date of year (0-365)
                  sim_out['DATE']: 1d array of str, Date
                  sim_out['output']: 2d array of float with size (noutput, nsim),
                                     the output
        :rtype: sim_out: a dictionary,

        """

        super(sim_ecosys_DukeForest, self).output()
        nsim = len(self.dir_list)
        sim_out = dict()
        sim_out['DOY'] = np.array([])
        sim_out['DATE'] = np.array([])
        sim_out['output'] = np.array([])

        i = 0
        hc_file_names = []
        hh_file_names = []
        # dc_file_names = []
        for ii in range(len(self.year)):
            hc_file_names.extend(glob.glob(self.dir_list[i] + '/' + '*0' +
                                        self.year[ii] + 'hc'))
            hh_file_names.extend(glob.glob(self.dir_list[i] + '/' + '*0' +
                                        self.year[ii] + 'hh'))
            # dc_file_names.extend(glob.glob(self.dir_list[i] + '/' + '*0' +
            #                             self.year[ii] + 'dc'))

        if len(hc_file_names) != len(self.year):
            raise Exception('hc files missing; simulation is not finished!')
        if len(hh_file_names) != len(self.year):
            raise Exception('hh files missing; simulation is not finished!')
        # if len(dc_file_names) != len(self.year):
        #     raise Exception('dc files missing; simulation is not finished!')

        hc_base_file_names = [os.path.basename(fn) for fn in hc_file_names]
        hh_base_file_names = [os.path.basename(fn) for fn in hh_file_names]
        # dc_base_file_names = [os.path.basename(fn) for fn in dc_file_names]

        for i in range(0, nsim):
            sim_out_tmp = np.array([])

            hc_file_names = [self.dir_list[i] + '/' + bfn
                             for bfn in hc_base_file_names]
            for fn in hc_file_names:
                if not os.path.isfile(fn):
                    raise Exception('No file ' + fn +
                                    '; simulation is not finished!')

                data_tmp = np.loadtxt(fn, skiprows=1)
                if i == 0:
                    sim_out['DOY'] = np.concatenate((sim_out['DOY'],
                                                     data_tmp[:, 0]))
                    sim_out['DATE'] = np.concatenate((sim_out['DATE'],
                                                      data_tmp[:, 1]))

                sim_out_tmp = np.concatenate((sim_out_tmp, data_tmp[:, 4]))

                # average the data every 7 days, that is every 7*24=168
                # sim_out_tmp= np.concatenate(
                #     (sim_out_tmp, np.mean(data_tmp[:, 4][:8736].reshape(-1, 168), 1)))

            hh_file_names = [self.dir_list[i] + '/' + bfn
                             for bfn in hh_base_file_names]
            for fn in hh_file_names:
                if not os.path.isfile(fn):
                    raise Exception('No file ' + fn +
                                    '; simulation is not finished!')

                data_tmp = np.loadtxt(fn, skiprows=1)
                if i == 0:
                    sim_out['DOY'] = np.concatenate((sim_out['DOY'],
                                                     data_tmp[:, 0]))
                    sim_out['DATE'] = np.concatenate((sim_out['DATE'],
                                                      data_tmp[:, 1]))

                sim_out_tmp = np.concatenate((sim_out_tmp, data_tmp[:, 11]))

                # average the data every 7 days, that is every 7*24=168
                # sim_out_tmp= np.concatenate(
                #     (sim_out_tmp, np.mean(data_tmp[:, 11][:8736].reshape(-1, 168), 1)))

            # dc_file_names = [self.dir_list[i] + '/' + bfn
            #                  for bfn in dc_base_file_names]
            # for fn in dc_file_names:
            #     if not os.path.isfile(fn):
            #         raise Exception('No file ' + fn +
            #                         '; simulation is not finished!')
            #
            #     data_tmp = np.loadtxt(fn, skiprows=1)
            #     if i == 0:
            #         sim_out['DOY'] = np.concatenate((sim_out['DOY'],
            #                                          data_tmp[:, 0]))
            #         sim_out['DATE'] = np.concatenate((sim_out['DATE'],
            #                                           data_tmp[:, 1]))
            #     sim_out_tmp = np.concatenate(
            #         (sim_out_tmp, data_tmp[:, 2]))

            if sim_out['output'].size == 0:
                sim_out['output'] = sim_out_tmp[:, np.newaxis]
            else:
                sim_out['output'] = np.concatenate((sim_out['output'],
                                                    sim_out_tmp[:, np.newaxis]),
                                                   axis=1)


        # indices = np.concatenate((range(0,2800), range(2900,8150), range(8250,8760)))
        indices = range(1000,1240)
        return sim_out['output'][indices,:], sim_out['DATE'][indices]
        # return sim_out['output'], sim_out['DATE']

    @staticmethod
    def get_output(foldername, year):
        """Go to the folder with name foldername and get the output from
        ecosys simulations

        :param foldername: string, the foldername
        :param year: a list of string, the year of which the observation data \
        will be used, for example ['2000'], ['2000', '2001']

        :returns: sim_out: a dictionary, the output of simulations.
                  sim_out['DOY']: 1d array of float, the date of year (0-365)
                  sim_out['DATE']: 1d array of str, Date
                  sim_out['HOUR']: 1d array of int, the hour (1-24)
                  sim_out['output']: 1d array of float, the output
        :rtype: sim_out: a dictionary

        """

        sim_out = dict()
        sim_out['DOY'] = np.array([])
        sim_out['DATE'] = np.array([])
        sim_out['output'] = np.array([])

        hc_file_names = []
        hh_file_names = []
        # dc_file_names = []
        for i in range(len(year)):
            hc_file_names.extend(glob.glob(foldername + '/' + '*0' + year[i] +
                                           'hc'))
            hh_file_names.extend(glob.glob(foldername + '/' + '*0' + year[i] +
                                           'hh'))
            # dc_file_names.extend(glob.glob(foldername + '/' + '*0' + year[i] +
            #                                'dc'))

        if len(hc_file_names) != len(year):
            raise Exception('hc files missing; simulation is not finished!')
        if len(hh_file_names) != len(year):
            raise Exception('hh files missing; simulation is not finished!')
        # if len(dc_file_names) != len(year):
        #     raise Exception('dc files missing; simulation is not finished!')

        for fn in hc_file_names:
            if not os.path.isfile(fn):
                raise Exception('No file ' + fn +
                                '; simulation is not finished!')
            data_tmp = np.loadtxt(fn, skiprows=1)

            sim_out['DOY'] = np.concatenate((sim_out['DOY'],
                                             data_tmp[:, 0]))
            sim_out['DATE'] = np.concatenate((sim_out['DATE'],
                                              data_tmp[:, 1]))

            sim_out['output'] = np.concatenate((sim_out['output'], data_tmp[:, 4]))

            # average the data every 7 days, that is every 7*24=168
            # sim_out['output'] = np.concatenate(
            #     (sim_out['output'], np.mean(data_tmp[:, 4][:8736].reshape(-1, 168), 1)))
        for fn in hh_file_names:
            if not os.path.isfile(fn):
                raise Exception('No file ' + fn +
                                '; simulation is not finished!')
            data_tmp = np.loadtxt(fn, skiprows=1)
            sim_out['DOY'] = np.concatenate((sim_out['DOY'],
                                             data_tmp[:, 0]))
            sim_out['DATE'] = np.concatenate((sim_out['DATE'],
                                              data_tmp[:, 1]))

            sim_out['output'] = np.concatenate((sim_out['output'], data_tmp[:, 11]))

            # average the data every 7 days, that is every 7*24=168
            # sim_out['output'] = np.concatenate(
            #     (sim_out['output'], np.mean(data_tmp[:, 11][:8736].reshape(-1, 168), 1)))

        # for fn in dc_file_names:
        #     if not os.path.isfile(fn):
        #         raise Exception('No file ' + fn +
        #                         '; simulation is not finished!')
        #     data_tmp = np.loadtxt(fn, skiprows=1)
        #     sim_out['DOY'] = np.concatenate((sim_out['DOY'],
        #                                      data_tmp[:, 0]))
        #     sim_out['DATE'] = np.concatenate((sim_out['DATE'],
        #                                       data_tmp[:, 1]))
        #     sim_out['output'] = np.concatenate(
        #         (sim_out['output'], data_tmp[:, 2]))

        # indices = np.concatenate((range(0,2800), range(2900,8150), range(8250,8760)))
        indices = range(1000,1240)
        return sim_out['output'][indices]
        # return sim_out['output']

    @staticmethod
    def objfunc_nlposterior(param, output, param_prior, prior_std, obs,
                            lik_std):
        """Compute the objective function values---negative log of posterior

        :param param: 2d array (nsample, nparam), input parameters
        :param output: 2d array (nobs=ntime*n_obs_pt, nsample), output values
        :param param_prior: 1d array (nparam), the prior parameters
        :param prior_std: 1d array (nparam), the std of the prior parameters
        :param obs: 1d array (nobs=ntime*n_obs_pt), the observations
        :param lik_std: 1d array (nobs=ntime*n_obs_pt), the Gaussian \
        likelihood standard deviations
        :returns: the negative log posterior objective function values
        :rtype: 1d array (nsample)

        """
        nsample = param.shape[0]
        objfunc = np.zeros(nsample, )
        for i in range(nsample):
            objfunc[i] = 0.5*np.square((obs-output[:, i])/lik_std).sum()
            objfunc[i] += 0.5*np.square(
                (param_prior-param[i, :])/prior_std).sum()

        return objfunc

    def create_files(self, param_in, postfix=[]):
        """the function does the following things:

        1) create all folders in the form of:
        if postfix=[]:
        current_dir/prefix_1, current_dir/prefix_2,...
        if postfix is not empty
        current_dir/prefix_postfix[0], current_dir/prefix_postfix[1], ...
        e.g., ~/sim_1, ~/sim_2,...

        3) In each folder, replace the parameters
        in the base file using the corresponding parameters from param_in

        4) set self.dir_list to be the list of folder names

        :param param_in: 2d array (nsample, nparam)
        :param postfix: a list of string, the postfix to be appended to \
        simulation folders
        :returns: the list of folder names as class variables
        :rtype: a list of string

        """

        # cur_dir = os.path.dirname(os.path.realpath(self.execlist['ecosys_file']))
        cur_dir = os.getcwd()

        # Check if the number of elements of postfix is the same as that of
        # if param_in is 1d array:
        if param_in.ndim == 1:
            param_in = param_in[np.newaxis, :]

        # param_in
        if postfix and len(postfix) != param_in.shape[0]:
            raise Exception('The number of postfixes does not equal that of '
                            'parameters.')

        fdir = []
        for i in range(0, param_in.shape[0]):
            if len(self.id) == 0:
                if not postfix:
                    fdir.append("{}/{}_{}".format(cur_dir, self.prefix,
                                                  str(i+1)))
                else:
                    fdir.append("{}/{}_{}".format(cur_dir, self.prefix,
                                                  postfix[i]))
            else:
                if not postfix:
                    fdir.append("{}/{}_{}_{}".format(cur_dir, self.prefix,
                                                     str(self.id), str(i+1)))
                else:
                    fdir.append("{}/{}_{}_{}".format(cur_dir, self.prefix,
                                                     str(self.id), postfix[i]))

            if os.path.isdir(fdir[i]):
                shutil.rmtree(fdir[i])

            os.mkdir(fdir[i])

            flist = fdir[i] + '/' + self.execlist['base_file']
            try:
                shutil.copyfile("{}/{}".format(self.execlist['sim_folder'],
                                           self.base_file), flist)
                # shutil.copyfile(self.base_file, flist)
            except IOError as e:
                # print(e.errno, e.strerror)
                # sys.exit('File copy failed: check ' + flist + '. ')
                raise Exception('File copy failed: check ' + flist + '. ' +
                                str(e))
            self.update_file(flist, self.var_list, param_in[i, :])

        self.dir_list = fdir

    @staticmethod
    def update_file(filename, var_list, param_in):
        """Replace the strings in var_list by parameter values given by param_in
        in the input file

        :param filename: string, the input file
        :param var_list: a list of string, the strings to be replaced in \
        input file. For ecosys, var_list is not used
        :param param_in: 1d array of size (nparam, ), input parameters
        :returns: None
        :rtype: None

        """

        with open(filename, 'r') as fin:
            input_content = fin.read()

        for i in range(0, len(var_list)):
            numstr = '{:.7E}'.format(param_in[i])
            input_content = re.sub(var_list[i], numstr, input_content)

        with open(filename, 'w') as fout:
            fout.write(input_content)

    def run_parallel(self, execution):
        """Run simulations in parallel with tigres. The function does:

        1) run run_ecosys_DukeForest_prep
        2) run all simulations in all folders

        :param execution: The tigres execution mode, \
        e.g., EXECUTION_LOCAL_THREAD
        :returns: None
        :rtype: None

        """
        import tigres

        self.execution = execution
        tigres.start(name='Run ecosys parallel for Duke Forest',
                     log_dest=os.path.splitext(__file__)[0]+'.log',
                     execution=execution)
        tigres.set_log_level(tigres.Level.INFO)

        # Create Tigres Task
        task_ecosys = tigres.Task('ecosys_task_prep', tigres.FUNCTION,
                                  impl_name=run_ecosys_DukeForest_prep,
                                  input_types=[dict, str])
        # Create Tigres Task array
        task_array = tigres.TaskArray('ecosys_tasks', tasks=[task_ecosys])
        # Create Tigres InputArray
        input_array = tigres.InputArray('folder_names', [])
        # Iterate through dir_list and add input values
        for i in range(0, len(self.dir_list)):
            input_values = tigres.InputValues('A folder {}'.format(
                self.dir_list[i]), [self.execlist, self.dir_list[i]])
            input_array.append(input_values)

        # Run the Tigres Parallel template for preparing the run tough
        # executables
        run_ecosys_commands = tigres.parallel('parallel run ecosys prep',
                                              task_array, input_array)

        # Execute the run ecosys tasks
        print("Prepared ecosys ... now preparing tasks for execution")
        exe_tasks = tigres.TaskArray('run_ecosys', [])
        for exe in run_ecosys_commands:
            exe_tasks.append(tigres.Task('run_ecosys_task', tigres.EXECUTABLE,
                                         impl_name=exe, input_types=[]))

        tigres.parallel('parallel run ecosys', exe_tasks,
                        tigres.InputArray('NULL', [tigres.InputValues()]))

        # Check the logs for the decode event
        # log_records = tigres.query(spec=["event = ecosys_task"])
        # for record in log_records:
        #     print(".. decoded {}".format(record.message))

        # Create DOT file (plain text graph)
        # tigres.dot_execution()

        # ## Check unfinished simulations and rerun from the last successful year
        # while not self.check_simulation():
        #     # Execute the run ecosys tasks
        #     print("Prepared ecosys ... now preparing tasks for execution")
        #     # Get run commands
        #     run_ecosys_commands = []
        #     for directory in self.dir_list:
        #         executable_sim = self.execlist['sim_name'] + ' < ' + \
        #                          self.execlist['ecosys_file'] + ' ' + '>logba'
        #         print(executable_sim)
        #         run_ecosys_commands.append("cd {}; {} 2>&1 ".format(directory, executable_sim))
        #         print(run_ecosys_commands)
        #
        #     print(run_ecosys_commands)
        #
        #     exe_tasks = tigres.TaskArray('run_ecosys', [])
        #     for exe in run_ecosys_commands:
        #         exe_tasks.append(tigres.Task('run_ecosys_task', tigres.EXECUTABLE,
        #                                      impl_name=exe, input_types=[]))
        #
        #     tigres.parallel('parallel run ecosys', exe_tasks,
        #                     tigres.InputArray('NULL', [tigres.InputValues()]))
        #
        #     # Check the logs for the decode event
        #     log_records = tigres.query(spec=["event = ecosys_task"])
        #     for record in log_records:
        #         print(".. decoded {}".format(record.message))
        #
        #     # Create DOT file (plain text graph)
        #     # tigres.dot_execution()

    def run_serial(self):
        """Run simulations in serial without tigres. The function does:

        1) run run_ecosys_DukeForest_prep
        2) run all simulations one by one in all folders

        :param: None
        :returns: None
        :rtype: None

        """
        for i in range(0, len(self.dir_list)):
            command = run_ecosys_DukeForest_prep(self.execlist, self.dir_list[i])
            subprocess.check_output(command, shell=True)

    def check_simulation(self):
        """Check which folders have unfinished simulations, and chane
        self.dir to these folders. Not needed after the bug was found

        :returns: self.dir_list
        :rtype: True/False. If all successful return True, else return False.

        """
        ## All the folders that have failed
        failure_cases = []
        ## The years when the simulation failed for all the failure folders
        failure_years = []

        ## Check which folders failed
        for folder in self.dir_list:
            res_file = glob.glob(folder + '/' + '01010' + self.year[-1] + 'hc')
            if not res_file:
                print('Simulation is not finished for folder ' + folder)
                failure_cases.append(folder)

        self.dir_list = failure_cases

        ## Check the years that failed
        for case in failure_cases:
            year = self.year[-1]
            file = glob.glob(case + '/' + '01010' + year + 'hc')
            while not file:
                year = str(int(year)-1)
                file = glob.glob(case + '/' + '01010' + year + 'hc')

            year = str(int(year)+1)
            failure_years.append(year)

        ## Change the opt files and the ecosys run file in each folder
        for case in failure_cases:
            idx = failure_cases.index(case)
            # change opt file
            opt_file = case + '/opt' + str(failure_years[idx])
            with open(opt_file, 'r') as f:
                content = f.readlines()
            if content[5] == 'NO\n':
                content[5] = 'YES\n'
            with open(opt_file, 'w') as f:
                f.write(''.join(content))

            # change ecosys run file
            run_file = case+'/'+self.execlist['ecosys_file']
            with open(run_file) as f:
                content = f.readlines()
            content[3] = str(int(self.year[-1])-int(failure_years[idx])+1) + ' 1\n'
            year_start = content[6][-5:-1]
            # idx_to_delete = np.arange(4, 4+15*(int(failure_years[idx])-int(year_start)))
            content = content[:4] + content[4+15*(int(failure_years[idx])-int(year_start)):]
            with open(run_file, 'w') as f:
                f.write(''.join(content))

        if failure_cases:
            return False
        else:
            return True

def run_ecosys_DukeForest_prep(execlist, dir_in):
    """prepare for running ecosys. The function does the following things:

    1) copy ecosys file (e.g., runba) to dir_in
    2) copy data file (none for ecosys) to dir_in
    3) return the command cd dir_in; call executable as a string

    :param execlist: dict, the execution list
    :param dir_in: string, the simulation directory
    :returns: the command to call executables in the form of a string
    :rtype: a string

    """
    try:
        shutil.copy("{sim_folder}/{ecosys_file}".format(**execlist),
                    dir_in)
    except IOError as e:
        raise Exception("File copy failed ({}): {}".format(
            execlist['ecosys_file'], e))

    for j in range(len(execlist['data_file'])):
        try:
            shutil.copy("{}/{}".format(execlist['sim_folder'],
                                       execlist['data_file'][j]), dir_in)
        except IOError as e:
            raise Exception("File copy failed ({}): {}".format(
                execlist['data_file'][j], e))

    ecosys_file = "{}/{}".format(execlist['sim_folder'],execlist['ecosys_file'])
    # cur_dir = os.path.dirname(dir_in)
    executable_sim = execlist['sim_name'] +  ' < ' + ecosys_file + ' ' + '>logba'
    return "cd {}; {} 2>&1 ".format(dir_in, executable_sim)
