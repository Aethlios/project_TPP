#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'sjagsch'

"""
'master_parallel.py'
usage: mpirun -n <#processes> python2.7 master_parallel.py [-p ~/path/to/FILE]

Master doc for Two-Photon (TP) calculations;
Using MPI-parallelization.

Abbreviations used:
        *0: bulk
        *cav: cavity
        tpa: TP absorption
        tpse: TP spontaneous emission
        tpsie: TP singly induced emission
        vv, cc, cv: intraband transitions in v(alence)/c(onduction) band;
                    interband transitions

Author: Stefan Jagsch
Mail:   stefanthomas.jagsch <at> gmail.com
Last updated:   April-16-2014

see also: 'variables.py'

"""
# terminal input handling
import argparse
# parallelization // mpi4py requires numpy objects to be passed
from mpi4py import MPI
# date & time
from time import strftime, localtime, clock
# filepath handling
from inspect import getfile, currentframe
from os import path
import glob
# used to exit program in exception
from sys import exit
# mathematics
import numpy as np
# variables/custom functions container
import variables as v
import functions as func


##################################################
masterpath = '/'.join([path.dirname(path.abspath(getfile(currentframe()))),
                      getfile(currentframe())])

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
comm.Barrier()


##################################################
# setting up argument parser for terminal input
def exists_and_file(inpath):
    """ checks whether input path is a file """
    if path.isfile(path.expanduser(inpath)):
        return path.dirname(path.expanduser(inpath))
    else:
        print 'Error when trying to locate file. Check path or working directory.'
        exit(1)

parser = argparse.ArgumentParser()
parser.add_argument('-p', help='Set file path to data if not working directory.',
                    dest='inpath', type=exists_and_file, default=masterpath,
                    metavar='~/path/to/FILE')
args = parser.parse_args()

##################################################
# write intro
if rank == 0:
    print
    print "Called:", masterpath
    tic = clock()
    print strftime("on %B-%d-%Y at %H:%M:%S.", localtime()), '\nLet\'s calculate...'

##################################################
# <data is transposed> when read with unpack=True
# infile has 11 columns where last one can be dumped
# format i(ndex)e, ih, E(ie), E(ih), (Re & Im of 3 dipole-components)
try:
    # bem: check if the read file is not a broken simlink (~)
    datavv = np.genfromtxt(glob.glob(args.inpath + '/*_vv.dat')[0], comments='#', skip_header=3,
                           unpack=True, usecols=range(0, 4),
                           autostrip=True)
    datacv = np.genfromtxt(glob.glob(args.inpath + '/*_cv.dat')[0], comments='#', skip_header=3,
                           unpack=True, usecols=range(0, 4),
                           autostrip=True)
    datacc = np.genfromtxt(glob.glob(args.inpath + '/*_cc.dat')[0], comments='#', skip_header=3,
                           unpack=True, usecols=range(0, 4),
                           autostrip=True)
except IndexError:
    print 'At least one of the files containing the dipole elements is missing.'
    exit('Abort...')
except ValueError:
    print 'Unexpected number of columns in at least one data file.'
    exit('Abort...')

for i in range(0, datacc[0].__len__()):
    pass

##################################################
# matrix element calculations
# before distribution; find zero-elements and dont distribute them
# else processes could be unbalanced.. (will be anyway?)

##################################################
# scaling factors
g0_tpa = v.pi / 2 * (v.P1 / (2 * v.hbar ** 2 * v.eps0 * v.c0 * v.A1 * v.n1)) * \
                    (v.P2 / (2 * v.hbar ** 2 * v.eps0 * v.c0 * v.A2 * v.n2))
gcav_tpa = v.pi / 2 * ((v.eta1 * v.P1 * v.Q1 * v.phi1) /
                       (2 * v.hbar ** 2 * v.eps0 * v.wc1 * v.Vmode1 * v.n1 ** 2)) * \
                      ((v.eta2 * v.P2 * v.Q2 * v.phi2) /
                       (2 * v.hbar ** 2 * v.eps0 * v.wc2 * v.Vmode2 * v.n2 ** 2))
# rates TPSE
g0_tpse = v.pi / 2 * ((v.w1 ** 3 * v.n1) / (3 * v.hbar * v.pi ** 2 * v.c0 ** 3 * v.eps0)) * \
                     ((v.w2 ** 3 * v.n2) / (3 * v.hbar * v.pi ** 2 * v.c0 ** 3 * v.eps0))
gcav_tpse = v.pi / 2 * ((2 * v.w1 * v.Q1 * v.phi1) /
                        (v.hbar * v.eps0 * v.n1 ** 2 * v.pi * v.wc1 * v.Vmode1)) * \
                       ((2 * v.w2 * v.Q2 * v.phi2) /
                        (v.hbar * v.eps0 * v.n2 ** 2 * v.pi * v.wc2 * v.Vmode2))
# rates TPSIE
g0_tpsie = v.pi / 2 * ((2 * v.w1 ** 3 * v.n1) /
                       (3 * v.hbar * v.pi ** 2 * v.c0 ** 3 * v.eps0)) * \
                      (v.P2 / (2 * v.hbar ** 2 * v.n2 * v.eps0 * v.c0 * v.A2))
gcav_tpsie = v.pi / 2 * ((4 * v.w1 * v.Q1 * v.phi1) /
                         (v.hbar * v.eps0 * v.n1 ** 2 * v.pi * v.wc1 * v.Vmode1)) * \
                        ((v.eta2 * v.P2 * v.Q2 * v.phi2) /
                         (2 * v.hbar ** 2 * v.eps0 * v.n2 ** 2 * v.w2 * v.Vmode2))
if rank == 0:
    print 'Scaling factors of the respective rates...'
    print 'TPA bulk:', g0_tpa
    print 'TPA double mode cavity:', gcav_tpa
    print 'TPSE bulk:', g0_tpse
    print 'TPSE double mode cavity:', gcav_tpse
    print 'TPSIE bulk:', g0_tpsie
    print 'TPSIE double mode cavity:', gcav_tpsie

##################################################
# collecting the data
#comm.Reduce(integral, total, op=MPI.SUM, root=0)
comm.Barrier()

##################################################
# write outro

if rank == 0:
    print strftime("Finished %B-%d-%Y at %H:%M:%S", localtime())
    print "Time elapsed:", clock() - tic, "s"
    print "Have a nice day."
    print

# EOF
##################################################