#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'sjagsch'

"""
'master.py'
usage: python2.7 master.py [-p ~/path/to/FILE]

Master doc for Two-Photon (TP) calculations;

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

for parallelized version:   'master_parallel.py'

todo: long list...different polarizations

"""
# terminal input handling
import argparse
# date & time
from time import strftime, localtime, clock
# filepath handling
from inspect import getfile, currentframe
import os
import glob
# used to exit program in exception
from sys import exit
# mathematics
import numpy as np
# variables container
import variables as v
import functions as func


##################################################
masterpath = [os.path.dirname(os.path.abspath(getfile(currentframe()))), getfile(currentframe())]


def exists_and_file(inpath):
    """ checks whether input path leads to a file """
    if os.path.isfile(os.path.expanduser(inpath)):
        return os.path.dirname(os.path.expanduser(inpath))
    else:
        print 'Error when trying to locate file. Check path or working directory.'
        exit(1)


def dipole_x_polarization():
    """
    projection of the dipole elements on the field polarization unit-vector
    """
    pass

##################################################
# setting up argument parser for terminal input

parser = argparse.ArgumentParser()
parser.add_argument('-p', help='Set file path to data if not working directory.',
                    dest='inpath', type=exists_and_file,
                    default=glob.glob(masterpath[0] + '/*_cc.dat')[0],
                    metavar='~/path/to/FILE')
args = parser.parse_args()

##################################################
# write intro
print
print "Called:", '/'.join(masterpath)
tic = clock()
print strftime("on %B-%d-%Y at %H:%M:%S.", localtime()), '\nLet\'s calculate...'

##################################################
# <data is transposed> when read with unpack=True
# The infile has 11 columns where the last one can be dumped, thus range(0,10).
# format i(ndex)e, ih, E(ie), E(ih), (Re & Im of 3 dipole-components x,y,z), dump
try:
    # bem: needs check if the read file is not a broken simlink (~)??
    datavv = np.genfromtxt(glob.glob(args.inpath + '/*_vv.dat')[0], comments='#', skip_header=3,
                           unpack=False, usecols=range(0, 10),
                           autostrip=True)
    datacv = np.genfromtxt(glob.glob(args.inpath + '/*_cv.dat')[0], comments='#', skip_header=3,
                           unpack=False, usecols=range(0, 10),
                           autostrip=True)
    datacc = np.genfromtxt(glob.glob(args.inpath + '/*_cc.dat')[0], comments='#', skip_header=3,
                           unpack=False, usecols=range(0, 10),
                           autostrip=True)
except IndexError:
    print 'At least one of the files containing the dipole elements is missing.'
    exit('Abort...')
except ValueError:
    print 'Unexpected number of columns in at least one data file.'
    exit('Abort...')

print len(datacc[:, 0]), 'transitions x', len(datacc[0, :]), 'values cc-data-block.'
print len(datacv[:, 0]), 'transitions x', len(datacv[0, :]), 'values cv-data-block.'
print len(datavv[:, 0]), 'transitions x', len(datavv[0, :]), 'values vv-data-block.'

# set up arrays to store data relevant for matrix element calculation
evalcc = np.zeros((len(datacc[:, 0]), 7))
evalcv = np.zeros((len(datacv[:, 0]), 7))
evalvv = np.zeros((len(datavv[:, 0]), 7))
# total number of possible intermediate states
# TODO (andrei) what about the spin?
tot_intermed_states = len(datacc[:, 0]) + len(datacv[:, 0]) + len(datavv[:, 0]) - 2

# TODO set up as a function
for i in xrange(len(datacc[:, 0])):
    dipole = np.complex128(datacc[i, 4:10:2] + 1j * datacc[i, 5:10:2])
    dipxpolar1 = np.dot(v.E1polar, dipole)  # scalar product: (polarization unit vector, complex dipole vector)
    dipxpolar2 = np.dot(v.E2polar, dipole)

    evalcc[i, 0], evalcc[i, 1] = datacc[i, 0], datacc[i, 1]         # indices
    evalcc[i, 2] = datacc[i, 3] - datacc[i, 2]                      # delta-energy E_j-E_i
    evalcc[i, 3], evalcc[i, 4] = dipxpolar1.real, dipxpolar1.imag   # dipole projected on field 1
    evalcc[i, 5], evalcc[i, 6] = dipxpolar2.real, dipxpolar2.imag   # dipole projected on field 2

for i in xrange(len(datacv[:, 0])):
    dipole = np.complex128(datacv[i, 4:10:2] + 1j * datacv[i, 5:10:2])
    dipxpolar1 = np.dot(v.E1polar, dipole)
    dipxpolar2 = np.dot(v.E2polar, dipole)

    evalcv[i, 0], evalcv[i, 1] = datacv[i, 0], datacv[i, 1]
    evalcv[i, 2] = datacv[i, 3] - datacv[i, 2]
    evalcv[i, 3], evalcv[i, 4] = dipxpolar1.real, dipxpolar1.imag
    evalcv[i, 5], evalcv[i, 6] = dipxpolar2.real, dipxpolar2.imag

for i in xrange(len(datavv[:, 0])):
    dipole = np.complex128(datavv[i, 4:10:2] + 1j * datavv[i, 5:10:2])
    dipxpolar1 = np.dot(v.E1polar, dipole)
    dipxpolar2 = np.dot(v.E2polar, dipole)

    evalvv[i, 0], evalvv[i, 1] = datavv[i, 0], datavv[i, 1]
    evalvv[i, 2] = datavv[i, 3] - datavv[i, 2]
    evalvv[i, 3], evalvv[i, 4] = dipxpolar1.real, dipxpolar1.imag
    evalvv[i, 5], evalvv[i, 6] = dipxpolar2.real, dipxpolar2.imag

del datacc, datacv, datavv

##################################################
# matrix element calculations - mind sign deltaE for emis/abs
M_absorb = np.complex128(0.0 + 1j * 0.0)
M_emiss = np.complex128(0.0 + 1j * 0.0)

# TODO: with |i/f>_emis, |i/f>_abs set: implement:
# abs: |i> == |v_i>
# 1st step: transits of type: |v_i> --> |v> != |v_i> --> |c_f>
# 2nd step: transits of type: |v_i> --> |c> != |c_f> --> |c_f>
# emis: |i> == |c_i>
# 1st step: transits of type: |c_i> --> |c> != |c_i> --> |v_f>
# 2nd step: transits of type: |c_i> --> |v> != |v_f> --> |v_f>

# absorption - ggf. complex conjugate bilden!!?
for i in xrange(len(evalvv[:, 0])):
    if evalvv[i, 1] == v.absinit:
        continue
    else:  # index v.absfinal and from which |v> have to fit (--> how are files structured? --> evalcv's 1st index!!)
        # first step with field 1
        M_absorb += ((evalvv[i, 3] + 1j * evalvv[i, 4]) *
                     (evalcv[(v.absfinal - 1) * 4 + i % 4, 5] + 1j * evalcv[(v.absfinal - 1) * 4 + i % 4, 6]) /
                     (evalvv[i, 2] - v.hbar * v.w1))
        # first step with field 2
        M_absorb += ((evalvv[i, 5] + 1j * evalvv[i, 6]) *
                     (evalcv[(v.absfinal - 1) * 4 + i % 4, 3] + 1j * evalcv[(v.absfinal - 1) * 4 + i % 4, 4]) /
                     (evalvv[i, 2] - v.hbar * v.w2))

for i in xrange(len(evalvv[:, 0])):
    if evalcv[i, 1] == v.absfinal:
        continue

print 'Matrix-element squared is:', (M_absorb*np.conj(M_absorb)).real
print 'Nonsense of course, since field not chosen properly...etc.'

##################################################
# scaling factors

# rates TPA
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
print 'Scaling factors of the respective rates...'
print 'TPA bulk:', g0_tpa
print 'TPA double mode cavity:', gcav_tpa
print 'TPSE bulk:', g0_tpse
print 'TPSE double mode cavity:', gcav_tpse
print 'TPSIE bulk:', g0_tpsie
print 'TPSIE double mode cavity:', gcav_tpsie

##################################################
# write results to file
np.savetxt('results.txt', evalcc, delimiter='',
           header='Two-Photon Transition rates:\n' +
                  strftime("Timestamp: %B-%d-%Y, %H:%M:%S\n", localtime()) +
                  'format: TPA(b/c), TPSE(b/c), TPSIE(b/c)',
           comments='#')

##################################################
# write outro
print strftime("Finished %B-%d-%Y at %H:%M:%S", localtime())
print "Time elapsed:", clock() - tic, "s"
print "Have a nice day."
print

# EOF
##################################################