#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'sjagsch'

"""
Module 'variables.py'

Contains:  settings for Two-Photon calculations,
            Field parameters for the involved fields,
            physical constants needed as supplied by scipy
            with more convenient names.
           functions used in master-script


Author: Stefan Jagsch
Mail:   stefanthomas.jagsch <at> gmail.com
Latest update: 2014-14-04

see also: 'master.py'

"""
import scipy.constants
import numpy as np


##################################################
# settings
absinit = 2           # initial and final state for abs/emis process
absfinal = 2
emisinit = 2
emisfinal = 2

##################################################
# field 1 parameter
w1 = 1.4e15             # ang. frequency (rad/s)
wc1 = 1.44e15           # resonance frequency (rad/s)
P1 = 10.0e-6            # incident power (W)
# equiv. ~1kW/cm^2; how about that?
A1 = 1.0e-12            # area (m^2), intensity I = P/A
Q1 = 1000.0             # quality factor
Vmode1 = 1.3e-20        # mode volume (m^3)
n1 = 2.429              # index of refraction (wikipedia)
E1polar = [1.0, 0.0, 0.0]   # field polarization vector (x, y, z)
E1max = 15.0            # field maximum (needed?)
E1local = 11.0          # local field at emitter (needed?)
phi1 = 0.5              # frequency missmatch factor (if needed)
eta1 = 0.05             # incoupling efficiency

# field 2 parameter
w2 = 2.6e15
wc2 = 2.588e15
P2 = 1.0
A2 = 1.0
Q2 = 1000.0
Vmode2 = 1.3e-20
n2 = 2.429
E2polar = [1.0, 0.0, 0.0]
E2max = 15.0
E2local = 11.0
phi2 = 0.5
eta2 = 0.05

##################################################
# physical constants (CODATA 2010)
me0 = scipy.constants.electron_mass
hbar = scipy.constants.hbar
e0 = scipy.constants.elementary_charge
c0 = scipy.constants.speed_of_light
eps0 = scipy.constants.epsilon_0
pi = scipy.constants.pi

# EOF
##################################################