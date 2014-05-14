#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'sjagsch'

"""
Module 'functions.py'

Contains:   functions used in master-script

Author: Stefan Jagsch
Mail:   stefanthomas.jagsch <at> gmail.com
Latest update: 2014-17-04

see also: 'master.py'

"""
import numpy as np

##################################################
# function definitions


def vec_angle_rad(x, y):
    """
    angle between two vectors in radians
    """
    scalar_prod = np.dot(x, y)
    abs_x = np.sqrt((x*x).sum())  # do not use linalg.norm.. bottleneck!
    abs_y = np.sqrt((y*y).sum())
    # angle in radians
    angle = np.arccos(scalar_prod / abs_x / abs_y)
    return angle


def vec_angle_deg(x, y):
    """
    angle between two vectors in degrees
    """
    scalar_prod = np.dot(x, y)
    abs_x = np.sqrt((x*x).sum())
    abs_y = np.sqrt((y*y).sum())
    # angle in radians
    angle = np.arccos(scalar_prod / abs_x / abs_y)
    # in degrees
    angle_deg = angle * 360 / 2 / np.pi
    return angle_deg

# EOF
##################################################