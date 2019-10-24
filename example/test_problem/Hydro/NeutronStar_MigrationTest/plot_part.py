#! /usr/bin/env python
#
# -*- coding: utf-8 -*-
#
#  Purpose:
#    Plot the 1D denisty profile and generate movies
#
#  Last Updated: 2019/10/24
#  He-Feng Hsieh


import matplotlib
matplotlib.use("Agg")

import re, os
import numpy as np
from glob import glob
from matplotlib import pyplot as plt


# setting
Gene_animation    = True  # generate animation


# get the unit given in Input__Parameter
parfile = "Input__Parameter"
parinfo = open(parfile).read()

unit_l = re.findall(r"UNIT_L\s*(-?\d+\.?\d*[eE]?\d*)", parinfo)[0]
unit_m = re.findall(r"UNIT_M\s*(-?\d+\.?\d*[eE]?\d*)", parinfo)[0]
unit_l = float(unit_l)
unit_m = float(unit_m)


# obtain all the Xline_* files
fn_in = glob("Xline_*" + "[0-9]" * 6)
fn_in.sort()

fn_out_fmt = "{}.png"


# plot
for idx, fn in enumerate(fn_in):
    x, rho = np.genfromtxt(fn, usecols = [3, 6], unpack = 1)

    x   *= unit_l / 1.e5       # km
    rho *= unit_m / unit_l**3  # cgs

    plt.figure()
    plt.scatter(x, rho)

    plt.xlabel("x (km)")
    plt.ylabel(r"Density (g/cm$^3$)")
    plt.title("Nout = {}".format(idx))

    plt.tight_layout()
    plt.savefig(fn_out_fmt.format(fn))

    plt.close()


# generate animation
if Gene_animation:
    # use the first output (idx = 0), and replace 000000 by %6d
    fn_in  = fn_in[0].replace("0" * 6, "%6d") + ".png"
    fn_out = "line_density.mp4"
    cmd = "ffmpeg -r 25 -i {} -vcodec h264 -pix_fmt yuv420p {}".format(fn_in, fn_out)

    os.system(cmd)

