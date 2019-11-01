#! /usr/bin/env python
#
# -*- coding: utf-8 -*-
#
#  Purpose:
#    Plot the evolution of central density
#
#  Last Updated: 2019/11/01
#  He-Feng Hsieh


import matplotlib
matplotlib.use("Agg")

import numpy as np
import matplotlib.pyplot as plt


# load data
fn_in = "Record__CentralDens"

time, rhoc = np.genfromtxt(fn_in, usecols = [0, 2], unpack = 1)
rhoc /= 1e14

# plot the evolution of central density
fig, ax = plt.subplots()

ax.semilogy(time, rhoc, c = "k")
ax.set_xlabel("Time (ms)")
ax.set_ylabel(r"Central density ($10^{14}$ g/cm$^3$)")

fig.tight_layout()
plt.savefig("Rhoc_evolve.png")

