#! /usr/bin/env python
#
# -*- coding: utf-8 -*-
#
#  Purpose:
#    Plot the slice of density profile at the center (w/o grid)
#    and generate the movies
#
#  Last Updated: 2019/10/24
#  He-Feng Hsieh


import matplotlib
matplotlib.use("Agg")

import os
import yt
from glob import glob
import matplotlib.pyplot as plt


# setting
Plot_Density      = True  # plot the slice of density at the center
Plot_Density_grid = True  # plot the slice of density at the center with grid
Plot_Rhoc         = True  # plot the evolution of central density
Gene_animation    = True  # generate animation


# obtain all the HDF5 output
fn_in = glob("Data_" + "[0-9]" * 6)
fn_in.sort()

fn_out_fmt      = "{}_Slice_x_density.png"
fn_out_fmt_grid = "{}_Slice_x_density_grid.png"

rhoc_all = list()  # store the central density, in format of (time, rhoc)


# load data
for fn in fn_in:
    ds = yt.load(fn)

    # retrieve the simulation time of lv = 0
    time = ds.parameters["Time"][0]

    if Plot_Density:
        slc = yt.SlicePlot(ds, "x", "density", center = "c")
        slc.annotate_title("Time = {:.2f} ms".format(time))
        slc.save(fn_out_fmt.format(fn))

    if Plot_Density_grid:
        try:
            slc.annotate_grids()
        except:
            slc = yt.SlicePlot(ds, "x", "density", center = "c")
            slc.annotate_title("Time = {:.2f} ms".format(time))
            slc.annotate_grids()

        slc.save(fn_out_fmt_grid.format(fn))

    if Plot_Rhoc:
        center = ds.domain_center.tolist()
        rhoc = ds.r[center]["density"]

        rhoc_all.append((time, rhoc))

    # free memory
    del slc, ds


# plot the evolution of central density
if Plot_Rhoc:
    rhoc_all.sort(key = lambda x: x[0])  # sort by simulation time

    time, rhoc = zip(*rhoc_all)

    fig, ax = plt.subplots()

    ax.semilogy(time, rhoc, c = "k")
    ax.set_xlabel("Time (ms)")
    ax.set_ylabel(r"Central density (g/cm$^3$)")

    fig.tight_layout()
    plt.savefig("Rhoc_evolve.png")


# generate animation
if Gene_animation:
    if Plot_Density:
        fn_in  = r"Data_%6d_Slice_x_density.png"
        fn_out = "Slice_density.mp4"
        cmd = 'ffmpeg -r 25 -i {} -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -vcodec h264 -pix_fmt yuv420p {}'.format(fn_in, fn_out)

        os.system(cmd)


    if Plot_Density_grid:
        fn_in  = r"Data_%6d_Slice_x_density_grid.png"
        fn_out = "Slice_density_grid.mp4"
        cmd = 'ffmpeg -r 25 -i {} -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -vcodec h264 -pix_fmt yuv420p {}'.format(fn_in, fn_out)

        os.system(cmd)

