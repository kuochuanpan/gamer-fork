# Plot the slice of density profile at the center (w/o grid)
# and generate the movie

import os
import yt
from glob import glob
import matplotlib.pyplot as plt


fn_list = glob("Data*")
rhoc_all = list()  # (time, rhoc)


for fn in fn_list:
    ds = yt.load(fn)

    # plot the slice (with grid)
    slc1 = yt.SlicePlot(ds, "x", "density", center = "c")
    slc1.annotate_grids()
    slc1.save()

    # rename the output file
    fn_out_old = fn + "_Slice_x_density.png"
    fn_out_new = fn + "_Slice_x_density_grid.png"
    os.system("mv {} {}".format(fn_out_old, fn_out_new))

    # plot the slice (without grid)
    slc2 = yt.SlicePlot(ds, "x", "density", center = "c")
    slc2.save()


    # collect the central density
    time = ds.parameters['Time'][0]
    center = ds.domain_center.tolist()
    rhoc = ds.r[center]["density"]

    rhoc_all.append((time, rhoc))

    # free memory
    del slc1, slc2, ds


# plot the evolution of central density
rhoc_all.sort(key = lambda x: x[0])  # sort by simulation time

time, rhoc = zip(*rhoc_all)

fig, ax = plt.subplots()

ax.semilogy(time, rhoc, c = "k")
ax.set_xlabel("Time (ms)")
ax.set_ylabel(r"Central density (g/cm$^3$)")

fig.tight_layout()
plt.savefig("Rhoc_evolve.png")


# generate animation
fn_in = r"Data_%6d_Slice_x_density.png"
fn_out = "Slice_profile.mp4"
cmd = 'ffmpeg -r 25 -i {} -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -vcodec h264 -pix_fmt yuv420p {}'.format(fn_in, fn_out)

os.system(cmd)


fn_in = r"Data_%6d_Slice_x_density_grid.png"
fn_out = "Slice_profile_grid.mp4"
cmd = 'ffmpeg -r 25 -i {} -vf "pad=ceil(iw/2)*2:ceil(ih/2)*2" -vcodec h264 -pix_fmt yuv420p {}'.format(fn_in, fn_out)

os.system(cmd)


