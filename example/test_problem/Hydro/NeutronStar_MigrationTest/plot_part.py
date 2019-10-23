# Plot the 1D denisty profile and generate movies

import re, os
import numpy as np
from glob import glob
from matplotlib import pyplot as plt

#plt.ion()
plt.ioff()

### load required parameters
parfile = "Input__Parameter"
parinfo = open(parfile).read()

output_part_y = re.findall(r"OUTPUT_PART_Y\s*(-?\d+\.?\d*[eE]?\d*)", parinfo)[0]
output_part_z = re.findall(r"OUTPUT_PART_Z\s*(-?\d+\.?\d*[eE]?\d*)", parinfo)[0]
output_part_y = float(output_part_y)
output_part_z = float(output_part_z)


unit_l = re.findall(r"UNIT_L\s*(-?\d+\.?\d*[eE]?\d*)", parinfo)[0]
unit_m = re.findall(r"UNIT_M\s*(-?\d+\.?\d*[eE]?\d*)", parinfo)[0]
unit_l = float(unit_l)
unit_m = float(unit_m)

file_fmt = "Xline_y{:.3f}_z{:.3f}_".format(output_part_y, output_part_z) + "{:06d}"


### plot function
def plotdens(nout):
    fn = file_fmt.format(nout)

    x, rho = np.genfromtxt(fn, usecols = [3, 6], unpack = 1)

    x   *= unit_l / 1.e5
    rho *= unit_m / unit_l**3

    plt.figure()
    plt.scatter(x, rho)

    plt.xlabel("x (km)")
    plt.ylabel(r"Density (g/cm$^3$)")
    plt.title("Nout = {}".format(nout))

    plt.tight_layout()
    plt.savefig(fn + ".png")

    plt.close()


fn_list = glob("Xline*")
fn_num = len(fn_list)

for i in range(fn_num):
    plotdens(i)


### generate animation
fn_in = "Xline_y{:.3f}_z{:.3f}_".format(output_part_y, output_part_z) \
      + r"%6d.png"
fn_out = "line_profile.mp4"
cmd = "ffmpeg -r 25 -i {} -vcodec h264 -pix_fmt yuv420p {}".format(fn_in, fn_out)

os.system(cmd)

