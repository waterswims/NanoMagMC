from IPython import get_ipython
ipython = get_ipython()
import mayavi.mlab as mlab
ipython.magic("gui qt")
# mlab.options.offscreen = True
import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
import subprocess
import os
from PyPDF2 import PdfFileMerger
from tqdm import tqdm


# File Name Stuff

file_name = sys.argv[1]
prefix = file_name.split(".h5")[0]
outname = prefix + ".pdf"

if len(file_name.split(".h5")) < 2:
    print("Please use a h5 file")
    sys.exit(-1)

# Read Data
print("Creating report from ", file_name)

f = h5py.File(file_name)
dataset = f["/Ts"]
Ts = dataset.value
dataset = f["/Hs"]
Hs = dataset.value
dataset = f["/Fulldat/mag_xs"]
mag_xs_all = dataset.value
mag_xs = np.mean(mag_xs_all, axis=2)
Nsamps = mag_xs_all.shape[2]
dataset = f["/Fulldat/mag_ys"]
mag_ys_all = dataset.value
mag_ys = np.mean(mag_ys_all, axis=2)
dataset = f["/Fulldat/mag_zs"]
mag_zs_all = dataset.value
mag_zs = np.mean(mag_zs_all, axis=2)
big_Ts = np.zeros_like(mag_zs)
for i in range(big_Ts.shape[0]):
    big_Ts[i, :] = Ts
chi_zs = (np.mean(mag_zs_all**2, axis=2) - mag_zs**2) / big_Ts
dataset = f["/Fulldat/mags"]
mags_all = dataset.value
mags = np.mean(mags_all, axis=2)
dataset = f["/Fulldat/energies"]
energies_all = dataset.value
energies = np.mean(energies_all, axis=2)
HCs = (np.mean(energies_all**2, axis=2) - energies**2) / (big_Ts**2)
dataset = f["/Fulldat/top_chars"]
top_chars_all = dataset.value
l_size = top_chars_all.shape[2]
top_chars = np.mean(top_chars_all, axis=3)

dataset = f["/Sing_Latt"]
pad = 0
if len(file_name.split("per_1")) < 2:
    if file_name.count("ha_s") > 0:
        pad = 2
    else:
        pad = 1
all_latt = dataset.value[:, :, pad:(l_size+pad), pad:(l_size+pad),
                         pad:(l_size+pad), :]
# all_latt /= Nsamps
f.close()
print("Reading Complete")

numTs = len(Ts)
numHs = len(Hs)

# Phase Diagrams

plt.figure()
plt.pcolor(mags[::-1, ::-1], cmap="plasma")
plt.colorbar()
plt.xlim(0, numTs)
plt.ylim(0, numHs)
plt.xticks(range(0, numTs, 10), Ts[::-10])
plt.yticks(range(0, numHs, 10), Hs[::-10])
plt.xlabel("T")
plt.ylabel("H")
plt.title("Absolute Magnetisation")
plt.savefig("mag.pdf")

plt.figure()
plt.pcolor(mag_zs[::-1, ::-1], cmap="plasma")
plt.colorbar()
plt.xlim(0, numTs)
plt.ylim(0, numHs)
plt.xticks(range(0, numTs, 10), Ts[::-10])
plt.yticks(range(0, numHs, 10), Hs[::-10])
plt.xlabel("T")
plt.ylabel("H")
plt.title("Magnetisation in Z-axis")
plt.savefig("magz.pdf")

plt.figure()
plt.pcolor(mag_xs[::-1, ::-1], cmap="plasma")
plt.colorbar()
plt.xlim(0, numTs)
plt.ylim(0, numHs)
plt.xticks(range(0, numTs, 10), Ts[::-10])
plt.yticks(range(0, numHs, 10), Hs[::-10])
plt.xlabel("T")
plt.ylabel("H")
plt.title("Magnetisation in X-axis")
plt.savefig("magx.pdf")

plt.figure()
plt.pcolor(mag_ys[::-1, ::-1], cmap="plasma")
plt.colorbar()
plt.xlim(0, numTs)
plt.ylim(0, numHs)
plt.xticks(range(0, numTs, 10), Ts[::-10])
plt.yticks(range(0, numHs, 10), Hs[::-10])
plt.xlabel("T")
plt.ylabel("H")
plt.title("Magnetisation in Y-axis")
plt.savefig("magy.pdf")

plt.figure()
plt.pcolor(chi_zs[::-1, ::-1], cmap="plasma")
plt.colorbar()
plt.xlim(0, numTs)
plt.ylim(0, numHs)
plt.xticks(range(0, numTs, 10), Ts[::-10])
plt.yticks(range(0, numHs, 10), Hs[::-10])
plt.xlabel("T")
plt.ylabel("H")
plt.title("Magnetic Susecptibility in Z-axis")
plt.savefig("chiz.pdf")

plt.figure()
chi_zs[chi_zs <= 0] = (chi_zs[chi_zs > 0]).min()
plt.pcolor(np.log(chi_zs[::-1, ::-1]), cmap="plasma")
plt.colorbar()
plt.xlim(0, numTs)
plt.ylim(0, numHs)
plt.xticks(range(0, numTs, 10), Ts[::-10])
plt.yticks(range(0, numHs, 10), Hs[::-10])
plt.xlabel("T")
plt.ylabel("H")
plt.title("log(Magnetic Susecptibility) in Z-axis")
plt.savefig("lnchiz.pdf")

plt.figure()
plt.pcolor(HCs[::-1, ::-1], cmap="plasma")
plt.colorbar()
plt.xlim(0, numTs)
plt.ylim(0, numHs)
plt.xticks(range(0, numTs, 10), Ts[::-10])
plt.yticks(range(0, numHs, 10), Hs[::-10])
plt.xlabel("T")
plt.ylabel("H")
plt.title("Heat Capacity")
plt.savefig("C.pdf")

plt.figure()
HCs[HCs <= 0] = (HCs[HCs > 0]).min()
plt.pcolor(np.log(HCs[::-1, ::-1]), cmap="plasma")
plt.colorbar()
plt.xlim(0, numTs)
plt.ylim(0, numHs)
plt.xticks(range(0, numTs, 10), Ts[::-10])
plt.yticks(range(0, numHs, 10), Hs[::-10])
plt.xlabel("T")
plt.ylabel("H")
plt.title("log(Heat Capacity)")
plt.savefig("lnC.pdf")

plt.figure()
plt.pcolor(energies[::-1, ::-1], cmap="plasma")
plt.colorbar()
plt.xlim(0, numTs)
plt.ylim(0, numHs)
plt.xticks(range(0, numTs, 10), Ts[::-10])
plt.yticks(range(0, numHs, 10), Hs[::-10])
plt.xlabel("T")
plt.ylabel("H")
plt.title("Internal Energy")
plt.savefig("ener.pdf")

plt.figure()
plt.pcolor(top_chars[::-1, ::-1, 0], cmap="plasma")
plt.colorbar()
plt.xlim(0, numTs)
plt.ylim(0, numHs)
plt.xticks(range(0, numTs, 10), Ts[::-10])
plt.yticks(range(0, numHs, 10), Hs[::-10])
plt.xlabel("T")
plt.ylabel("H")
plt.title("Topological charge of the first layer")
plt.savefig("TC0.pdf")

plt.figure()
plt.pcolor(top_chars[::-1, ::-1, top_chars.shape[2]//2], cmap="plasma")
plt.colorbar()
plt.xlim(0, numTs)
plt.ylim(0, numHs)
plt.xticks(range(0, numTs, 10), Ts[::-10])
plt.yticks(range(0, numHs, 10), Hs[::-10])
plt.xlabel("T")
plt.ylabel("H")
plt.title("Topological charge of the middle layer")
plt.savefig("TCmid.pdf")

plt.figure()
plt.pcolor(top_chars[::-1, ::-1, top_chars.shape[2]-1], cmap="plasma")
plt.colorbar()
plt.xlim(0, numTs)
plt.ylim(0, numHs)
plt.xticks(range(0, numTs, 10), Ts[::-10])
plt.yticks(range(0, numHs, 10), Hs[::-10])
plt.xlabel("T")
plt.ylabel("H")
plt.title("Topological charge of the lower layer")
plt.savefig("TClast.pdf")

# Plot Lattices
print("Phase diagrams complete")
images = ["mag.pdf", "magx.pdf", "magy.pdf", "magz.pdf", "chiz.pdf",
          "lnchiz.pdf", "ener.pdf", "C.pdf", "lnC.pdf", "TC0.pdf",
          "TCmid.pdf", "TClast.pdf"]

slice_show = [0, l_size//2, l_size-1]

fig = mlab.figure(size=(800, 596))

for H_ind in tqdm(range(numHs)):
    if H_ind%5 != 0:
        continue
    for T_ind in range(numTs):
        if T_ind%5 != 0:
            continue
        mlab.clf()

        ## Quiver of the slices
        Xs, Ys, Zs = np.mgrid[0:l_size, 0:l_size, 0:len(slice_show)]
        xspin3d = np.zeros_like(Xs, dtype=np.float)
        yspin3d = np.zeros_like(Xs, dtype=np.float)
        zspin3d = np.zeros_like(Xs, dtype=np.float)
        for k in range(len(slice_show)):
            Zs[:, :, k] = slice_show[k]
            xspin3d[:, :, k] = all_latt[H_ind, T_ind, :, :, slice_show[k], 0]
            yspin3d[:, :, k] = all_latt[H_ind, T_ind, :, :, slice_show[k], 1]
            zspin3d[:, :, k] = all_latt[H_ind, T_ind, :, :, slice_show[k], 2]
        pts = mlab.quiver3d(Xs, Ys, Zs, xspin3d, yspin3d, zspin3d, mode="cone",
                            resolution=50, scale_factor=2, scalars=zspin3d,
                            opacity=1)
        pts.glyph.color_mode = 'color_by_scalar'
        cbar = mlab.scalarbar(pts, title=r"Mz", orientation="vertical")
        pts.module_manager.scalar_lut_manager.lut_mode="viridis"
        pts.module_manager.scalar_lut_manager.label_text_property.color = (0.0,
            0.0, 0.0)
        pts.module_manager.scalar_lut_manager.title_text_property.color = (0.0,
            0.0, 0.0)

        ## Isosurfaces
        surf_skip = 0
        Xs, Ys, Zs = np.mgrid[surf_skip:l_size-surf_skip,
                              surf_skip:l_size-surf_skip,
                              0:l_size]
        full_xspin3d = all_latt[H_ind, T_ind, surf_skip:l_size-surf_skip,
                                surf_skip:l_size-surf_skip,
                                0:l_size, 0]
        full_yspin3d = all_latt[H_ind, T_ind, surf_skip:l_size-surf_skip,
                                surf_skip:l_size-surf_skip,
                                0:l_size, 1]
        full_zspin3d = all_latt[H_ind, T_ind, surf_skip:l_size-surf_skip,
                                surf_skip:l_size-surf_skip,
                                0:l_size, 2]
        cont = mlab.contour3d(Xs, Ys, Zs, full_zspin3d, contours=[0],
                              color=(0, 0.75, 0.75), opacity=0.5)

        ## Background Color
        fig = mlab.gcf()
        fig.scene.background = (1, 1, 1)

        mlab.title("T = {:.2f}, H = {:.2f}".format(Ts[T_ind], Hs[H_ind]), color=(0,0,0), size=0.7)

        name = "T_{}-H_{}.png".format(T_ind, H_ind)
        mlab.savefig(name)

        pdfname = "T_{}-H_{}.pdf".format(T_ind, H_ind)
        images.append(pdfname)
        sysout = subprocess.run(["convert", name, pdfname], stdout=subprocess.PIPE)
        # sysout = os.system("convert {} {}".format(name, pdfname))
        # print(sysout)
        os.remove(name)


# Merge PDFs and cleanup

merger = PdfFileMerger()

for pdf in images:
    merger.append(pdf)

merger.write(outname)

for pdf in images:
    os.remove(pdf)
