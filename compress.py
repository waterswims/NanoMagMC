import numpy as np
import h5py
import sys

for file_name in sys.argv[1:]:
    prefix = file_name.split(".h5")[0]
    if len(file_name.split(".h5")) < 2:
        print("Please use a h5 file")
        sys.exit(-1)
    new_file_name = prefix + "(Compressed).h5"

    f = h5py.File(file_name, "r")
    f2 = h5py.File(new_file_name, "w")

    for name in ["Complete","Hs","Ts"]:
        full_name = name
        dataset = f[full_name]
        dat = dataset.value
        dset = f2.create_dataset(full_name, dat.shape, chunks=tuple(dat.shape),
                                 compression="gzip", shuffle=True,
                                 scaleoffset=3,
                                 compression_opts=9, data=dat)
        print(full_name, " complete")

    for name in ["Av_Latt", "Sing_Latt"]:
    # for name in ["Sing_Latt"]:
        test_dataset = f["/{}/T_0-H_0".format(name)]
        test_dat = test_dataset.value
        test_dataset = f["Fulldat/mags"]
        test_dat2 = test_dataset.value
        giant_latt = np.zeros([test_dat2.shape[0], test_dat2.shape[1],
                               test_dat.shape[0], test_dat.shape[1],
                               test_dat.shape[2], 3])
        for Hi in range(test_dat2.shape[0]):
            for Ti in range(test_dat2.shape[1]):
                full_name = "/{}/T_{}-H_{}".format(name, Ti, Hi)
                dataset = f[full_name]
                giant_latt[Hi, Ti, :, :, :, :] = dataset.value
        print(name, " reading complete")
        dset = f2.create_dataset(name, giant_latt.shape, shuffle=True,
                                 chunks=tuple(giant_latt.shape),
                                 scaleoffset=3,
                                 compression="gzip", compression_opts=9,
                                 data=giant_latt)
        print(name, " printing complete")

    for name in list(f['Fulldat']):
        if len(name.split("sub_mag")) > 1:
            continue
        full_name = "/Fulldat/{}".format(name)
        dataset = f[full_name]
        dat = dataset.value
        dset = f2.create_dataset(full_name, dat.shape, chunks=tuple(dat.shape),
                                 compression="gzip", shuffle=True,
                                 scaleoffset=3,
                                 compression_opts=9, data=dat)
        print(full_name, " complete")

    f2.close()
    f.close()
