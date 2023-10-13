import numpy as np
import h5py as h5
import unyt
import matplotlib.pyplot as plt
import yt

if __name__ == "__main__":
     ds = yt.load("m8_1800_r_mag_gap_data_z0.h5")
     data = ds.data

     
     #changes
     plt.scatter(data['M500c'], data['mag_gap']
     breakpoint()
