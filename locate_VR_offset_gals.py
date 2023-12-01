'''
Script to check if offcenter CoM in weird gals
coincides with another entry in VR
'''
import numpy as np
import h5py as h5
import cat_reader as cat
import unyt
import matplotlib.pyplot as plt
import yt
from utils import snapshot, region, locate_gals,data_bin, paths
import unyt
from yt.funcs import get_pbar

group_path = "/cosma8/data/dp004/flamingo/Runs/"
run = "HYDRO_FIDUCIAL"
catalogue = "/SOAP"
box = 'L1000N3600'
data_path = group_path + box + "/" + run

if __name__ == "__main__":
    z = 0
    res = '3600'
    #SOAP cluster idx
    idx = 9037213
    snapshot_ID = snapshot(z, data_path, catalogue, res)
    halo_props_name = "/halo_properties_" + str(snapshot_ID).zfill(4) + ".hdf5"
    catalogue_path = data_path + catalogue + halo_props_name
    cc = cat.catalogue(catalogue_path,'50')
    
    cop_list = cc.CoP
    center = cc.CoM_200m[idx]
    
    indices = []
    pbar = yt.get_pbar('Comparing centers...', len(cop_list))
    #distance to corner of box
    for i,cop in enumerate(cop_list):
        pbar.update(i)
        if (np.abs(cop[0] - center[0]) < 0.1) & \
        (np.abs(cop[1] - center[1]) < 0.1) & \
        (np.abs(cop[2] - center[2]) < 0.1):
           indices.append(i) 

    breakpoint()
    
