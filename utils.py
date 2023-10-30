import h5py as h5
import unyt
from scipy import spatial
import numpy as np
import yt
from yt.funcs import get_pbar
def snapshot(z, data_path, catalogue, res):
    start = 0
    print('...locating redshift..')
    if res == "3600":
        start = 78
    elif (res == "1800") or (res == "0900"):
        start = 77

    for i in range(start,0,-1):
        if (res == "3600") and i == 58:
            continue
        snapshot_ID = i
        halo_props_name = "/halo_properties_" + str(snapshot_ID).zfill(4) + ".hdf5"
        catalogue_path = data_path + catalogue + halo_props_name
        h5file = h5.File(catalogue_path, 'r')
        # Read in halo properties                                                                                                                                                                           
        groupName = "/SWIFT/Cosmology/"
        h5group = h5file[groupName]
        attrName =  'Redshift'
        redshift = h5group.attrs[attrName]
        h5file.close()
        print(redshift)
        if res == "3600":
            #higher resolution has z = 0.95 and 1.05 etc...                                                                                                                                                 
            if round(float(redshift),1) == z:
                return i
        elif (res == "1800") or (res == "0900"):
            if round(float(redshift),2) == z:
                return i
        else:
            continue
    raise Exception('Redshift not found...')