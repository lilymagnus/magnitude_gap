import yt
import unyt
import numpy as np
import h5py as h5

yt.enable_parallelism()

group_path = "/cosma8/data/dp004/flamingo/Runs/"
box = "L1000N1800"
run = "HYDRO_FIDUCIAL"
run_name = "HF"
data_path = group_path + box + "/" + run
catalogue = "/SOAP"
delta="500_crit" # density contrast [200_mean, 500_crit etc]   

if __name__ == "__main__":
    subh_dict_list = dict.fromkeys(np.array(subhalo_host_ids))
    
    storage{}
    i = 0
    
    bad_dict = {}
    for my_storage, key in yt.parallel_objects(sub_dict_list, storage=storage, dynamic=False):
        temp_dict = {}
