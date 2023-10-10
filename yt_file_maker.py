import yt
import unyt
import numpy as np
import h5py as h5
from collections import defaultdict
yt.enable_parallelism()

group_path = "/cosma8/data/dp004/flamingo/Runs/"
box = "L1000N1800"
run = "HYDRO_FIDUCIAL"
run_name = "HF"
data_path = group_path + box + "/" + run
catalogue = "/SOAP"
delta="500_crit" # density contrast [200_mean, 500_crit etc]   

if __name__ == "__main__":
    

    subh_dict_list = list(dict.fromkeys(np.array(subhalo_host_ids)))
    storage{}
    i = 0
    
    bad_dict = {}
    for my_storage, i in yt.parallel_objects(sub_dict_list, storage=storage, dynamic=False):
        temp_dict = {}
        key = 
        pbar.update(key)
        if len(subhalos)< 4:
            my_storage.result_id = key
            my_storage.resutl = bad_dict
        else:
            my_storage.result_id = key
            my_storage.resutl = temp_dict
            
            
    if yt.is_roo():
        data = defaultdict(list)
        for key, value in storage.items():
            data[].append(value[])

            
        fake_ds = {'H0':, 'om_L':, 'om_M':}
        yt.save_as_dataset((fake_ds, 'name.h5', data=data) 
