import yt
import unyt
import numpy as np
import h5py as h5
from collections import defaultdict
yt.enable_parallelism()
from yt.funcs import get_pbar

group_path = "/cosma8/data/dp004/flamingo/Runs/"
box = "L1000N1800"
run = "HYDRO_FIDUCIAL"
run_name = "HF"
data_path = group_path + box + "/" + run
catalogue = "/SOAP"
delta="500_crit" # density contrast [200_mean, 500_crit etc]   

if __name__ == "__main__":
    #### read in data ####
    print('--------------------------------')
    print('reading in data')
    snapshot_ID = 77
    halo_props_name = "/halo_properties_" + str(snapshot_ID).zfill(4) + ".hdf5"

    catalogue_path = data_path + catalogue + halo_props_name
    h5file = h5.File(catalogue_path, 'r')

    # Read in halo properties                                                                                                                                                                             
    groupName = "/SWIFT/Cosmology/"
    h5group = h5file[groupName]
    attrName =  'H0 [internal units]'
    H0 = h5group.attrs[attrName]
    attrName="Omega_m"
    OmegaMatter=h5group.attrs[attrName]
    attrName="Omega_lambda"
    OmegaLambda=h5group.attrs[attrName]
    
    h5file = h5.File(catalogue_path, 'r')
    h5dset = h5file['/SO/'+ delta + '/TotalMass']
    M500c = h5dset[...]
    h5dset = h5file['VR/ID']
    VR_ID = h5dset[...]
    h5dset = h5file['InclusiveSphere/50kpc/StellarMass']
    stellar_mass = h5dset[...]
    h5dset = h5file['InclusiveSphere/50kpc/StellarLuminosity']
    subhalo_luminosities = h5dset[...]
    
    #all subhalos have the same host id, halos that are hosts/empty = -1
    h5dset = h5file['VR/HostHaloID']
    host_id = h5dset[...]
    h5dset = h5file['InclusiveSphere/50kpc/TotalMass']
    aperture_mass = h5dset[...]

    h5file.close()
    #get the r-band lums
    gal_lums = subhalo_luminosities[:,2]
    
    #remove low mass halos and get only subhalo host ids
    host_id_cut = host_id[(host_id != -1)]

    pbar = yt.get_pbar('Fetching subhalo data', len(host_id_cut))
    sub_dict_list = list(dict.fromkeys(np.array(host_id_cut)))
    storage = {}
    i = 0
    breakpoint()
    bad_dict = {'M500c':-1,'host_id':-2,'lum_bgg':-1,'lum_4thbgg':-1}
    for my_storage, key in yt.parallel_objects(sub_dict_list, storage=storage, dynamic=False):
        temp_dict = {'M500c':[],'host_id':[],'lum_bgg':[],'lum_4thbgg':[], 'mag_gap':[]}
        
        host_key = key
        pbar.update(i)
        
        host_idx = np.where(VR_ID == host_key)[0]
        subh_lums = gal_lums[host_id == host_key]

        if len(subh_lums) < 4:
            my_storage.result_id = host_key
            my_storage.resutl = bad_dict
        else:
            temp_dict['M500c'] = M500c[host_idx]
            temp_dict['host_id'] = host_key
            #luminosity of bgg is just luminosity of host halo
            temp_dict['lum_bgg'] = gal_lums[host_idx]
            temp_dict['lum_4thbgg'] = subh_lums[-4]
            temp_dict['mag_gap'] = np.abs((-2.5 * np.log10(gal_lums[host_idx])) - 
                                          (-2.5* np.log10(subh_lums[-4])))
            my_storage.result_id = key
            my_storage.result = temp_dict
            
            
    if yt.is_root():
        data = defaultdict(list)
        for key, value in storage.items():
            data['M500c'].append(value['M500c'])
            data['host_id'].append(value['host_id'])
            data['lum_bgg'].append(value['lum_bgg'])
            data['lum_4thbgg'].append(value['lum_4thbgg'])
            data['mag_gap'].append(value['mag_gap'])

            
        fake_ds = {'H0': H0, 'om_L':OmegaLambda, 'om_M':OmegaMatter}

        yt.save_as_dataset(fake_ds, '1Gpc_r_mag_gap_data_z0.h5', data= data) 
