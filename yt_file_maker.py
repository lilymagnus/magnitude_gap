'''
Script to store data in .h5 file using yt parallelisation
'''
import yaml
import yt
import unyt
import numpy as np
import h5py as h5
from collections import defaultdict
yt.enable_parallelism()
from yt.funcs import get_pbar

group_path = "/cosma8/data/dp004/flamingo/Runs/"
res = '0900'
box = "L1000N" + res
run = "HYDRO_FIDUCIAL"
run_name = "HF"
data_path = group_path + box + "/" + run
catalogue = "/SOAP"
delta_500="500_crit" # density contrast [200_mean, 500_crit etc]   
delta_200 = '200_crit'

def check_z(z):
    start = 0 
    print('...locating redshift..')
    if res == "3600":
        start = 78
    elif res == "1800":
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
        #print(redshift)
        if res == "3600":
            #higher resolution has z = 0.95 and 1.05 etc...
            if round(float(redshift),1) == z:
                return i
        elif res == "1800":
            if round(float(redshift),2) == z:
                return i
        else:
            continue
    raise Exception('Redshift not found...')

if __name__ == "__main__":
    
    z = 2
    #### read in data ####
    print('--------------------------------')
    print('--------------------------------')
    print('reading in data')
    snapshot_ID = check_z(z)
    print(snapshot_ID)
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
    h5dset = h5file['/SO/'+ delta_500 + '/TotalMass']
    M500c = h5dset[...]

    h5dset = h5file['/SO/'+ delta_200 + '/TotalMass']
    M200c = h5dset[...]
    
    h5dset = h5file['VR/ID']
    VR_ID = h5dset[...]
    
    h5dset = h5file['/SO/200_crit/SORadius']
    R200c = h5dset[...] #Mpc 

    h5dset = h5file['InclusiveSphere/50kpc/StellarMass']
    stellar_mass = h5dset[...]
    
    h5dset = h5file['InclusiveSphere/50kpc/StellarLuminosity']
    subhalo_luminosities = h5dset[...]
    
    #all subhalos have the same host id, halos that are hosts/empty = -1
    h5dset = h5file['VR/HostHaloID']
    host_id = h5dset[...]
    
    h5dset = h5file['InclusiveSphere/50kpc/TotalMass']
    aperture_mass = h5dset[...]
    
    h5dset = h5file['/VR/CentreOfPotential']
    CoP = h5dset[...]
    
    h5file.close()
    #get the r-band lums
    gal_lums = subhalo_luminosities[:,2]
    
    #mass cut the halos so we only have field halo ids, masses and lums
    VR_ID_hosts = VR_ID[M500c > 1e13]
    M500c_hosts = M500c[M500c > 1e13]
    gal_lums_hosts = gal_lums[M500c > 1e13]
    R200c_hosts = R200c[M500c > 1e13]
    CoP_hosts = CoP[M500c > 1e13]

    #remove host halos and get only subhalo host ids and lums
    host_id_cut = host_id[(host_id != -1)]
    gal_lums_cut = gal_lums[(host_id != -1)] 
    VR_ID_subs = VR_ID[(host_id != -1)]
    CoP_subs = CoP[(host_id != -1)]
    #get remove duplicates from the list 
    sub_dict_list = list(dict.fromkeys(np.array(host_id_cut)))
    pbar = yt.get_pbar('Fetching subhalo data', len(sub_dict_list))
    storage = {}
    
    for my_storage, i in yt.parallel_objects(range(len(sub_dict_list)), storage=storage, dynamic=False):
        temp_dict = {'M500c':[],'host_id':[],'lum_bgg':[],'lum_4thbgg':[], 'mag_gap':[], 'bgg_VR_id':[]}
        host_key = sub_dict_list[i]
        pbar.update(i)

        #if halo id not found, then halo did not pass the mass cut
        if len(np.where(VR_ID_hosts == host_key)[0]) == 0:
            my_storage.result_id = host_key
            my_storage.result = 0
            continue

        host_idx = np.where(VR_ID_hosts == host_key)[0]
        host_cop = CoP_hosts[host_idx]
        subh_lums = gal_lums_cut[(host_id_cut == host_key) & (gal_lums_cut > 0)]
        subh_coords = CoP_subs[(host_id_cut == host_key) & (gal_lums_cut > 0)]
                                     
        #only include subs that are within R200c
        dist_to_host_cop = subh_coords - host_cop
        radius = unyt.unyt_quantity(R200c_hosts[host_idx], 'Mpc')
        r_sqrd = np.sum(dist_to_host_cop**2, axis=1)
        subs_lums_R200 = subh_lums[np.where(r_sqrd < (radius)**2)[0]]
        #breakpoint()
        #ignore host halos that are too small (so have been removed from the VR list) 
        #ignore subhalo lists that are too small
        if (len(subs_lums_R200) < 5):
            my_storage.result_id = host_key
            my_storage.result = 0
        else:
            temp_dict['M500c'] = M500c_hosts[host_idx]            
            temp_dict['host_id'] = host_key
            
            #luminosity of bgg is just luminosity of host halo
            bgg_4th = np.sort(subs_lums_R200)[-3]
            bgg = gal_lums_hosts[host_idx]
            
            #if BCG < 4th-BCG
            if bgg < bgg_4th:
                bgg = np.sort(subs_lums_R200)[-1]
                bgg_4th = np.sort(subs_lums_R200)[-4]
                #if bgg < bgg4th then save the id of the brightest subhalo, else the 
                #brightest subhalo will correspond to the central 
                sub_ids = VR_ID_subs[(host_id_cut == host_key) & (gal_lums_cut > 0)]
                sub_ids_R200 = sub_ids[(np.where(r_sqrd < (radius)**2)[0])]
                bgg_id = sub_ids_R200[np.argsort(subs_lums_R200)][-1]
                temp_dict['bgg_VR_id'] = bgg_id
            else:
                temp_dict['bgg_VR_id'] = host_key
            
            temp_dict['lum_bgg'] = float(bgg)
            temp_dict['lum_4thbgg'] = float(bgg_4th)
            
            #bgg is more negative, so mag gap is (-bgg) + (4th bgg)
            mag_gap = -(-2.5 * np.log10(bgg)) + (-2.5* np.log10(bgg_4th))
            
            if (np.isnan(mag_gap)):
                my_storage.result_id = host_key
                my_storage.result = 0
                continue
            else:
                temp_dict['mag_gap'] = float(mag_gap)
                my_storage.result_id = host_key
                my_storage.result = temp_dict
            
            
    if yt.is_root():
        data = defaultdict(list)
        for key, value in storage.items():
            if value == 0:
                continue
            else:
                data['M500c'].append(value['M500c'])
                data['host_id'].append(value['host_id'])
                data['lum_bgg'].append(value['lum_bgg'])
                data['lum_4thbgg'].append(value['lum_4thbgg'])
                data['mag_gap'].append(value['mag_gap'])
                data['bgg_VR_id'].append(value['bgg_VR_id'])

        my_units = {'M500c':'Msun','host_id':'','lum_bgg':'Lsun','lum_4thbgg':'Lsun', 'mag_gap':'', 'bgg_VR_id':''}
        
        for field in data:
            data[field] = unyt.unyt_array(data[field],my_units[field])

        fake_ds = {'H0': H0, 'om_L':OmegaLambda, 'om_M':OmegaMatter}        
        yt.save_as_dataset(fake_ds, 'r_mag_gap/' + res + '/m8_' + res + '_z' + str(z) + '_500c_1e13_50kpc.h5', data= data) 
