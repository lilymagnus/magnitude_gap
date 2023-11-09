import numpy as np
import h5py as h5
import unyt
import matplotlib.pyplot as plt
import yt
import scipy
from yt.funcs import get_pbar

group_path = "/cosma8/data/dp004/flamingo/Runs/"
run = "HYDRO_FIDUCIAL"
catalogue = "/SOAP"
delta="500_crit"

def mass_bin(st_mass, r_gap):
     #sort and bin data in terms of halo masses                                                                                                                                                             
     array = tuple(zip(st_mass, r_gap))
     #sorting the tuple in order of mass values                                                                                                                                                             
     sorted_list = sorted(array, key=lambda x:x[0])
     mass_sort, mag_sort = zip(*sorted_list)

     mass_list = list(np.log10(mass_sort))
     mag_list = list(mag_sort)

     bins = np.linspace(np.min(mass_list), np.max(mass_list), 20)

     n,_=np.histogram(mass_list,bins=bins)

     r_gap_median = []
     for i in range(len(n)):
          r_gap_section = list(np.array(mag_list[0:n[i]]).flatten())
          r_gap_median.append(np.median(r_gap_section))
          del mag_list[0:n[i]]

     return r_gap_median, mass_list, bins


def SM_r_mag_data(res):
    box = "L1000N" + res
    data_path = group_path + box + "/" + run

    ds = yt.load("r_mag_gap/" + res + '/m8_' + res +"_z0_500c_1e13_50kpc.h5")
    if res == '3600':
        snapshot_ID = 78
    else:
        snapshot_ID = 77

    halo_props_name = "/halo_properties_" + str(snapshot_ID).zfill(4) + ".hdf5"
    catalogue_path = data_path + catalogue + halo_props_name

    snapshot_name = "/flamingo_" + str(snapshot_ID).zfill(4)
    snapshot_path = data_path + "/snapshots" + snapshot_name + snapshot_name + ".hdf5"

    #use halo_ids to find index of cluster in full catalogue list = index in snapshot                                                                                                                       
    h5file = h5.File(catalogue_path, 'r')

    #full list of all halo ids                                                                                                                                                                             
    h5dset = h5file['VR/ID']
    VR_ID = h5dset[...]
    h5dset = h5file['InclusiveSphere/50kpc/StellarMass']
    gal_mass = h5dset[...]
    #take r-band only                                                                                                                                                                                    
    h5dset = h5file['InclusiveSphere/50kpc/StellarLuminosity']
    gal_lums = h5dset[...]

    gal_lums = gal_lums[:,2]

    h5file.close()
    
    
    data = ds.data
    halo_id = data['host_id']
    
    stellar_mass = []
    r_mag = []
    #pbar = yt.get_pbar('Fetching subhalo data', len(halo_id))
    for i, _id in enumerate(halo_id):
        if _id != data['bgg_VR_id'][i]:
             halo_idx = np.where(VR_ID == data['bgg_VR_id'][i])[0]
        else:
             halo_idx  = np.where(VR_ID == _id)[0]
              
        stellar_mass.append(gal_mass[halo_idx])
        r_mag.append(-2.5*np.log10(gal_lums[halo_idx]))

    print(res + ' done')
    return stellar_mass, r_mag

if __name__ == "__main__":
    res_list = ['1800','3600']
    
    for res in res_list:
        stellar_mass, r_mag = SM_r_mag_data(res)
        r_gap_median, mass_list, bins = mass_bin(stellar_mass, r_mag)
        x = np.linspace(np.min(mass_list), np.max(mass_list), len(bins)-1)
        plt.scatter(x, r_gap_median, label=res)
         
    plt.xscale('log')
    plt.xlabel('Stellar Mass [Msun]')
    plt.ylabel('r-magnitude')
    plt.legend()
    plt.savefig('res_test.png')
    #breakpoint()
    
