import yt 
import numpy as np
import cat_reader as cat
import h5py as h5
import unyt
import matplotlib.pyplot as plt
ds = yt.load("r_mag_gap/3600/m8_3600_z0_500c_1e13_50kpc.h5")
data = ds.data

group_path = "/cosma8/data/dp004/flamingo/Runs/"
box = "L1000N3600"
run = "HYDRO_FIDUCIAL"
catalogue = "/SOAP"
data_path = group_path + box + "/" + run
snapshot_ID = 78
halo_props_name = "/halo_properties_" + str(snapshot_ID).zfill(4) + ".hdf5"
delta_500="500_crit" # density contrast [200_mean, 500_crit etc]                                                                                                                  
delta_200 = '200_crit'
catalogue_path = data_path + catalogue + halo_props_name

def plot_offset():
    h5file = h5.File(catalogue_path, 'r')
    h5dset = h5file['/SO/200_crit/SORadius']
    R200c = h5dset[...] #Mpc
    
    h5dset = h5file['/VR/CentreOfPotential']
    CoP = h5dset[...] #Mpc                                                                                                                                                                               

    h5dset = h5file['VR/ID']
    VR_ID = h5dset[...]
    h5file.close()
    
    bgg_id = data['bgg_VR_id']
    host_id = data['host_id']
    M500c = data['M500c']

    offset = []
    mass = []
    #if bgg_id != host_id, there is an offset from center
    for idx, _id  in enumerate(bgg_id):
        if host_id[idx] == 1:
            continue
        if host_id[idx] != _id:
            mass.append(M500c[idx])
            
            #find location of bgg
            bgg_loc = CoP[np.where(VR_ID == _id)[0]]
            host_cop = CoP[np.where(VR_ID == host_id[idx])[0]] 
            
            r =np.sqrt(np.sum((bgg_loc-host_cop)**2))
            r = unyt.unyt_quantity(r,'Mpc')
            print(np.sqrt(r)/R200c[np.where(VR_ID == host_id[idx])[0]])
            offset.append(np.sqrt(r)/R200c[np.where(VR_ID == host_id[idx])[0]])
    

    plt.scatter(mass,offset)
    plt.xscale('log')
    plt.xlabel('Halo Mass [Msun]')
    plt.ylabel('CoP offset [r/R200c]')
    plt.savefig('offset.png')

def plot_4th_offset(cc):
    #ids for 4th BCG
    gal_ids = data['bgg4th_VR_id']
    clust_ids = data['host_id']
    mass = data['M500c']
    cop_gals  = cc.CoP[np.array(gal_ids-1, dtype=int)]
    cop_clust = cc.CoP[np.array(clust_ids-1, dtype=int)]

    dist_to_host_cop = cop_gals - cop_clust
    r_sqrd = np.sum(dist_to_host_cop**2, axis=1)
    plt.scatter(mass,np.sqrt(r_sqrd))
    plt.xscale('log')
    plt.xlabel('Halo Mass [Msun]')
    plt.ylabel('CoP offset [Mpc]')
    plt.savefig('offset.png')

if __name__ == "__main__":
    
    cc = cat.catalogue(catalogue_path, delta_500, delta_200)
    plot_4th_offset(cc)
