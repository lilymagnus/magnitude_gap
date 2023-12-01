'''
For each cluster, take the galaxies within R200 (or whatever you use) and 
make a histogram in log luminosity bins. Then divide the y-axis by the 
volume of the cluster sphere to get a number density.
'''

import numpy as np
import h5py as h5
import cat_reader as cat
import unyt
import matplotlib.pyplot as plt
import yt
from utils import snapshot, star_map, locate_gals,data_bin, paths
import unyt
from velociraptor.tools.luminosity_functions import create_luminosity_function

group_path = "/cosma8/data/dp004/flamingo/Runs/"
run = "HYDRO_FIDUCIAL"
catalogue = "/SOAP"
delta="500_crit" # density contrast [200_mean, 500_crit etc]                                                                                                                                                                                                                  
delta_200 = '200_crit'
 
def mass_function(mass,boxsize,cumulative=False):     
    if np.min(mass) == 0:
        raise Exception('Remove galaxies with no stars!')
    mass= np.sort(mass)
    bins = np.logspace(np.log10(mass[0]), np.log10(mass[-1]), 15)                                                                                                                        
    hist, bin_edges = np.histogram(mass, bins=bins)  
    if cumulative == True:
        hist = hist[::-1].cumsum()[::-1]                                                                                                                                                              
    hist_update = []                                                                                                                                                                              
    for ihist, hist_data in enumerate(hist):                                                                                                                                                     
        point = hist_data/ ((np.log10(bin_edges[ihist+1])-np.log10(bin_edges[ihist])))                                                                                                           
        hist_update.append(point / boxsize)   
    return hist_update, bin_edges

def luminosity_function(lums, boxsize, cumulative=False):
    if np.min(lums) == 0:
        raise Exception('Remove galaxies with no stars!')
    lums = np.sort(lums)
    bins =  np.logspace(np.log10(lums[0]), np.log10(lums[-1]), 15)
    #normalization_factor = 1.0 / (bin_width * box_volume)
    hist, bin_edges = np.histogram(lums,bins)
    if cumulative == True:
        hist = hist[::-1].cumsum()[::-1]
    hist_update = []
    for ihist, hist_data in enumerate(hist):
        bin_width = (np.log10(bin_edges[ihist+1])-np.log10(bin_edges[ihist]))
        hist_update.append(hist_data / (bin_width * boxsize))
    
    return hist_update, bin_edges
  
if __name__ == "__main__":
    z = 0
    res_list = ['3600','1800']
    for res in res_list:
        catalogue_path, snapshot_path, data = paths(res, z, group_path, catalogue, run, type_='mag_500c')
        data_median, mass_list, bin_item = data_bin(np.log10(data['M500c']), data['host_id'], bin_num=20, stats=False)

        cc = cat.catalogue(catalogue_path, delta, delta_200, apert='50')
    
        #gal_lums = cc.subhalo_luminosities[:,5]
        halo_id_list= bin_item[-1]
        #subhalo mass function for largest cluster
        #cluster_idx = int(halo_id_list[np.argmax(cc.M500c[np.array(halo_id_list, dtype=int)-1])])-1
        #luminosity_function(cc, cluster_idx, res, z)
        #mass_function(mass,idx_cluster,cc,cumulative=False)
        
        
        LF = []
        MF = []
        lums = []
        st_mass = []
        #find index of clusters for the largest mass bin halos                                                                                                                                                                                                                     
        for i in range(0,5):
            halo_id = halo_id_list[i]
            print(halo_id)
            idx_cluster=int(halo_id-1)
            
            idx_gals, gal_lums, pos, stellar_mass, tot_mass = locate_gals(cc, idx_cluster, z)
            lums.extend(gal_lums[gal_lums > 0])
            st_mass.extend(stellar_mass[gal_lums > 0])
            
            #lum_func, bin_edges = luminosity_function(gal_lums[gal_lums > 0], ((4/3) * np.pi * cc.R200c[idx_cluster] ** 3), cumulative=False)
            #LF.append(lum_func)
            #mass_func, bin_edges = mass_function(tot_mass,  ((4/3) * np.pi * cc.R200c[idx_cluster] ** 3), cumulative=False)  
            #MF.append(mass_func) 

        #plt.loglog(bin_edges[:-1], np.mean(LF, axis=0),label=res)
        #plt.loglog(bin_edges[:-1], np.mean(MF, axis=0),label=res) 
        plt.scatter(st_mass,lums,s=10,label=res)
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'Stellar Mass [$M^*_{\odot}$]')
    plt.ylabel(r'Luminosity [$L_{\odot}$]')
    plt.legend()
    #plt.title('Subhalo Mass function') 
    #plt.title('Luminosity function')
    plt.xlabel(r'Stellar Mass [$M^*_{\odot}$]') 
    plt.ylabel(r'Luminosity [$L_{\odot}$]')
    #plt.ylabel(r'dn/dM $[Mpc^{-3}]$')
    plt.savefig('star_lum.png')
    breakpoint()
       
