'''
Script to locate the galaxies with 
weird luminosities + stellar masses
and plot one of these metallicity vs age
'''
from matplotlib.ticker import FuncFormatter
import numpy as np
import h5py as h5
import cat_reader as cat
import unyt
import matplotlib.pyplot as plt
import yt
from yt.utilities.cosmology import Cosmology
from utils import snapshot, star_map, locate_gals,data_bin, paths
import unyt
from swiftsimio import mask
from swiftsimio import load

group_path = "/cosma8/data/dp004/flamingo/Runs/"
run = "HYDRO_FIDUCIAL"
catalogue = "/SOAP"
delta="500_crit" # density contrast [200_mean, 500_crit etc]                                                                                                                                            
delta_200 = '200_crit'

def _t_from_z(z, pos):

    co = Cosmology(omega_matter=0.304611,omega_lambda=0.693922,hubble_constant=0.681)
    val = co.t_from_z(z).to('Gyr')
    if val > 1:
        return "%d" % np.round(val)
    else:
        return "%.1f" % val

def _z_from_t(t, pos):
    co = Cosmology(omega_matter=0.304611,omega_lambda=0.693922,hubble_constant=0.681)
    return "%d" % np.round(co.z_from_t(co.quan(t, "Myr")))

def weird_gal_finder(z, cc, data):
    ids_gals = [] 
    lums = []
    st_mass = []
    #find index of clusters for the largest mass bin halos                                                                                                                                           \
    data_median, mass_list, bin_item = data_bin(np.log10(data['M500c']), data['host_id'], bin_num=20, stats=False)
    #gal_lums = cc.subhalo_luminosities[:,5]                                                                                                                                                          
    halo_id_list= bin_item[-1]

    for i in range(0,5):                                                                                                                                                                              
        halo_id = halo_id_list[i]                                                                                                                                                                     
        idx_cluster=int(halo_id-1)                                                                                                                                                                    
        
        ids_subs, gal_lums, pos, stellar_mass, tot_mass = locate_gals(cc, idx_cluster, z)                                                                                                            
        lums.extend(gal_lums[gal_lums > 0])                                                                                                                                                         
        ids_gals.extend(ids_subs[gal_lums > 0])
        st_mass.extend(stellar_mass[gal_lums > 0]) 
    
    breakpoint()
    return np.array(ids_gals)[np.argsort(lums)[0:9]], np.array(ids_gals)[np.argsort(lums)[20:29]]  


def mass_function(mass):
    bins = np.linspace(np.min(mass), np.max(mass), 20)
    hist, bin_edges = np.histogram(np.sort(mass), bins=bins)
    hist = hist[::-1].cumsum()[::-1]
    hist_update = []
    for ihist, hist_data in enumerate(hist):
        point = (hist_data / (bin_edges[ihist+1]-bin_edges[ihist]))
        hist_update.append(point)
 
    return hist_update, bin_edges

def plot_metals_from_swift():
    z = 0
    res = '1800'
    catalogue_path, snapshot_path, data = paths(res, z, group_path, catalogue, run)
    cc = cat.catalogue(catalogue_path, delta, delta_200, apert='50')
    idx_gals = wierd_gal_finder(z, cc, data)
    idx = idx_gals[0]
    
    xChoice=cc.CoP[idx,0]
    yChoice=cc.CoP[idx,1]
    zChoice=cc.CoP[idx,2]

    xCen = unyt.unyt_quantity(xChoice,'kpc')
    yCen = unyt.unyt_quantity(yChoice,'kpc')
    zCen = unyt.unyt_quantity(zChoice,'kpc')
    maxRegion = unyt.unyt_quantity(50,'kpc')
    maskRegion = mask(snapshot_path)

    #spatially mask the snapshot data around the cluster                                                                                                                                                 
    region=[[xCen-maxRegion,xCen+maxRegion],
            [yCen-maxRegion,yCen+maxRegion],
            [zCen-maxRegion,zCen+maxRegion]]

    maskRegion.constrain_spatial(region)
    data = load(snapshot_path,mask=maskRegion)

    metals = data.stars.metal_mass_fractions.value
    a = data.stars.birth_scale_factors
    z = (1/a) - 1

    #hist_update, bin_edges = mass_function(metals)
    median,_,_ = data_bin(z, metals, bin_num=20, stats=False)
    co = Cosmology(omega_matter=cc.OmegaMatter,omega_lambda=cc.OmegaLambda,hubble_constant=cc.h)
    zx = np.linspace(np.min(z), np.max(z), 20)
    fig, my_axes = plt.subplots()
    my_axes.plot(zx[:-1], median)
    xlim = ( zx[0], zx[-1])
    my_axes.set_xlim(zx[0], zx[-1])
    my_axes.xaxis.set_label_text("z")
    tx = my_axes.twiny()
    tx.xaxis.tick_top()
    my_times = co.arr([13, 5, 3, 1.5, 1, 0.5], 'Gyr')
    t_ticks_in_z = co.z_from_t(my_times)
    tx.xaxis.set_ticks(t_ticks_in_z)
    tx.xaxis.set_major_formatter(FuncFormatter(_t_from_z))
    tx.set_xlim(xlim)
    tx.xaxis.set_label_text("t [Gyr]")
    plt.savefig('metal.png')  
    breakpoint()

if __name__ == "__main__":
    z = 0
    res = '3600'
    catalogue_path, snapshot_path, data = paths(res, z, group_path, catalogue, run)
    cc = cat.catalogue(catalogue_path, delta, delta_200, apert='50')
    
    w_ids, n_ids = weird_gal_finder(z, cc, data)
    
    fig, my_axes = plt.subplots()
    for j, _id in enumerate(w_ids):
        t = cc.age[int(_id-1)]
        Z = cc.metallicity[int(_id-1)]
        my_axes.scatter(unyt.unyt_array(t,'Myr').to('Gyr'),Z,color='blue')
        #print(unyt.unyt_quantity(t,'Myr').to('Gyr'))
        #prinw_idt(Z) 

    my_axes.xaxis.set_label_text("t [Gyr]")
    my_axes.yaxis.set_label_text("Metallicity [Zsun]")
    plt.savefig('metals.png')
    breakpoint()
    
