import numpy as np
import h5py as h5
import unyt
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import yt
from utils import snapshot, DM_map, star_map, data_bin, paths
import unyt
from swiftsimio import mask
from swiftsimio import load
import cat_reader as cat
from swiftsimio.visualisation.smoothing_length_generation import generate_smoothing_lengths
from swiftsimio.visualisation.projection import scatter

group_path = "/cosma8/data/dp004/flamingo/Runs/"
run = "HYDRO_FIDUCIAL"
catalogue = "/SOAP"
size = 'L1000N'
res = '3600'

def plot_density_map(clust_idx,CoP, haloRadius, flamingoFile, z, center_list):
    dx, dy, dz, h, m, rChoice = DM_map(clust_idx,CoP, haloRadius, flamingoFile, z)
    maxRegion = unyt.unyt_quantity(0.5*rChoice,'Mpc')
    xChoice=CoP[clust_idx,0]
    yChoice=CoP[clust_idx,1]
    zChoice=CoP[clust_idx,2]

    mapRes = 600
    fig,axs = plt.subplots(figsize=(6,6))
    '''
    #for stars
    map=scatter(x=dx,y=dy,h=h,m=m,res=mapRes)
    map+=1.e-10
    map/=map.max()
    mapMin=-4
    mapMax=0
    '''
    #for DM
    map=scatter(x=dx,y=dy,h=h,m=m,res=mapRes)
    map+=1.e-10
    map/=map.max()
    mapMin=-3
    mapMax=0

    image=axs.pcolormesh(np.log10(map),cmap='magma',vmin=mapMin,vmax=mapMax)

    axs.axis('off')
    axs.set_box_aspect(1)
    axs.add_patch(Circle((0.5*mapRes,0.5*mapRes),rChoice*mapRes/(1.*maxRegion.value),fill=False,color='yellow',linewidth=2))
    
    if len(center_list) != 0:
        for cop in center_list:
            r = unyt.unyt_quantity(50,'kpc')
            x = cop[0] - xChoice
            y = cop[1] - yChoice
            xloc = (x +maxRegion.value)/(2.*maxRegion.value)
            yloc = (y +maxRegion.value)/(2.*maxRegion.value)
            axs.add_patch(Circle((yloc*mapRes,xloc*mapRes),r.to('Mpc')*mapRes/(2.*maxRegion.value),fill=False,color='red',linewidth=2))
                     
    fig.tight_layout()
    #axs.set_title(CoP[clust_idx,:])                                                                                                                                                                                                                                            
    plt.savefig('star_density.png',dpi=300)
    plt.close()
    breakpoint()

if __name__ == "__main__":
    z = 0
    catalogue_path, snapshot_path, ds = paths(res, z, group_path, catalogue, run, size, mass='200m',type_='mag')
    cc = cat.catalogue(catalogue_path, apert='50')
    data = ds.data

    clust_idx = 9037213 
    close_clust_VR = 8933463
    idxs = [9070882, 9071026, 9071121, 9071138, 9071596, 9071799, 9071922]
    
    print(cc.fof_subhalo_mass[9037213])
    plot_density_map(clust_idx,cc.CoM_200m, cc.R200m, snapshot_path, z, cc.CoP[idxs])
