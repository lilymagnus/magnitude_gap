import numpy as np
import h5py as h5
import unyt
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import yt
from utils import snapshot, star_map, data_bin
import unyt
from swiftsimio import mask
from swiftsimio import load
import cat_reader as cat
from swiftsimio.visualisation.smoothing_length_generation import generate_smoothing_lengths
from swiftsimio.visualisation.projection import scatter
#choose halo with small and large magnitude gap

def halo_finder(gap, z, res):
    #find halos according to mag gap
    #look through presaved dataset 
    ds = yt.load("r_mag_gap/" + res +"/m8_" + res +"_z"+str(z)+"_500c_1e13_50kpc.h5")
    data = ds.data
    bgg_4th = -2.5 * np.log10(data['lum_4thbgg'])
    bgg = -2.5 * np.log10(data['lum_bgg'])
    halo_id = data['host_id']
    mass = data['M500c']
    print('finding halo id...')
    for i in range(len(halo_id)):
        mag_gaps =  np.abs(bgg[i] - bgg_4th[i])
        if np.abs(mag_gaps - gap) <= 0.05:
            print(bgg[i])
            print( bgg_4th[i])
            print(mag_gaps)
            return halo_id[i]

    raise Exception ('no magnitude gap near what has been input')

def plot_density_map(clust_idx, idx_gals,CoP, haloRadius, flamingoFile, z, mag_gals, gal_lums, gal_pos):    
    dx, dy, dz, h, m, rChoice = star_map(clust_idx,CoP, haloRadius, flamingoFile, z)
    maxRegion = unyt.unyt_quantity(1*rChoice,'Mpc')
    xChoice=CoP[clust_idx,0]
    yChoice=CoP[clust_idx,1]
    zChoice=CoP[clust_idx,2]

    mapRes = 600    
    fig,axs = plt.subplots(figsize=(6,6))
    
    map=scatter(x=dx,y=dy,h=h,m=m,res=mapRes) 
    map+=1.e-10
    map/=map.max()
    mapMin=-4
    mapMax=0
    
    image=axs.pcolormesh(np.log10(map),cmap='magma',vmin=mapMin,vmax=mapMax)
    
    axs.axis('off')
    axs.set_box_aspect(1) 
    axs.add_patch(Circle((0.5*mapRes,0.5*mapRes),rChoice*mapRes/(2.*maxRegion.value),fill=False,color='yellow',linewidth=2))
    #VR ID is just the index+1 since the VR_ID list starts at 1 
    cs = ['red', 'white','green','yellow']
    for i in range(4):
        rsub= unyt.unyt_quantity(50,'kpc')
        xsub= gal_pos[i][0] - xChoice
        ysub= gal_pos[i][1] - yChoice
        xsub = (xsub +maxRegion.value)/(2.*maxRegion.value)
        ysub = (ysub +maxRegion.value)/(2.*maxRegion.value)

        #scatter image is inverted so flip the subhalo x and y coords
        axs.add_patch(Circle(( ysub*mapRes,  xsub*mapRes), rsub.to('Mpc')*mapRes/(2.*maxRegion.value),fill=False,linewidth=2, color=cs[i]))
        axs.annotate('%.2f' % mag_gals[i], (((0.01 + ysub)*mapRes), ((+0.02 + xsub)*mapRes)), color=cs[i])
 
    bgg_subx, bgg_suby, bgg_rsub = 0.5, 0.5,  unyt.unyt_quantity(50,'kpc')
    axs.add_patch(Circle((bgg_subx*mapRes, bgg_suby*mapRes),bgg_rsub.to('Mpc')*mapRes/(2.*maxRegion.value),fill=False,color='blue',linewidth=2))
    axs.annotate('%.2f' % (-2.5*np.log10(gal_lums[clust_idx])), (((0.01 + bgg_suby)*mapRes), ((-0.02 + bgg_subx)*mapRes)),color='blue')
    axs.annotate(np.round(CoP[clust_idx,:],2), (((0.5)*mapRes), ((1)*mapRes)),color='blue')
    fig.tight_layout()  
    #axs.set_title(CoP[clust_idx,:])
    plt.savefig('star_density.png',dpi=300)
    plt.close()
    breakpoint()
    
group_path = "/cosma8/data/dp004/flamingo/Runs/"
run = "HYDRO_FIDUCIAL"
catalogue = "/SOAP"
delta="500_crit" # density contrast [200_mean, 500_crit etc]
delta_200 = '200_crit'

def paths(res, z):
    box = "L1000N" + res
    data_path = group_path + box + "/" + run
     
    snapshot_ID = snapshot(z, data_path, catalogue, res)
    halo_props_name = "/halo_properties_" + str(snapshot_ID).zfill(4) + ".hdf5"
    catalogue_path = data_path + catalogue + halo_props_name
    
    snapshot_name = "/flamingo_" + str(snapshot_ID).zfill(4)
    snapshot_path = data_path + "/snapshots" + snapshot_name + snapshot_name + ".hdf5"
    
    ds = yt.load("r_mag_gap/" + res +"/m8_" + res +"_z"+str(z)+"_500c_1e13_50kpc.h5")
    data = ds.data

    return catalogue_path, snapshot_path, data

def locate_gals(cc, idx_cluster, num):
    halo_id = cc.VR_ID[idx_cluster]

    print(cc.CoP[idx_cluster])

    #find index of subhalos brightest galaxy                                                                                                                                                                                                                                   
    sub_lums = gal_lums[cc.host_ID == halo_id]
    sub_VRs = cc.VR_ID[cc.host_ID == halo_id]
    sub_CoP = cc.CoP[cc.host_ID == halo_id]

    idx_gals = []
    mag_gals = []
    pos = []

    cluster_radius = cc.R200c[idx_cluster]
    cluster_CoP = cc.CoP[idx_cluster]
    dist_to_host_cop = cc.CoP  - cluster_CoP
    r_sqrd = np.sum(dist_to_host_cop**2, axis=1)

    subs_ids = cc.VR_ID[np.where(r_sqrd < cluster_radius**2)[0]]
    lums = gal_lums[np.where(r_sqrd < cluster_radius**2)[0]]
    sub_CoP = cc.CoP[np.where(r_sqrd < cluster_radius**2)[0]]
    
    for i in range(num):
        j = i + 1
        VR_gals = sub_VRs[np.argsort(lums)[-j]]
        idx_gals.append(int(VR_gals - 1))
        mag_gals.append(-2.5*np.log10(np.sort(lums)[-j]))
        pos.append(sub_CoP[np.argsort(lums)[-j]])
        print(sub_CoP[np.argsort(lums)[-j]])

    #breakpoint()
    return idx_gals, mag_gals, pos

if __name__ == "__main__":
    z = 0
    catalogue_path, snapshot_path, data = paths('1800', z) 
    data_median, mass_list, bin_item = data_bin(data['M500c'], data['host_id'], bin_num=20, stats=False)

    cc = cat.catalogue(catalogue_path, delta, delta_200)
    gal_lums = cc.subhalo_luminosities[:,2]
    halo_id_list= bin_item[-1]
    
    #find index of clusters for the largest mass bin halos
    for i in range(0,6):
        halo_id = halo_id_list[i]
        print(halo_id)
        idx_cluster=int(halo_id-1)
        idx_gals, mag_gals, pos = locate_gals(cc, idx_cluster, num=4)
        print(cc.CoP[idx_cluster])
        print(cc.M500c[idx_cluster])
        plot_density_map(idx_cluster,idx_gals, cc.CoP, cc.R200c, snapshot_path, z, mag_gals, gal_lums, pos)
