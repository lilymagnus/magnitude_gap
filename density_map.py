import numpy as np
import h5py as h5
import unyt
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import yt
from utils import snapshot
import unyt
from swiftsimio import mask
from swiftsimio import load
import cat_reader as cat
from swiftsimio.visualisation.smoothing_length_generation import generate_smoothing_lengths
from swiftsimio.visualisation.projection import scatter
#choose halo with small and large magnitude gap

def halo_finder(gap, z):
    #look through presaved dataset 
    ds = yt.load("r_mag_gap/3600/m8_3600_z"+str(z)+"_500c_1e13_50kpc.h5")
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

def plot_density_map(clust_idx, idx_gals,CoP, haloRadius, flamingoFile, z, mag_gals, gal_lums):
    #massChoice=cat.haloMass[clust_idx]
    rChoice=haloRadius[clust_idx]*(1.+z)
    xChoice=CoP[clust_idx,0]
    yChoice=CoP[clust_idx,1]
    zChoice=CoP[clust_idx,2]
    print(clust_idx)
    xCen = unyt.unyt_quantity(xChoice,'Mpc')
    yCen = unyt.unyt_quantity(yChoice,'Mpc')
    zCen = unyt.unyt_quantity(zChoice,'Mpc')
    maxRegion = unyt.unyt_quantity(1*rChoice,'Mpc')
    maskRegion = mask(flamingoFile)

    #spatially mask the snapshot data around the cluster                                                                                                                                                     
    region=[[xCen-maxRegion,xCen+maxRegion],
            [yCen-maxRegion,yCen+maxRegion],
            [zCen-maxRegion,zCen+maxRegion]]
    
    maskRegion.constrain_spatial(region)
    data = load(flamingoFile,mask=maskRegion)
    
    dx=data.stars.coordinates.value[:,0]-xChoice
    dy=data.stars.coordinates.value[:,1]-yChoice
    dz=data.stars.coordinates.value[:,2]-zChoice
    h= generate_smoothing_lengths(
        data.stars.coordinates,
        data.metadata.boxsize,
        kernel_gamma=1.8,
        neighbours=57,
        speedup_fac=2,
        dimension=3
    )
    m=data.stars.masses.value
    #t=data.gas.temperatures.value
    #d=data.gas.densities.value
    '''
    h = generate_smoothing_lengths(data.dark_matter.coordinates,
                                   data.metadata.boxsize,
                                   1.8,
                                   neighbours=57)
    
    dx=data.dark_matter.coordinates.value[:,0]-xChoice                                                                                                                                                   
    dy=data.dark_matter.coordinates.value[:,1]-yChoice                                                                                                                                                  
    dz=data.dark_matter.coordinates.value[:,2]-zChoice                                                                                                                                                  
    m=data.dark_matter.masses.value                                                                                                                                                                      
    '''
    ind=np.where((dx>-maxRegion.value)&(dx<maxRegion.value)&
             (dy>-maxRegion.value)&(dy<maxRegion.value)&
             (dz>-maxRegion.value)&(dz<maxRegion.value))[0]
    #put within bounds [0,1] 
    dx=(dx[ind]+maxRegion.value)/(2.*maxRegion.value)
    dy=(dy[ind]+maxRegion.value)/(2.*maxRegion.value)
    dz=(dz[ind]+maxRegion.value)/(2.*maxRegion.value)
    h=h[ind]/(2.*maxRegion.value)
    m=m[ind]
    #t=t[ind]
    #d=d[ind]

    mapRes = 600    
    fig,axs = plt.subplots(figsize=(6,6))
    
    map=scatter(x=dx,y=dy,h=h,m=m,res=mapRes) 
    map+=1.e-10
    map/=map.max()
    mapMin=-4
    mapMax=0
    mapLabel=r'$\log_{10}(\Sigma_{\rm *}/\Sigma_{\rm *,max})$'
    #map=scatter(x=dx,y=dy,h=h,m=m,res=mapRes)    
    image=axs.pcolormesh(np.log10(map),cmap='magma',vmin=mapMin,vmax=mapMax)
    
    axs.axis('off')
    axs.set_box_aspect(1) 
    axs.add_patch(Circle((0.5*mapRes,0.5*mapRes),rChoice*mapRes/(2.*maxRegion.value),fill=False,color='yellow',linewidth=2))
    #VR ID is just the index+1 since the VR_ID list starts at 1 
    cs = ['red','green','yellow']
    for i in range(4):
        rsub= unyt.unyt_quantity(50,'kpc')
        xsub=CoP[idx_gals[i],0] - xChoice
        ysub=CoP[idx_gals[i],1] - yChoice
        xsub = (xsub +maxRegion.value)/(2.*maxRegion.value)
        ysub = (ysub +maxRegion.value)/(2.*maxRegion.value)

        #scatter image is inverted so flip the subhalo x and y coords
        axs.add_patch(Circle(( ysub*mapRes,  xsub*mapRes), rsub.to('Mpc')*mapRes/(2.*maxRegion.value),fill=False,linewidth=2, color='orange'))
        axs.annotate('%.2f' % mag_gals[i], (((0.01 + ysub)*mapRes), ((+0.02 + xsub)*mapRes)), color='orange')
 
    bgg_subx, bgg_suby, bgg_rsub = 0.5, 0.5,  unyt.unyt_quantity(50,'kpc')
    axs.add_patch(Circle((bgg_subx*mapRes, bgg_suby*mapRes),bgg_rsub.to('Mpc')*mapRes/(2.*maxRegion.value),fill=False,color='blue',linewidth=2))
    axs.annotate('%.2f' % (-2.5*np.log10(gal_lums[clust_idx])), (((0.01 + bgg_suby)*mapRes), ((-0.02 + bgg_subx)*mapRes)),color='blue')
    
    fig.tight_layout()  
    plt.savefig('star_density.png',dpi=300)
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
    
    return catalogue_path, snapshot_path

if __name__ == "__main__":
    z = 0
    #cp_18 = paths('1800', 0)
    catalogue_path, snapshot_path = paths('1800', 0) 
    cc_36 = cat.catalogue(catalogue_path, delta, delta_200)
    gal_lums = cc_36.subhalo_luminosities[:,2]

    #identify the VR ID and index of the host halo
    #halo_id = halo_finder(1.5,z)
    #halo_id = cc_36.VR_ID[np.argmax(cc_36.M500c)]
    #halo_id = 71444974
     
    cent = [567.7116363, 148.0196063, 138.0276663] - cc_36.CoP
    idx_cluster = np.argmin(np.sum(cent**2, axis=1) )
    #idx_cluster = np.argmin(np.abs(cc_36.M500c - 10**14))
    
    if cc_36.host_ID[idx_cluster] != -1:
        idx_cluster = cc_36.host_ID[idx_cluster] - 1

    halo_id = cc_36.VR_ID[idx_cluster]
    
    print(cc_36.CoP[idx_cluster])
    #find index of subhalos brightest galaxy
    sub_lums = gal_lums[cc_36.host_ID == halo_id]
    sub_VRs = cc_36.VR_ID[cc_36.host_ID == halo_id]
    sub_CoP = cc_36.CoP[cc_36.host_ID == halo_id]
    
    idx_gals = []
    mag_gals = []
    pos = []
    
    cluster_radius = cc_36.R200c[idx_cluster]
    cluster_CoP = cc_36.CoP[idx_cluster]
    dist_to_host_cop = cc_36.CoP  - cluster_CoP 
    r_sqrd = np.sum(dist_to_host_cop**2, axis=1)
    
    subs_ids = cc_36.VR_ID[np.where(r_sqrd < cluster_radius**2)[0]]
    lums = gal_lums[np.where(r_sqrd < cluster_radius**2)[0]]
    sub_CoP = cc_36.CoP[np.where(r_sqrd < cluster_radius**2)[0]]
    ''' 
    for i in range(20):                                                                                                                                                                                  
        j = i +1                                                                                                                                                                                         
        VR_gals = sub_VRs[np.argsort(lums)[-j]]                                                                                                                                                      
        idx_gals.append(np.where(VR_ID == VR_gals)[0])                                                                                                                                                   
        mag_gals.append(-2.5*np.log10(np.sort(sub_lums)[-j])) 
        pos.append(sub_CoP[np.argsort(lums)[-j]])
    
    dist = [669.428,808.23,707.261] - CoP
    rad = np.sum(dist**2, axis=1)
    R = unyt.unyt_quantity(50,'kpc')
    
    p = CoP[np.where(rad < (R.to('Mpc'))**2)[0]]
    
    '''
    breakpoint()
    for i in range(4):
        j = i +1
        VR_gals = sub_VRs[np.argsort(lums)[-j]]
        idx_gals.append(np.where(cc_36.VR_ID == VR_gals)[0])
        mag_gals.append(-2.5*np.log10(np.sort(lums)[-j]))
        pos.append(sub_CoP[np.argsort(lums)[-j]])
    
    plot_density_map(idx_cluster,idx_gals, cc_36.CoP, cc_36.R200c, snapshot_path, z, mag_gals, gal_lums)
