import numpy as np
import h5py as h5
import unyt
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import yt
import unyt
from swiftsimio.visualisation.slice import slice_gas
from swiftsimio import mask
from swiftsimio import load
from swiftsimio.visualisation.smoothing_length_generation import generate_smoothing_lengths
from swiftsimio.visualisation.projection import scatter
#choose halo with small and large magnitude gap
def halo_finder(gap):
    #look through presaved dataset 
    ds = yt.load("r_mag_gap_data_z0.h5")
    data = ds.data

    #remove all bad data
    bgg_4th = data['4thBGG_mag'][data['host_mass'] > 0]
    bgg = data['halo_mag'][data['host_mass'] > 0]
    halo_id = data['halo_id'][data['host_mass'] > 0]

    print('finding halo id...')
    for i in range(len(halo_id)):
        mag_gaps =  np.abs(bgg[i] - bgg_4th[i])
        if np.abs(mag_gaps - gap) <= 0.05:
            return halo_id[i]

    raise Exception ('no magnitude gap near what has been input')


def get_sub_date(idx_4th_BCG, haloRadius, z, xc, yc, maxRegion):
    #need to normalize and recenter wrt to center of cluster
    rsub=haloRadius[idx_4th_BCG]*(1.+z)
    xsub=CoP[idx_4th_BCG,0] - xc
    ysub=CoP[idx_4th_BCG,1] - yc    
    breakpoint()
    return (xsub +maxRegion.value)/(2.*maxRegion.value),(ysub +maxRegion.value)/(2.*maxRegion.value), rsub
     
def plot_density_map(clust_idx, idx_4th_BCG,CoP, haloRadius, flamingoFile, z):
    #massChoice=cat.haloMass[clust_idx]
    rChoice=haloRadius[clust_idx]*(1.+z)
    xChoice=CoP[clust_idx,0]
    yChoice=CoP[clust_idx,1]
    zChoice=CoP[clust_idx,2]
    print(clust_idx)
    xCen = unyt.unyt_quantity(xChoice,'Mpc')
    yCen = unyt.unyt_quantity(yChoice,'Mpc')
    zCen = unyt.unyt_quantity(zChoice,'Mpc')
    maxRegion = unyt.unyt_quantity(2*rChoice,'Mpc')
    maskRegion = mask(flamingoFile)

    #spatially mask the snapshot data around the cluster                                                                                                                                                     
    region=[[xCen-maxRegion,xCen+maxRegion],
            [yCen-maxRegion,yCen+maxRegion],
            [zCen-maxRegion,zCen+maxRegion]]
    
    maskRegion.constrain_spatial(region)
    data = load(flamingoFile,mask=maskRegion)
    """
    dx=data.gas.coordinates.value[:,0]-xChoice
    dy=data.gas.coordinates.value[:,1]-yChoice
    dz=data.gas.coordinates.value[:,2]-zChoice
    h=data.gas.smoothing_lengths.value
    m=data.gas.masses.value
    t=data.gas.temperatures.value
    d=data.gas.densities.value
    """
    h = generate_smoothing_lengths(data.dark_matter.coordinates,
                                   data.metadata.boxsize,
                                   1.8,
                                   neighbours=57)
    
    dx=data.dark_matter.coordinates.value[:,0]-xChoice                                                                                                                                                   
    dy=data.dark_matter.coordinates.value[:,1]-yChoice                                                                                                                                                  
    dz=data.dark_matter.coordinates.value[:,2]-zChoice                                                                                                                                                  
    m=data.dark_matter.masses.value                                                                                                                                                                      
    
    ind=np.where((dx>-maxRegion.value)&(dx<maxRegion.value)&
             (dy>-maxRegion.value)&(dz<maxRegion.value)&
             (dz>-maxRegion.value)&(dz<maxRegion.value))[0]
    #put within bounds [0,1] 
    dx=(dx[ind]+maxRegion.value)/(2.*maxRegion.value)
    dy=(dy[ind]+maxRegion.value)/(2.*maxRegion.value)
    dz=(dz[ind]+maxRegion.value)/(2.*maxRegion.value)
    h=h[ind]/(2.*maxRegion.value)
    m=m[ind]
    #t=t[ind]
    #d=d[ind]

    mapRes = 512    
    fig,axs = plt.subplots(figsize=(6,6))
    map=scatter(x=dx,y=dy,h=h,m=m,res=mapRes)    
    image=axs.pcolormesh(np.log10(map),cmap='magma')
    axs.axis('off')
    axs.set_box_aspect(1)
    axs.add_patch(Circle((0.5*mapRes,0.5*mapRes),rChoice*mapRes/(2.*maxRegion.value),fill=False,color='yellow',linewidth=2))
    subx, suby, rsub = get_sub_date(idx_4th_BCG, haloRadius, z, xChoice, yChoice, maxRegion)
    breakpoint()
    axs.add_patch(Circle((subx*mapRes, suby*mapRes),rsub*mapRes/(2.*maxRegion.value),fill=False,color='green',linewidth=2))
    fig.tight_layout()
    plt.savefig('density.png',dpi=300)
    
group_path = "/cosma8/data/dp004/flamingo/Runs/"
box = "L1000N1800"
run = "HYDRO_FIDUCIAL"
run_name = "HF"
data_path = group_path + box + "/" + run
catalogue = "/SOAP"
delta="500_crit" # density contrast [200_mean, 500_crit etc]

if __name__ == "__main__":
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
    
    h5dset = h5file['VR/HostHaloID']
    #list of host halo ids for each subhalo, -1 for hosts
    host_ID = h5dset[...]
    
    h5dset = h5file['InclusiveSphere/50kpc/StellarLuminosity']
    gal_lums = h5dset[...]
    #take r-band only
    gal_lums = gal_lums[:,2]
    
    h5dset = h5file['/SO/' + delta + '/SORadius']
    haloRadius = h5dset[...]
    
    h5dset = h5file['/VR/CentreOfPotential']
    CoP = h5dset[...]
    h5file.close()
    
    #identify the VR ID and index of the host halo
    halo_id = halo_finder(gap=2.5)
    idx_cluster = np.where(VR_ID == halo_id)[0]
     
    #VR IDs of subhalos
    sub_VR = VR_ID[host_ID == halo_id] - 1 

    #find index of subhalo with 4th brightest galaxy
    sub_lums = gal_lums[host_ID == halo_id]
    idx_4th_BCG = np.argsort(sub_lums)[-4]
    breakpoint()
    plot_density_map(idx_cluster,sub_VR[idx_4th_BCG], CoP, haloRadius, snapshot_path, z=0)
