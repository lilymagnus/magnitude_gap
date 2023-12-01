import h5py as h5
import unyt
from scipy import spatial
import numpy as np
import yt
from yt.funcs import get_pbar
from swiftsimio.visualisation.smoothing_length_generation import generate_smoothing_lengths
from swiftsimio import mask
from swiftsimio import load

def snapshot(z, data_path, catalogue, res):
    start = 0
    print('...locating redshift..')
    if (res == "3600") or (res =='5040'):
        start = 78
    elif (res == "1800") or (res == "0900"):
        start = 77

    for i in range(start,0,-1):
        if (res == "3600") and i == 58:
            continue
        if (res == '5040')  and (i not in [78,76,72,68,58,48,38,18,14,10]):
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
        #snaps divided as z = 0.95 and 1.05 etc...                                                                                                                              
        if round(float(redshift),2) == z:
            return i

    raise Exception('Redshift not found...')



def DM_map(clust_idx,CoP, haloRadius, flamingoFile, z):
    #massChoice=cat.haloMass[clust_idx]                                                                                                                                                                                                                                        
    rChoice=haloRadius[clust_idx]*(1.+z)
    xChoice=CoP[clust_idx,0]
    yChoice=CoP[clust_idx,1]
    zChoice=CoP[clust_idx,2]

    xCen = unyt.unyt_quantity(xChoice,'Mpc')
    yCen = unyt.unyt_quantity(yChoice,'Mpc')
    zCen = unyt.unyt_quantity(zChoice,'Mpc')
    maxRegion = unyt.unyt_quantity(0.5*rChoice,'Mpc')
    maskRegion = mask(flamingoFile)

    #spatially mask the snapshot data around the cluster                                                                                                                                                                                                                       
    region=[[xCen-maxRegion,xCen+maxRegion],
            [yCen-maxRegion,yCen+maxRegion],
            [zCen-maxRegion,zCen+maxRegion]]

    maskRegion.constrain_spatial(region)
    data = load(flamingoFile,mask=maskRegion)

    dx=data.dark_matter.coordinates.value[:,0]-xChoice
    dy=data.dark_matter.coordinates.value[:,1]-yChoice
    dz=data.dark_matter.coordinates.value[:,2]-zChoice
    
    h = generate_smoothing_lengths(
        data.dark_matter.coordinates,
        data.metadata.boxsize,
        kernel_gamma=1.8,
        neighbours=57,
        speedup_fac=2,
        dimension=3
    )

    m=data.dark_matter.masses.value

    ind=np.where((dx>-maxRegion.value)&(dx<maxRegion.value)&
             (dy>-maxRegion.value)&(dy<maxRegion.value)&
             (dz>-maxRegion.value)&(dz<maxRegion.value))[0]
    #put within bounds [0,1]                                                                                                                                                                                                                                                   
    dx=(dx[ind]+maxRegion.value)/(2.*maxRegion.value)
    dy=(dy[ind]+maxRegion.value)/(2.*maxRegion.value)
    dz=(dz[ind]+maxRegion.value)/(2.*maxRegion.value)
    h=h[ind]/(2.*maxRegion.value)
    m=m[ind]

    return dx, dy, dz, h, m, rChoice

def star_map(clust_idx,CoP, haloRadius, flamingoFile, z):    
    #massChoice=cat.haloMass[clust_idx]                                                                                                                                                                                                                                        
    rChoice=haloRadius[clust_idx]*(1.+z)
    xChoice=CoP[clust_idx,0]
    yChoice=CoP[clust_idx,1]
    zChoice=CoP[clust_idx,2]

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
    
    ind=np.where((dx>-maxRegion.value)&(dx<maxRegion.value)&
             (dy>-maxRegion.value)&(dy<maxRegion.value)&
             (dz>-maxRegion.value)&(dz<maxRegion.value))[0]
    #put within bounds [0,1]                                                                                                                                                                                                                                                    
    dx=(dx[ind]+maxRegion.value)/(2.*maxRegion.value)
    dy=(dy[ind]+maxRegion.value)/(2.*maxRegion.value)
    dz=(dz[ind]+maxRegion.value)/(2.*maxRegion.value)
    h=h[ind]/(2.*maxRegion.value)
    m=m[ind]

    return dx, dy, dz, h, m, rChoice

def data_bin(data2, data1, bin_num, stats):
    #sort and bin data in terms of other data set
    #e.g bin mass data (data2) in terms of mag gap (data1)

    array = tuple(zip(data2, data1)) 
    #sorting the tuple in order of data2 values                                                                                                                                                     
    sorted_list = sorted(array, key=lambda x:x[0])
    d2_sort, d1_sort = zip(*sorted_list)
    
    d2_list = list(d2_sort)
    d1_list = list(d1_sort)
                                                                                                                                                                                                                                                           
    bins = np.linspace(np.min(d2_list), np.max(d2_list), bin_num)
    n,_=np.histogram(d2_list,bins=bins)
    
    median = []
    bin_item = []
    boot_low = []
    boot_high = []
    perc_16 = []
    perc_84 = []
    pbar = yt.get_pbar('Binning data data', len(n))
    for i in range(len(n)):
        pbar.update(i)
        section = list(np.array(d1_list[0:n[i]]).flatten())
        #m = list(np.array(m_copy[0:n[i]]).flatten())                                                                                                                                                                                                                        
        median.append(np.median(section))
        bin_item.append(section)
        if stats == True:
            boot = scipy.stats.bootstrap((section,), np.median,  
            confidence_level=0.95, method='percentile').confidence_interval           
            perc_16.append(np.percentile(section,16))             
            perc_84.append(np.percentile(section,84))             
            boot_low.append(boot[0])                                 
            boot_high.append(boot[1])                                                                                                                                                                                                                                           
        del d1_list[0:n[i]]

    if stats == False:
        return median, d2_list, bin_item

def paths(res, z, group_path, catalogue, run, size, mass, type_):
    box = size + res
    data_path = group_path + box + "/" + run

    snapshot_ID = snapshot(z, data_path, catalogue, res)
    halo_props_name = "/halo_properties_" + str(snapshot_ID).zfill(4) + ".hdf5"
    catalogue_path = data_path + catalogue + halo_props_name

    snapshot_name = "/flamingo_" + str(snapshot_ID).zfill(4)
    snapshot_path = data_path + "/snapshots" + snapshot_name + snapshot_name + ".hdf5"
    
    #type: acc_200m or mag_500c
    if (res == '1800') or (res == '5040') :
        ds = yt.load("saved_data/" + size + '/' + res + "/m9_" + type_ + '_' + mass +"_z"+str(z)+"_1e14_50kpc.h5")
    else:                         
        ds = yt.load("saved_data/" + size + '/' + res + "/m8_" + type_ + '_'+ mass + "_z"+str(z)+"_1e14_50kpc.h5")

    return catalogue_path, snapshot_path, ds

def region(cent, pos_list, z, radius):
    #only include subs that are within R200c                       
    dist_to_host_cop = cent - pos_list #Mpc               
    R = radius * (1 + z) #Mpc                                                                                                                                                                      
    r_sqrd = np.sum(dist_to_host_cop**2, axis=1)
    return np.where(r_sqrd < (R)**2)[0]


def locate_gals(cc, idx_cluster, z):
    halo_id = cc.VR_ID[idx_cluster]

    print(cc.CoP[idx_cluster])
    gal_lums = cc.subhalo_luminosities[:,2]

    #find index of subhalos brightest galaxy                                                                                                                                                                                                                                   
    sub_lums = gal_lums[cc.host_ID == halo_id]
    sub_VRs = cc.VR_ID[cc.host_ID == halo_id]
    sub_CoP = cc.CoP[cc.host_ID == halo_id]
    sub_star_mass = cc.stellar_mass[cc.host_ID == halo_id]
    sub_mass = cc.bound_subhalo_mass[cc.host_ID == halo_id]
    

    sub_ids = sub_VRs[region(cc.CoP[idx_cluster], sub_CoP, z, cc.R200c[idx_cluster])]
    lums = sub_lums[region(cc.CoP[idx_cluster], sub_CoP, z, cc.R200c[idx_cluster])]
    sub_pos = sub_CoP[region(cc.CoP[idx_cluster], sub_CoP, z, cc.R200c[idx_cluster])]
    star_mass =  sub_star_mass[region(cc.CoP[idx_cluster], sub_CoP, z, cc.R200c[idx_cluster])]
    tot_mass = sub_mass[region(cc.CoP[idx_cluster], sub_CoP, z, cc.R200c[idx_cluster])]
    
    idx = sub_ids[np.argsort(lums[lums>0])][np.where(np.sort(np.log10(lums[lums>0])) < 7.14)[0]]
    
    if len(idx) != 0:
        idx = idx - 1
        idx = int(idx)
        print(cc.VR_ID[idx])
        print(cc.host_ID[idx])
        print(cc.stellar_mass[idx])
        print(cc.CoP[idx])
        print(np.log10(cc.subhalo_luminosities[idx,2]))
    #breakpoint()                                                                                                                                                                                                                                 
    return sub_ids, lums, sub_pos, star_mass, tot_mass
