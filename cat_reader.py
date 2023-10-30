
import numpy as np
import h5py as h5

class catalogue:
    def __init__(self,filename, delta_500, delta_200):
        h5file=h5.File(filename,'r')
        groupName="/SWIFT/Cosmology/"
        #choose the group, in this case Cosmology is the type of cosmology framework setup                                                                                                                 
        h5group=h5file[groupName]

        #setup the attributed of the simualtion --> LCDM parameters                                                                                                                                         
        attrName="Critical density [internal units]"
        self.rhoCrit=h5group.attrs[attrName]
        attrName="H [internal units]"
        attrName="Omega_lambda"
        self.OmegaLambda=h5group.attrs[attrName]
        self.Hz=h5group.attrs[attrName]
        attrName="H0 [internal units]"
        self.H0=h5group.attrs[attrName]
        self.Ez=self.Hz/self.H0
        attrName="Omega_b"
        self.OmegaBaryon=h5group.attrs[attrName]
        attrName="Omega_m"
        self.OmegaMatter=h5group.attrs[attrName]
        attrName="Redshift"
        self.redshift=h5group.attrs[attrName]
        attrName="h"
        self.hubbleParameter=h5group.attrs[attrName]

        #open the rest of the grouped simulations
        h5file = h5.File(filename, 'r')
        h5dset = h5file['/SO/'+ delta_500 + '/TotalMass']
        self.M500c = h5dset[...]
         
        h5dset = h5file['/SO/'+ delta_200 + '/TotalMass']
        self.M200c = h5dset[...]
        
        h5dset = h5file['VR/ID']
        self.VR_ID = h5dset[...]
        
        h5dset = h5file['/SO/200_crit/SORadius']
        self.R200c = h5dset[...] #Mpc                                                                                                                                                                      
        
        h5dset = h5file['InclusiveSphere/50kpc/StellarMass']
        self.stellar_mass = h5dset[...]
        
        h5dset = h5file['InclusiveSphere/50kpc/StellarLuminosity']
        self.subhalo_luminosities = h5dset[...]
        
        #all subhalos have the same host id, halos that are hosts/empty = -1                                                                                                                                    
        h5dset = h5file['VR/HostHaloID']
        self.host_ID = h5dset[...]
        
        h5dset = h5file['InclusiveSphere/50kpc/TotalMass']
        self.aperture_mass = h5dset[...]
        
        h5dset = h5file['/VR/CentreOfPotential']
        self.CoP = h5dset[...]
        
        h5file.close()
