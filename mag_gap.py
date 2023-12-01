import numpy as np
import h5py as h5
import unyt
import matplotlib.pyplot as plt
import yt
import scipy 
from scipy.optimize import curve_fit

group_path = "/cosma8/data/dp004/flamingo/Runs/"
run = "HYDRO_FIDUCIAL"
catalogue = "/SOAP"

def mass_bin(hmass, mag_gap):
     #print(np.max(hmass))
     #sort and bin data in terms of halo masses                                                                                                                                                          
     array = tuple(zip(hmass, mag_gap))
     #sorting the tuple in order of mass values                                                                                                                                                          
     sorted_list = sorted(array, key=lambda x:x[0])
     mass_sort, mag_sort = zip(*sorted_list)
     
     mass_list = list(np.log10(mass_sort))
     mag_list = list(mag_sort)
     #m_copy = mass_list
     bins = np.linspace(np.min(mass_list), np.max(mass_list), 15)

     n,_=np.histogram(mass_list,bins=bins)
     
     mag_gap_median = []
     boot_low = []
     boot_high = []
     perc_16 = []
     perc_84 = []
     for i in range(len(n)):
          mag_gap_section = list(np.array(mag_list[0:n[i]]).flatten())
          #m = list(np.array(m_copy[0:n[i]]).flatten())
          mag_gap_median.append(np.median(mag_gap_section))
          #boot = scipy.stats.bootstrap((mag_section,), np.median,                                                                                                                                       
          #  confidence_level=0.95, method='percentile').confidence_interval                                                                                             
          perc_16.append(np.percentile(mag_gap_section,16))                                                                                                                                              
          perc_84.append(np.percentile(mag_gap_section,84))                                                                                                                                                 
          #boot_low.append(boot[0])                                                                                                                                                  
          #boot_high.append(boot[1])                                                                                                                                                                    
           
          #     breakpoint()
          del mag_list[0:n[i]]
          
          
     return mag_gap_median, mass_list, bins, perc_16, perc_84

if __name__ == "__main__":
     mgap = []
     z =0
     cs = ['blue','purple','orange','red', 'green', 'lightgreen']
     
     i = 0
     for res in [('3600','L1000N','m8'),('1800','L1000N','m9'), ('5040','L2800N','m9')]:
          
          ds = yt.load("saved_data/" + res[1] + "/" + res[0] + '/' + res[2] + "_mag_500c_z%s" % z + "_1e14_50kpc.h5")
          data = ds.data     
          hmass = data['M500c']          
          M14 = data['M14']
          '''
          mag, m, bs, p16, p84 = mass_bin(hmass,M14)
          x = np.linspace(np.min(np.log10(hmass)), np.max(np.log10(hmass)), len(bs)-1)
          plt.scatter(x,mag, color=cs[i], marker='x',label=res[1] + res[0],zorder=10)
          plt.fill_between(x,p16,p84, color=cs[i], alpha=0.4)
          '''
          bgg = -2.5 * np.log10(data['lum_bgg'])
          bgg4th = -2.5 * np.log10(data['lum_4thbgg'])
          mag, m, bs,_,_ = mass_bin(hmass, bgg4th)
          mag_median, mass_list, bins,_,_ = mass_bin(hmass, bgg)
          
          x = np.linspace(np.min(mass_list), np.max(mass_list), len(bins)-1)
          plt.scatter(x,mag, color=cs[i], marker='x',label= res[1] + res[0] + ': 4th bgg')
          i += 1
          plt.scatter(x,mag_median, color=cs[i], marker='x',label= res[1] + res[0] +': bgg')
          i += 1
          
          
     #plt.yscale('log')
     plt.ylabel('r-band magnitude gap')
     #plt.ylabel('Luminosity [Lsun]')
     plt.xlabel('Halo Mass [Msun]')
     plt.legend(loc='upper left', bbox_to_anchor=(1.04, 1))
     plt.tight_layout()
     plt.savefig("hmass_mag.png")

     breakpoint()
