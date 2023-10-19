import numpy as np
import h5py as h5
import unyt
import matplotlib.pyplot as plt
import yt
import scipy 

if __name__ == "__main__":
     #ds = yt.load("r_mag_gap/m8_1800_z0_500c_1e13_50kpc.h5")
     ds = yt.load('r_mag_gap/m8_1800_z0_200c_1e14_100kpc.h5')
     data = ds.data
     
     M200c = data['M200c']
     mag_gap = data['mag_gap']
     bgg = -2.5 * np.log10(data['lum_bgg'])
     bgg_4th = -2.5 * np.log10(data['lum_4thbgg'])
     '''
     #plt.scatter(np.log10(M200c), bgg, color='r', s=10)
     plt.scatter(mag_gap, bgg, color='r', s=10, label='bcg')
     plt.scatter(mag_gap, bgg_4th, color='b', s=10, label='4th bcg')
     plt.xlabel('Halo Mass (200crit)')
     plt.ylabel('r-band magnitude')
     plt.legend()
     plt.savefig('bcg_mag_mass.png')
     breakpoint()
     '''
     #sort and bin data in terms of halo masses

     array = tuple(zip(M200c, bgg, bgg_4th))
     #sorting the tuple in order of mass values
     sorted_list = sorted(array, key=lambda x:x[0])
     
     mass_sort, bgg_sort, bgg_4th_sort = zip(*sorted_list)

     mass_list = list(np.log10(mass_sort))
     bgg_list = list(bgg_sort)
     bgg_4th_list = list(bgg_4th_sort)

     bins = np.linspace(np.min(mass_list), np.max(mass_list), 20)
     
     n,_=np.histogram(mass_list,bins=bins)

     bgg_median = []
     bgg_4th_median = []
     boot_low = []
     boot_high = []
     perc_16 = []
     perc_84 = []
     for i in range(len(n)):
          bgg_section = list(np.array(bgg_list[0:n[i]]).flatten())
          bgg_4th_section = list(np.array(bgg_4th_list[0:n[i]]).flatten())
          bgg_median.append(np.median(bgg_section))
          bgg_4th_median.append(np.median(bgg_4th_section))
          #boot = scipy.stats.bootstrap((mag_section,), np.median, 
          #                                  confidence_level=0.95, method='percentile').confidence_interval
          #perc_16.append(np.percentile(mag_section,16))
          #perc_84.append(np.percentile(mag_section,84))
          #boot_low.append(boot[0])
          #boot_high.append(boot[1])
          
          del bgg_list[0:n[i]]
          del bgg_4th_list[0:n[i]]

     x = np.linspace(np.min(mass_list), np.max(mass_list), len(bins)-1)
     #plt.fill_between(x, perc_16, perc_84,alpha=0.5)
     plt.scatter(x, bgg_median, label='bgg')
     plt.scatter(x, bgg_4th_median, label='4th bgg')
     #asymmetric_error = np.array(list(zip((np.array(median)-np.array(boot_low)), 
     #                                     (np.array(boot_high)-np.array(median))))).T
     #plt.errorbar(x, median, yerr=asymmetric_error, fmt='.', ecolor = 'red')
     plt.ylabel('r-band magnitude')
     plt.xlabel('Halo Mass [Msun]')
     plt.legend()
     plt.savefig("200c_bcg_vs_hmass.png")

     breakpoint()
