import numpy as np
import h5py as h5
import unyt
import matplotlib.pyplot as plt
import yt
import scipy 
from scipy.optimize import curve_fit

def mass_bin(hmass, mag_gap):
     print(np.max(hmass))
     #sort and bin data in terms of halo masses                                                                                                                                                          
     array = tuple(zip(hmass, mag_gap))
     #sorting the tuple in order of mass values                                                                                                                                                          
     sorted_list = sorted(array, key=lambda x:x[0])
     mass_sort, mag_sort = zip(*sorted_list)
     
     mass_list = list(np.log10(mass_sort))
     mag_list = list(mag_sort)
     #m_copy = mass_list
     bins = np.linspace(np.min(mass_list), np.max(mass_list), 20)

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
          #perc_16.append(np.percentile(mag_section,16))                                                                                                                                              
          #perc_84.append(np.percentile(mag_section,84))                                                                                                                                                 
          #boot_low.append(boot[0])                                                                                                                                                  
          #boot_high.append(boot[1])                                                                                                                                                                    
          
          #     breakpoint()
          del mag_list[0:n[i]]
          
          
     return mag_gap_median, mass_list, bins

def func(x, a, b):
     #choose sigmoid curve
     return 1/(1 + a*np.exp(-x*b))

if __name__ == "__main__":
     mgap = []
     z =0
     cs = ['blue','purple','orange','red']
     #for z in range(0,3):
     i = 0
     for res in ['3600','1800']:
          ds = yt.load("r_mag_gap/" + res + "/m8_"+ res + "_z%s" % z + "_500c_1e13_50kpc.h5")
          data = ds.data     
          hmass = data['M500c']
          mag_gap = data['M14']
          
          bgg = -2.5 * np.log10(data['lum_bgg'])
          bgg_4th = -2.5 * np.log10(data['lum_4thbgg'])
          bgg_2nd = -2.5 * np.log10(data['lum_2ndbgg'])
          '''
          plt.scatter(np.log10(hmass), bgg, color='r', s=10, label='bcg')
          plt.scatter(np.log10(hmass), bgg_4th, color='b', s=10, label='4th bcg')
          #plt.scatter(mag_gap, bgg, color='r', s=10, label='bcg')
          #plt.scatter(mag_gap, bgg_4th, color='b', s=10, label='4th bcg')
          plt.xlabel('Halo Mass')
          #plt.xlabel('M14')
          plt.ylabel('r-band magnitude')
          plt.legend()
          plt.savefig('bcg_mag_mass.png')
          breakpoint()
          '''
          mag, m, bs = mass_bin(hmass, bgg)
          mag_median, mass_list, bins = mass_bin(hmass, mag_gap)
          x = np.linspace(np.min(mass_list), np.max(mass_list), len(bins)-1)
          #plt.fill_between(x, perc_16, perc_84,alpha=0.5)
          plt.scatter(x,mag_median, marker='x',label=res)
          i += 1
          #plt.scatter(x,mag, color=cs[i], marker='x',label=res + ', 4th bgg')
          i += 1
          mgap.append(mag_median)
     #plt.plot(x, func(x, param[0], param[1]))
     
     #asymmetric_error = np.array(list(zip((np.array(median)-np.array(boot_low)), 
     #                                     (np.array(boot_high)-np.array(median))))).T
     #plt.errorbar(x, median, yerr=asymmetric_error, fmt='.', ecolor = 'red')
     #del ds
     
     #param, _ = curve_fit(func, x, mgap, p0=[2,3])
     plt.ylabel('r-band magnitude')
     plt.xlabel('Halo Mass [Msun]')
     plt.legend()
     plt.savefig("redshift_mag_gap.png")

     breakpoint()
