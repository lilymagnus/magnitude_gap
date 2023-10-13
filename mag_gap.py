import numpy as np
import h5py as h5
import unyt
import matplotlib.pyplot as plt
import yt
import scipy 

if __name__ == "__main__":
     ds = yt.load("m8_1800_r_mag_gap_data_z0.h5")
     data = ds.data
     M500c = data['M500c']
     mag_gap = data['mag_gap']

     #sort and bin data in terms of halo masses
     array = tuple(zip(M500c, mag_gap))
     #sorting the tuple in order of mass values
     sorted_list = sorted(array, key=lambda x:x[1])
     mass_sort, mag_gap_sort = zip(*sorted_list)

     mass_list = list(mass_sorted)
     mag_gap_list = list(mag_gap_sort)
     
     bins = np.logspace(mass_sort[0], mass_sort[-1], 40)
     
     n,_=np.histogram(mass_list,bins=bins)

     median = []
     boot_low = []
     boot_high
     for i in range(len(n)):
          mag_gap_list = mass_list[0:n[i]]
          mean.append(np.median(mag_gap_list))
          boot = scipy.stats.bootstrap(mag_gap_list, np.median, confidence_level=0.95,
                         random_state=1, method='percentile').confidence_interval
          boot_low.append(boot[0])
          boot_high.append(boot[1])
          del mass_list[0:n[i]]
          
     
     plt.fill_between(bins, np.percentile(median,16), np.percentile(median,84),alpha=0.5)
     plt.scatter(bins, median) 
     plt.savefig("mag_gap.png")

     breakpoint()
