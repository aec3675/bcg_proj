# coding: utf-8
import argparse

#defining the arguments -- works
parser = argparse.ArgumentParser()
parser.add_argument("data_file", help='input which csv file to use')
parser.add_argument("save_data_dir", help='input where to save the data produced')
args = parser.parse_args()
print('Using csv file :', args.data_file)
print('Saving data to :', args.save_data_dir)



from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os, sys
import numpy as np
import pandas as pd
sys.path.append('/content/gdrive/My Drive/data/')
#!ls '/content/gdrive/My Drive/data/'
import temp_plotting as pl
import shapes_lib as sl

###################################################################################
#NEED TO CHANGE THIS BEFORE RUNNING#
bcgdir = './SPA300/BCG_MAPS_STORED/ALL_CONVERG_MAPS/'
kappa_list = pd.read_csv('/content/drive/My Drive/data/kappa_maps_list.csv', sep=',')
bcg_list = kappa_list['filename_x']
nbcg = len(bcg_list)
print('# of Kappa maps =', nbcg)#,bcg_list)

data = pd.read_csv(args.data_file)
k = data.columns
# print(k)

def get_params(data,cl):
  # #image moments, unweighted column names
  # mp1 = [data[cl][kk] for kk in k[10:15]]
  # mp=mp1+[]
  # mp[2]=mp1[2]/(1+mp1[3]**2)
  # mp[3]=mp1[3]*mp[2]
  # mp[4]=mp[4]+90
  
  # #sersic fit column names
  # sp1 = [data[cl][kk] for kk in k[15:20]] 
  # sp=sp1+[]
  # sp[2]=sp1[2]/np.sqrt(sp1[3])
  # sp[3]=sp1[2]*np.sqrt(sp1[3])
  # sp[4]=sp[4]+90

  #contour 1.0
  cp1 = [data[kk][cl] for kk in k[0:5]] 
  cp_1=cp1+[]
  cp_1[0]=cp1[1]
  cp_1[1]=cp1[0]
  cp_1[2]=cp1[2]/np.sqrt(cp1[3])
  cp_1[3]=cp1[2]*np.sqrt(cp1[3])
  cp_1[4]=cp_1[4]+90

  #contour 1.1
  cp2 = [data[kk][cl] for kk in k[5:10]] 
  cp_11=cp2+[]
  cp_11[0]=cp2[1]
  cp_11[1]=cp2[0]
  cp_11[2]=cp2[2]/np.sqrt(cp1[3])
  cp_11[3]=cp2[2]*np.sqrt(cp1[3])
  cp_11[4]=cp_11[4]+90

  #contour 0.9
  cp3 = [data[kk][cl] for kk in k[10:15]] 
  cp_09=cp3+[]
  cp_09[0]=cp3[1]
  cp_09[1]=cp3[0]
  cp_09[2]=cp3[2]/np.sqrt(cp1[3])
  cp_09[3]=cp3[2]*np.sqrt(cp1[3])
  cp_09[4]=cp_09[4]+90

  return cp_1, cp_11, cp_09

# Get image
bcg_im = bcgdir+bcg_list[0]
#print(bcg_im)
bcg_im = pl.get_image(bcg_im,downscale=4,expand=1)

# Image size in pixels and Mpc
tsy,tsx = bcg_im.shape
mpc = tsx/5. # pix/Mpc
xctr=int(tsx/2); yctr=int(tsy/2)

# Cluster properties
cat='./SPA300/halo_properties/radii_masses_z0.2198.txt'
clname,r500a,r200a = np.loadtxt(cat,unpack=True,dtype=str,usecols=(0,2,3))

# Print some info
kx=k[4]; ky1=k[9]; ky2=k[14]   # table keys to be plotted

# between 1.0 and 1.1
xdiff1 = np.sqrt((data[k[0]]-data[k[5]])**2 + 
                (data[k[1]]-data[k[6]])**2)
padiff1 = ((data[k[4]]+360)%180) - ((data[k[9]]+360)%180) 

#between 1.0 and 0.9
xdiff2 = np.sqrt((data[k[0]]-data[k[10]])**2 + 
                (data[k[1]]-data[k[11]])**2)
padiff2 = ((data[k[4]]+360)%180) - ((data[k[14]]+360)%180)

#difference bw image center and contour
cntrdiff1 = np.sqrt((data[k[0]]-xctr)**2 + 
                (data[k[1]]-yctr)**2)
cntrdiff11 = np.sqrt((data[k[5]]-xctr)**2 + 
                (data[k[6]]-yctr)**2)
cntrdiff09 = np.sqrt((data[k[10]]-xctr)**2 + 
                (data[k[11]]-yctr)**2)

#getting ellipticity (q) measurements
ellip1 = data[k[3]]
ellip11 = data[k[8]]
ellip09 = data[k[13]]

#finding ellip differences
#bw 1 and 1.1
ellip_a = np.abs(ellip1-ellip11) 
#bw 1 and 0.9
ellip_b = np.abs(ellip1-ellip09)


bad_centroid_11 = np.where(xdiff1>20)
bad_centroid_09 = np.where(xdiff2>20)
bad_imgcntr_1 = np.where(cntrdiff1>60)
bad_imgcntr_11 = np.where(cntrdiff11>60)
bad_imgcntr_09 = np.where(cntrdiff09>60)
bad_ellip_1 = np.where(ellip1<0.5)
bad_ellip_11 = np.where(ellip11<0.5)
bad_ellip_09 = np.where(ellip09<0.5)
bad_ellip_2_1 = np.where(ellip1<0.3)
bad_ellip_2_11 = np.where(ellip11<0.3)
bad_ellip_2_09 = np.where(ellip09<0.3)
print('-----')
# print('1.1-BAD POS ANGLES', np.where(np.fabs(data[kx]%180-data[ky1]%180)>20))
print('1.1-BAD CENTROID', bad_centroid_11)
# print('0.9-BAD POS ANGLES', np.where(np.fabs(data[kx]%180-data[ky2]%180)>20))
print('0.9-BAD CENTROID', bad_centroid_09)
print('-----')
print('BAD IMG CNTR-1.0', bad_imgcntr_1)
print('BAD IMG CNTR-1.1', bad_imgcntr_11)
print('BAD IMG CNTR-0.9', bad_imgcntr_09)
print('-----')
print('TOO ELLIP (0.5 CUTOFF) -1.0', bad_ellip_1)
print('TOO ELLIP (0.5 CUTOFF) -1.1', bad_ellip_11)
print('TOO ELLIP (0.5 CUTOFF) -0.9', bad_ellip_09)
print('-----')
print('TOO ELLIP (0.3 CUTOFF) -1.0', bad_ellip_2_1)
print('TOO ELLIP (0.3 CUTOFF) -1.1', bad_ellip_2_11)
print('TOO ELLIP (0.3 CUTOFF) -0.9', bad_ellip_2_09)


#function that pulls out actual index number and appends to flag_list (for ellip cutoff 0.5)
big_bad_flag_list_1 = []
def flag_bad_1(arraylist):
  for badarray in arraylist:
    for num in badarray[0]:
      big_bad_flag_list_1.append(num)
      
#function that pulls out actual index number and appends to flag_list (for ellip cutoff 0.3)
big_bad_flag_list_2 = []
def flag_bad_2(arraylist):
  for badarray in arraylist:
    for num in badarray[0]:
      big_bad_flag_list_2.append(num)

#pulling out bad cl for all bad lists using func
flag_bad_1([bad_centroid_11,bad_centroid_09,bad_imgcntr_09,bad_imgcntr_1,bad_imgcntr_11,bad_ellip_1,bad_ellip_11,bad_ellip_09])

#pulling out bad cl for all bad lists using func
flag_bad_2([bad_centroid_11,bad_centroid_09,bad_imgcntr_09,bad_imgcntr_1,bad_imgcntr_11,bad_ellip_2_1,bad_ellip_2_11,bad_ellip_2_09])


#finding which cl in big_bad are repeated enough to be classified as badbad (for ellip cutoff 0.5)
bad_cluster_list_1 = []
for i in range(0,nbcg):
  freq = big_bad_flag_list_1.count(i)
  if freq >= 3:
    bad_cluster_list_1.append(i)
print('bad cluster indexes (0.5 cutoff):', bad_cluster_list_1)

#finding which cl in big_bad are repeated enough to be classified as badbad (for ellip cutoff 0.3)
bad_cluster_list_2 = []
for i in range(0,nbcg):
  freq = big_bad_flag_list_2.count(i)
  if freq >= 3:
    bad_cluster_list_2.append(i)
print('bad cluster indexes (0.3 cutoff):', bad_cluster_list_2)

#creating flag column for data frame (ellip<0.5)
flaglist1 = np.zeros(nbcg)
for cl in bad_cluster_list_1:
  for j in range(len(flaglist1)):
    if cl == j:
      flaglist1[j] = 1
      
#creating flag column for data frame (ellip<0.3)
flaglist2 = np.zeros(nbcg)
for cl in bad_cluster_list_2:
  for j in range(len(flaglist2)):
    if cl == j:
      flaglist2[j] = 1
# print(flaglist)

data['flags (ellip<0.5)'] = flaglist1
data['flags (ellip<0.3)'] = flaglist2

#exporting data+flags as new csv file
data.to_csv(args.save_data_dir+'FLAGGED_kappa_measurements.csv', index=False)


