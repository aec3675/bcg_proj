# coding: utf-8
import configparser
import argparse

#defining the arguments -- works
parser = argparse.ArgumentParser()
parser.add_argument("config_path", help='input which config file to use')
parser.add_argument("save_data_dir", help='input where to save the data produced')
args = parser.parse_args()
print('Using Config File :', args.config_path)
print('Saving data to :', args.save_data_dir)

#reading in config file -- works
config = configparser.ConfigParser()
config.read(args.config_path) #needs to be in same dir as savedir
#print(config.sections())
bcgdir = config['files']['bcgdir']
percentile = int(config['params']['percentile'])
w_guess = int(config['params']['w_guess'])
print('BCGDIR =', bcgdir, type(bcgdir))
print('PTILE =', percentile, type(percentile))
print('W GUESS =', w_guess, type(w_guess))



#COPYING AND PASTING (BASICALLY) THE BCG ELLIP FILE
import numpy as np
from astropy.io import fits
import skimage
from skimage.measure import *
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import patches
import os, sys
from astropy.modeling import models
import photutils
import time
import statmorph
from astropy.visualization import simple_norm
sys.path.append('/content/gdrive/My Drive/data/Ricardo_data') #added in my file path
sys.path.append('gdrive/My Drive/data/Ricardo_data') #added in my file path
#!ls '/content/gdrive/My Drive/data/'
import shapes_lib as sl
import plotting as pl
from astropy.io import ascii
from astropy.table import Table
import multiprocessing
from multiprocessing import Pool

print('I IMPORTED STUFFS')

#reading in all the bcg files
bcg_list = np.array(os.listdir(bcgdir))
nbcg = len(bcg_list)
print('# of BCGs =', nbcg) #,bcg_list)
this_bcg=bcgdir+bcg_list[0]

print('I READ IN THE FILES')


#defining the measure all function which finds all the fits/contours
def measure_all(i):
  # Read in image
  bcg_im3 = fits.getdata(bcgdir+bcg_list[i].replace('maps','bcg'))
  bcg_im = np.zeros((bcg_im3.shape[0]+2,bcg_im3.shape[1]+2),dtype=bcg_im3.dtype)
  bcg_im[1:-1,1:-1]=bcg_im3

  # Image size in pixels and Mpc
  tsy,tsx = bcg_im.shape
  mpc = tsx/0.4 # pix/Mpc
  xctr=int(tsx/2); yctr=int(tsy/2)

  t0=time.time()
  the_res=[] ; colnames=[]

  # Compute moments -- input param w_guess
  # Weighted moments full image----------------------------------------------
  q20,q11,q02 = sl.weighted_higher_order_moments(bcg_im,w_guess,2,xctr,yctr)
  q00 = sl.weighted_higher_order_moments(bcg_im,w_guess,0,xctr,yctr)[0]
  Tchi = np.around(np.sqrt((q20+q02)/q00),decimals=2)
  size_a,size_b,pos_angle = sl.axes_lengths_from_chi_moments(Tchi,(q20-q02)/(q20+q02),2*q11/(q20+q02))
  # Save results
  the_res += [xctr,yctr,Tchi,size_b/size_a,pos_angle+90]
  colnames += ['moment_wghtd_'+str(w_guess)+'_xc','moment_wghtd_'+str(w_guess)+'_yc','moment_wghtd_'+str(w_guess)+'_T','moment_wghtd_'+str(w_guess)+'_q','mom_wgh_'+str(w_guess)+'_pos_ang']

  mom_p=((xctr,yctr),Tchi, 
         np.around((q20-q02)/(q20+q02),decimals=4),
         np.around((2*q11)/(q20+q02),decimals=4))

  # Fit Sersic-----------------------------------------------------
  fit_morph = sl.fit_sersic(bcg_im,verbose=False,plot_mod=False)
  q_sersic= (1-fit_morph.sersic_ellip)/(1+fit_morph.sersic_ellip)
  # Save results
  if q_sersic==-1:
    the_res+=[-1,-1,-1,q_sersic,-1]
  else:
    the_res+=[fit_morph.sersic_xc,fit_morph.sersic_yc,fit_morph.sersic_rhalf,
                        q_sersic,fit_morph.sersic_theta*180/np.pi + 90]
  colnames+=['sersic_fit_xc','sersic_fit_yc','sersic_fit_r','sersic_fit_q','sersic_fit_pos_ang']


  # Measure (longest) contour and its properties -- ((percentile variable used here))-----------
  longcon,longcon_props = sl.get_longest_contour_props(bcg_im,contour_intensity=np.percentile(bcg_im,percentile))
  a_c,b_c = longcon_props[0][0].major_axis_length,longcon_props[0][0].minor_axis_length
  theta_c = 180 - longcon_props[0][0].orientation*180/np.pi
  # Save results
  if a_c==0: 
    the_res+=[longcon_props[0][0].centroid[0],longcon_props[0][0].centroid[1],np.sqrt(b_c*a_c),-1, -99,]
  else: 
    the_res+=[longcon_props[0][0].centroid[0],longcon_props[0][0].centroid[1],np.sqrt(b_c*a_c),b_c/a_c, theta_c]
  colnames+=['contour_'+str(percentile)+'ptile_xc','contour_'+str(percentile)+'ptile_yc','contour_'+str(percentile)+'ptile_r','contour_'+str(percentile)+'ptile_q','contour_'+str(percentile)+'ptile_pos_ang']
  

  #make visual check images -----------------------------------------------
  pl.plot_img(bcg_im,
              mom_p,fit_morph,
              contour=longcon, contour_props=longcon_props,
              zoomins=[0,100,150],
              save=args.save_data_dir+'TEST_shapes_%s.pdf' %(bcg_list[i]),
              show=False)
  print(i,bcg_list[i], time.time()-t0)

  return the_res,colnames



#defining the writing routine -- (maybe switch to pd df.to_csv eventually?)
def write_it(results,measured_bcg_list,colnames=None):
  # Ensured above that all position angles are using the same definition (I think/hope)
  if colnames is None:
    colnames=['moment_xc','moment_yc','moment_T','moment_q','mom_pos_ang',
              'contour_xc','contour_yc','contour_r','contour_q','contour_pos_ang',
              'sersic_fit_xc','sersic_fit_yc','sersic_fit_r','sersic_fit_q','sersic_fit_pos_ang']
  data = Table({cn: results[:,i] for i,cn in enumerate(colnames)})
  data['cluster']=measured_bcg_list
  clnr=[bb.split('_')[1]+'.'+bb.split('_')[3][0] for bb in measured_bcg_list]
  data['clusternr']=clnr
  ascii.write(data,args.save_data_dir+'TEST_BCG_shape_measurements.csv',overwrite=True,format='csv')


#expands measure_all to all bcg indices
def do_all(bcg_indices):
  results=None
  measured_bcg_list=np.array(nbcg*['O'*len(bcg_list[0])],dtype=str)
  for i in bcg_indices:
    res,cols = measure_all(i)
    if results is None:
      results=np.zeros((nbcg,len(res)))
    results[i] = res
    measured_bcg_list[i] = bcg_list[i]
  return results,measured_bcg_list,cols

print('I DEFINED THE FUNCTIONS')
print('IMMA START COMPUTING NOW -- PLS WAIT ')

#actually computing all the measurements -- (not totally sure whats going on here)
print(os.cpu_count())
poolnr = os.cpu_count()
with Pool(poolnr) as pp:
	all_res1, all_res2 = pp.map(do_all, [np.arange(int(nbcg/2)),np.arange(int(nbcg/2),nbcg)]) 

#printing results if u want
# print(all_res1)
# all_res1[0][int(nbcg/2):] = all_res2[0][int(nbcg/2):]
# all_res1[1][int(nbcg/2):] = all_res2[1][int(nbcg/2):]
# print(all_res1[1])

print('FINSIHED COMPUTING! WRITING RESULTS OUT NOW! THANKS FOR WAITING')
#writing results to file
write_it(all_res1[0],all_res1[1],colnames=all_res1[2])