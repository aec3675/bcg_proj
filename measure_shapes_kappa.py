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
kappadir = config['files']['kappadir']
Rdelta = config['params']['radius'] #can also be 'rd200'
# w_guess = int(config['params']['w_guess'])

print('KAPPA DIR =', kappadir, type(kappadir))
print('R DELTA =', Rdelta, type(Rdelta))
# print('W GUESS =', w_guess, type(w_guess))




#copying and pasting from ricardos code w some changes
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os, sys
import numpy as np
import pandas as pd
sys.path.append('/content/gdrive/My Drive/data/Ricardo_data')
#!ls '/content/gdrive/My Drive/data/'
import temp_plotting as pl #THIS IS MINE NOT RICARDOS
import shapes_lib as sl
from astropy.io import ascii
from astropy.table import Table
import multiprocessing
from multiprocessing import Pool
import time
# !! Downscaling the image to save runtime !!
from skimage.transform import downscale_local_mean

print('IMPORTED STUFFS')

#######################################################################################################
#NEED TO CHANGE THIS BEFORE RUNNING #
# Just take the clusters where there are both BCG and kappa image
kappa_list = pd.read_csv('/content/drive/My Drive/data/kappa_maps_list.csv', sep=',')
bcg_list = kappa_list['filename_x']
nbcg = len(bcg_list)
print('# of Kappa maps =', nbcg)#,bcg_list)


#defining the function to measure all the diff fits
def measure_all(i):
  # Read in image
  print(kappadir+bcg_list[i])
  bcg_im2 = fits.getdata(kappadir+bcg_list[i])
  # !! Downscaling the image to save runtime !!
  # float is to escape some error
  bcg_im3 = downscale_local_mean(np.array(bcg_im2,dtype=np.float32),(4,4))
  # Add an edge of zeros for the contour estimation
  bcg_im = np.zeros((bcg_im3.shape[0]+2,bcg_im3.shape[1]+2),dtype=bcg_im3.dtype)
  bcg_im[1:-1,1:-1]=bcg_im3

  # Image size in pixels and Mpc
  tsy,tsx = bcg_im.shape
  mpc = tsx/5. # pix/Mpc
  xctr=int(tsx/2); yctr=int(tsy/2)

  # Cluster radii in pixels -- PARAM INPUT R_DELTA
  if Rdelta == 'r500':
    r500=float(r500a[clname=='CL'+bcg_list[i].split('_')[1]][0]) *mpc/1000.
    rdelta = r500
  if Rdelta == 'r200':
    r200=float(r200a[clname=='CL'+bcg_list[i].split('_')[1]][0]) *mpc/1000.
    rdelta = r200
  if Rdelta == 'r500_05':
    r500=float(r500a[clname=='CL'+bcg_list[i].split('_')[1]][0]) *mpc/1000.
    rdelta = r500 / 2
    r500_05 = rdelta

  # Save results
  the_res=[]
  colnames=[]

  # Annulus at Rdelta -- possible param for annulus width??
  annulus_width=0.3*mpc   # Half annulus width actually (in pixels)
  masky,maskx = np.ogrid[-1*yctr:yctr,-1*xctr:xctr]
  mask = maskx**2+masky**2 <= (rdelta+annulus_width)**2
  mask[maskx**2+masky**2 < (rdelta-annulus_width)**2]=False
  annulus_im = bcg_im*1
  annulus_im[~mask]=0.

  t0=time.time()

  # #Compute moments ------------------------------- 
  # # Unweighted moments full image
  # q20,q11,q02 = sl.weighted_higher_order_moments(bcg_im,1.e6,2,xctr,yctr)
  # q00 = sl.weighted_higher_order_moments(bcg_im,1.e6,0,xctr,yctr)[0]
  # Tchi = np.around(np.sqrt((q20+q02)/q00),decimals=2)
  # size_a,size_b,pos_angle = sl.axes_lengths_from_chi_moments(Tchi,(q20-q02)/(q20+q02),2*q11/(q20+q02))
  # # Save results
  # the_res += [xctr,yctr,np.around(np.sqrt((q20+q02)/q00),decimals=2),size_b/size_a,pos_angle+90]
  # colnames += ['moment_unwghtd_xc','moment_unwghtdyc','moment_unwghtd_T','moment_unwghtd_q','moment_unwghtd_pos_ang']


  # # Weighted moments full image
  # q20,q11,q02 = sl.weighted_higher_order_moments(bcg_im,rdelta,2,xctr,yctr)
  # q00 = sl.weighted_higher_order_moments(bcg_im,rdelta,0,xctr,yctr)[0]
  # Tchi = np.around(np.sqrt((q20+q02)/q00),decimals=2)
  # size_a,size_b,pos_angle = sl.axes_lengths_from_chi_moments(Tchi,(q20-q02)/(q20+q02),2*q11/(q20+q02))
  # # Save results
  # the_res += [xctr,yctr,Tchi,size_b/size_a,pos_angle+90]
  # colnames += ['moment_wghtd_'+Rdelta+'_xc','moment_wghtd_'+Rdelta+'_yc','moment_wghtd_'+Rdelta+'_T','moment_wghtd_'+Rdelta+'_q','moment_wghtd_'+Rdelta+'_pos_ang']


  # # Unweighted moments annulus around Rdelta
  # q20,q11,q02 = sl.weighted_higher_order_moments(annulus_im,1.e6,2,xctr,yctr)
  # q00 = sl.weighted_higher_order_moments(annulus_im,1.e6,0,xctr,yctr)[0]
  # print('Q',q20,q11,q02,q00,q20+q02/q00, (q20-q02)/(q20+q02),2*q11/(q20+q02))
  # Tchi = np.around(np.sqrt((q20+q02)/q00),decimals=2)
  # size_a,size_b,pos_angle = sl.axes_lengths_from_chi_moments(Tchi,(q20-q02)/(q20+q02),2*q11/(q20+q02))
  # # Save results
  # the_res += [xctr,yctr,Tchi,size_b/size_a,pos_angle+90]
  # colnames += ['moment_unwghtd_annulus_xc','moment_unwghtd_annulus_yc','moment_unwghtd_annulus_T','moment_unwghtd_annulus_q','moment_unwghtd_annulus_pos_ang']

  # # Plot the annulus measurements
  # mom_p=((xctr,yctr),Tchi, 
  #        np.around((q20-q02)/(q20+q02),decimals=4),
  #        np.around((2*q11)/(q20+q02),decimals=4))

  tm = time.time()


# # Fit Sersic -------------------------------------

#   fit_morph = sl.fit_sersic(bcg_im,mask=~mask,verbose=False,plot_mod=False)
#   q_sersic= (1-fit_morph.sersic_ellip)/(1+fit_morph.sersic_ellip)
#   #q_sersic = -1

#   # Save results
#   if q_sersic==-1:
#     the_res+=[-1,-1,-1,q_sersic,-1]
#   else:
#     the_res+=[fit_morph.sersic_xc,fit_morph.sersic_yc,fit_morph.sersic_rhalf,
#                         q_sersic,fit_morph.sersic_theta*180/np.pi + 90]
#   colnames+=['sersic_fit_xc','sersic_fit_yc','sersic_fit_r','sersic_fit_q','sersic_fit_pos_ang']

  ts = time.time()


  #Measure contour ---------------------------------
  # Measure contour and its properties for contour closest in size to Rdelta
  #at r500
  lcon,lcon_p,mincont = sl.get_contour_at_radius(bcg_im,rdelta,begin=30,end=90)
  a_c,b_c = lcon_p[0][0].major_axis_length,lcon_p[0][0].minor_axis_length
  theta_c = 180 - lcon_p[0][0].orientation*180/np.pi

  # Save results
  if a_c==0: 
    the_res+=[lcon_p[0][0].centroid[0],lcon_p[0][0].centroid[1],np.sqrt(b_c*a_c),-1, -99,]
  else: 
    the_res+=[lcon_p[0][0].centroid[0],lcon_p[0][0].centroid[1],np.sqrt(b_c*a_c),b_c/a_c, theta_c]
  colnames+=['contour_%s_xc' %(Rdelta),'contour_%s_yc' %(Rdelta),'contour_%s_r' %(Rdelta),
             'contour_%s_q' %(Rdelta),'contour_%s_pos_ang' %(Rdelta)]
  #Make visual check images --------------------------
  # Plot and save the image
  pl.plot_img(bcg_im,
              # mom_p,
              # model_morph=fit_morph,
              contour=lcon, contour_props=lcon_p,
              zoomins=[0,0,0],
              r500_05 = (float(r500a[clname=='CL'+bcg_list[i].split('_')[1]][0]) *mpc/1000.)/2,
              r500=float(r500a[clname=='CL'+bcg_list[i].split('_')[1]][0]) *mpc/1000.,
              r200=float(r200a[clname=='CL'+bcg_list[i].split('_')[1]][0]) *mpc/1000.,
              save=args.save_data_dir+'_1.0_kappa_%s.pdf' %(bcg_list[i]),
              show=False)
  tc1 = time.time()

  #at r500*1.1
  lcon,lcon_p,mincont1 = sl.get_contour_at_radius(bcg_im,rdelta*1.1,begin=max(20,mincont-10),end=min(90,mincont+10))
  a_c,b_c = lcon_p[0][0].major_axis_length,lcon_p[0][0].minor_axis_length
  theta_c = 180 - lcon_p[0][0].orientation*180/np.pi

  # Save results
  if a_c==0: 
    the_res+=[lcon_p[0][0].centroid[0],lcon_p[0][0].centroid[1],np.sqrt(b_c*a_c),-1, -99,]
  else: 
    the_res+=[lcon_p[0][0].centroid[0],lcon_p[0][0].centroid[1],np.sqrt(b_c*a_c),b_c/a_c, theta_c]
  colnames+=['contour_1.1_%s_xc' %(Rdelta),'contour_1.1_%s_yc' %(Rdelta),'contour_1.1_%s_r' %(Rdelta),
             'contour_1.1_%s_q' %(Rdelta),'contour_1.1_%s_pos_ang' %(Rdelta)]

  #Make visual check images --------------------------
  # Plot and save the image
  pl.plot_img(bcg_im,
              # mom_p,
              # model_morph=fit_morph,
              contour=lcon, contour_props=lcon_p,
              zoomins=[0,0,0],
              r500_05 = (float(r500a[clname=='CL'+bcg_list[i].split('_')[1]][0]) *mpc/1000.)/2,
              r500=float(r500a[clname=='CL'+bcg_list[i].split('_')[1]][0]) *mpc/1000.,
              r200=float(r200a[clname=='CL'+bcg_list[i].split('_')[1]][0]) *mpc/1000.,
              save=args.save_data_dir+'_1.1_kappa_%s.pdf' %(bcg_list[i]),
              show=False)
  tc2 = time.time()

  #at r500*0.9
  lcon,lcon_p,mincont = sl.get_contour_at_radius(bcg_im,rdelta*0.9,begin=max(20,mincont-10),end=min(90,mincont+10))
  a_c,b_c = lcon_p[0][0].major_axis_length,lcon_p[0][0].minor_axis_length
  theta_c = 180 - lcon_p[0][0].orientation*180/np.pi

  # Save results
  if a_c==0: 
    the_res+=[lcon_p[0][0].centroid[0],lcon_p[0][0].centroid[1],np.sqrt(b_c*a_c),-1, -99,]
  else: 
    the_res+=[lcon_p[0][0].centroid[0],lcon_p[0][0].centroid[1],np.sqrt(b_c*a_c),b_c/a_c, theta_c]
  colnames+=['contour_0.9_%s_xc' %(Rdelta),'contour_0.9_%s_yc' %(Rdelta),'contour_0.9_%s_r' %(Rdelta),
             'contour_0.9_%s_q' %(Rdelta),'contour_0.9_%s_pos_ang' %(Rdelta)]
  
  #Make visual check images --------------------------
  # Plot and save the image
  pl.plot_img(bcg_im,
              # mom_p,
              # model_morph=fit_morph,
              contour=lcon, contour_props=lcon_p,
              zoomins=[0,0,0],
              r500_05 = (float(r500a[clname=='CL'+bcg_list[i].split('_')[1]][0]) *mpc/1000.)/2,
              r500=float(r500a[clname=='CL'+bcg_list[i].split('_')[1]][0]) *mpc/1000.,
              r200=float(r200a[clname=='CL'+bcg_list[i].split('_')[1]][0]) *mpc/1000.,
              save=args.save_data_dir+'_0.9_kappa_%s.pdf' %(bcg_list[i]),
              show=False)
  tc3 = time.time()


  # #Make visual check images --------------------------
  # # Plot and save the image
  # pl.plot_img(bcg_im,
  #             # mom_p,
  #             # model_morph=fit_morph,
  #             contour=lcon, contour_props=lcon_p,
  #             zoomins=[0,0,0],
  #             r500=float(r500a[clname=='CL'+bcg_list[i].split('_')[1]][0]) *mpc/1000.,
  #             r200=float(r200a[clname=='CL'+bcg_list[i].split('_')[1]][0]) *mpc/1000.,
  #             save=args.save_data_dir+'TEST_kappa_%s.pdf' %(bcg_list[i]),
  #             show=False)
  # tp = time.time()

  # End
  print(i, bcg_list[i], rdelta, mincont)
  # print(i, bcg_list[i], 'time',tm-t0,ts-tm,tc-ts,tp-tc,time.time()-t0)
  print(i, bcg_list[i], 'time',tc1-t0,tc2-tc1,tc3-tc2,time.time()-t0)
  return the_res,colnames


#defining the function to write out results
def write_it(results,measured_bcg_list,colnames=None):
  if colnames is None:
    colnames=['moment_xc','moment_yc','moment_T','moment_q','moment_pos_ang',
    		'contour_xc','contour_yc','contour_r','contour_q','contour_pos_ang',
            'sersic_fit_xc','sersic_fit_yc','sersic_fit_r','sersic_fit_q','sersic_fit_pos_ang']
  data = Table({cn: results[:,ii] for ii,cn in enumerate(colnames)})
  data['cluster']=measured_bcg_list
  clnr=[bb.split('_')[1]+'.'+bb.split('_')[3][0] for bb in measured_bcg_list]
  data['clusternr']=clnr
  ascii.write(data,args.save_data_dir+str(Rdelta)+'_contours_kappa_shape_measurements.csv',overwrite=True,format='csv')

#defining funciton that does measure_all for all indices
def do_all(bcg_indices):
  results=None
  measured_bcg_list=np.array(nbcg*['O'*len(bcg_list[0])],dtype=str)
  for i in bcg_indices:
    res,cols = measure_all(i)
    if results is None:
      results=np.zeros((nbcg,len(res)))
    results[i] = res
    measured_bcg_list[i] = bcg_list[i]
  print('MEASURED BCG LIST:', measured_bcg_list)
  return results,measured_bcg_list,cols

print('DEFINED ALL THE FUNCTIONS')

# Cluster properties
clname,r500a,r200a = np.loadtxt('./SPA300/halo_properties/radii_masses_z0.2198.txt',
                                unpack=True,dtype=str,usecols=(0,2,3))

print('IMMA COMPUTE STUFF -- PLS WAIT')

#nbcg=4
#print(os.cpu_count(),nbcg)
poolnr = os.cpu_count()
tstart=time.time()
with Pool(poolnr) as pp:
	all_res1, all_res2 = pp.map(do_all, [np.arange(int(nbcg/2)),np.arange(int(nbcg/2),nbcg)])


all_res1[0][int(nbcg/2):] = all_res2[0][int(nbcg/2):]
all_res1[1][int(nbcg/2):] = all_res2[1][int(nbcg/2):]

print('ALL DONE MEASURING!!')

write_it(all_res1[0],all_res1[1],colnames=all_res1[2])

print('\n total time:',time.time()-tstart)