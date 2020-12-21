import numpy as np
import math
import astropy.io.fits as fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import patches
import skimage
from skimage.measure import *
from skimage.transform import downscale_local_mean
from astropy.modeling import models
import photutils
import time
import statmorph
from astropy.visualization import simple_norm
import sys
sys.path.append('/content/gdrive/My Drive/data/')
import shapes_lib as sl



#======================= Functions ==================================
def testprint(word):
  return print(word)
 
def get_image(loc,downscale=0,expand=1):
    bcg_im2 = fits.getdata(loc)
    if downscale>0:
        # !! Downscaling the image to save runtime !!
        # float is to escape some error
        bcg_im3 = downscale_local_mean(np.array(bcg_im2,dtype=np.float32),
                                        (downscale,downscale))
    else:
        bcg_im3 = bcg_im2
    if expand>0:  
        # Add an edge of zeros for the contour estimation
        bcg_im = np.zeros((bcg_im3.shape[0]+(expand*2),bcg_im3.shape[1]+(expand*2)),dtype=bcg_im3.dtype)
        bcg_im[expand:-1*expand,expand:-1*expand]=bcg_im3
        return bcg_im
    else:
        return bcg_im2
    
    
 
def make_ellipse(epar,colour='k',linestyle='-'):
  # epar: (centroidx,centroidy),minor_axis_length,major_axis_length,orientation
  plot_ellipse = patches.Ellipse(epar[0],epar[1],epar[2],epar[3],
                            fill = False,edgecolor=colour,linestyle=linestyle)
  print('Ellipse params  ',colour, epar[0], epar[1],epar[2],epar[3]) 
  return plot_ellipse
 

def make_ellipse2(epar,colour='k',linestyle='-'):
  # epar: centroidx,centroidy,minor_axis_length,major_axis_length,orientation
  plot_ellipse = patches.Ellipse((epar[0],epar[1]),epar[2],epar[3],epar[4],
                            fill = False,edgecolor=colour,linestyle=linestyle)
  print('Ellipse params  %s  %.1f,%.1f  %.1f %.1f %.1f' %(colour,epar[0], epar[1],epar[2],epar[3],epar[4])) 
  return plot_ellipse
 


def plot_multi_contour(bcg_im, contour_params1=None, contour_params2=None, contour_params3=None,
                        r500=0,r200=0,
                        zoomins=[0,600,800],
                        save=None,show=True):
                        
  """Plot multiple contours on one image"""
  
  if type(bcg_im)==str:
    bcg_im = get_image(loc,downscale=0,expand=1)
  vmin,vmax = np.percentile(bcg_im[bcg_im>0],[5,95])
  
  # Display the (zoomed in) image and plot all contours found
  s=bcg_im.shape
  fig, axes = plt.subplots(1,len(zoomins), figsize=(15,15))
  for ii in range(len(axes)):
    axes[ii].imshow(bcg_im[zoomins[ii]:s[0]-zoomins[ii],zoomins[ii]:s[0]-zoomins[ii]],
                    interpolation='nearest',origin='lower',
                    norm=LogNorm(vmin=vmin,vmax=vmax)) 

  #plotting 1st contour
  bestfit_contour_ellipse = make_ellipse2(contour_params1,colour='g')
  axes[0].add_artist(bestfit_contour_ellipse)

  #second contour
  bestfit_contour_ellipse = make_ellipse2(contour_params2,colour='y')
  axes[0].add_artist(bestfit_contour_ellipse)

  #third contour
  bestfit_contour_ellipse = make_ellipse2(contour_params3,colour='r')
  axes[0].add_artist(bestfit_contour_ellipse)

  # If provided, show ellipse (circle) for R500 and R200 (in pixels!)
  # Ellipse uses diameter, so Rdelta*2
  if r500>0:
    epar=[bcg_im.shape[0]/2,bcg_im.shape[1]/2, 2*r500,2*r500,0.]
    r500_ellipse = make_ellipse2(epar,colour='grey',linestyle='--')
    axes[0].add_artist(r500_ellipse)
  if r200>0:
    epar=[bcg_im.shape[0]/2,bcg_im.shape[1]/2,2*r200,2*r200,0.]  
    r200_ellipse = make_ellipse2(epar,colour='grey',linestyle='--')
    axes[0].add_artist(r200_ellipse)
    
  if not save is None:
      fig.savefig(save,format=save.split('.')[-1])
  if show: 
      fig.show()
  else:
      fig.clf()




def plot_img2(bcg_im,chi_moments_params=None,model_params=None,
                contour_params=None,
                r500=0,r200=0,
                zoomins=[0,600,800],
                save=None,show=True):
  """
  Plot image and contour
  And print contour properties
  """
  if type(bcg_im)==str:
      bcg_im = get_image(loc,downscale=0,expand=1)
  vmin,vmax = np.percentile(bcg_im[bcg_im>0],[5,95])
  #print('Intensities: %.3e  %.3e' %(vmin,vmax), np.percentile(bcg_im,[5,99]))
  

  # Display the (zoomed in) image and plot all contours found
  s=bcg_im.shape
  fig, axes = plt.subplots(1,len(zoomins), figsize=(15,15))
  for ii in range(len(axes)):
    axes[ii].imshow(bcg_im[zoomins[ii]:s[0]-zoomins[ii],zoomins[ii]:s[0]-zoomins[ii]],
                    interpolation='nearest',origin='lower',
                    norm=LogNorm(vmin=vmin,vmax=vmax))  


  # If not provided, find contours
  if contour_params is None:
    contour, contour_props = sl.get_longest_contour_props(bcg_im,np.percentile(bcg_im,60))
    contour_params = [contour_props[0][0].centroid[1], 
                        contour_props[0][0].centroid[0], 
                        contour_props[0][0].minor_axis_length,
                        contour_props[0][0].major_axis_length, 
                       180 - contour_props[0][0].orientation*180/np.pi]
    axes[0].plot(contour[:, 1]-zoomins[0], contour[:, 0]-zoomins[0], linewidth=2)
  bestfit_contour_ellipse = make_ellipse2(contour_params,colour='k')
  axes[0].add_artist(bestfit_contour_ellipse)
  
  # If provided, show ellipse for moments measurements
  if not chi_moments_params is None:
      chi_moments_ellipse = make_ellipse2(chi_moments_params,colour='r')
      axes[0].add_artist(chi_moments_ellipse)
    
  # If provided, show ellipse for Sersic fitting measurements
  if not model_params is None:
      sersic_fit_ellipse = make_ellipse2(model_params,colour='g')
      axes[0].add_artist(sersic_fit_ellipse)
  
  # If provided, show ellipse (circle) for R500 and R200 (in pixels!)
  # Ellipse uses diameter, so Rdelta*2
  if r500>0:
    epar=[bcg_im.shape[0]/2,bcg_im.shape[1]/2, 2*r500,2*r500,0.]
    r500_ellipse = make_ellipse2(epar,colour='grey',linestyle='--')
    axes[0].add_artist(r500_ellipse)
  if r200>0:
    epar=[bcg_im.shape[0]/2,bcg_im.shape[1]/2,2*r200,2*r200,0.]  
    r200_ellipse = make_ellipse2(epar,colour='grey',linestyle='--')
    axes[0].add_artist(r200_ellipse)
    
    
  if not save is None:
      fig.savefig(save,format=save.split('.')[-1])
  if show: 
      fig.show()
  else:
      fig.clf()


def plot_img(bcg_im,chi_moments_params=None,model_morph=None,
                contour=None, contour_props=None,
                r500_05=0,r500=0,r200=0,
                zoomins=[0,600,800],
                save=None,show=True):
  """
  Plot image and contour
  And print contour properties
  """
  if type(bcg_im)==str:
    bcg_im = fits.getdata(bcg_im)
  vmin,vmax = np.percentile(bcg_im[bcg_im>0],[1,99])
  #print('Intensities: %.3e  %.3e' %(vmin,vmax), np.percentile(bcg_im,[5,99]))
  
  
  # If not provided, find contours
  if contour is None:
    contour, contour_props = sl.get_longest_contour_props(bcg_im,np.percentile(bcg_im,60))
  epar = [contour_props[0][0].centroid[::-1], 
            contour_props[0][0].minor_axis_length,
            contour_props[0][0].major_axis_length, 
            180 - contour_props[0][0].orientation*180/np.pi]
  bestfit_contour_ellipse = make_ellipse(epar,colour='k')
  '''          
  bestfit_contour_ellipse = patches.Ellipse(contour_props[0][0].centroid[::-1], 
                            contour_props[0][0].minor_axis_length, 
                            contour_props[0][0].major_axis_length, 
                            180 - contour_props[0][0].orientation*180/np.pi, 
                            fill = False,edgecolor='k')
  print('contour params', contour_props[0][0].centroid[::-1],
                            contour_props[0][0].minor_axis_length, 
                            contour_props[0][0].major_axis_length, 
                            180 - contour_props[0][0].orientation*180/np.pi)
  '''
  
  # If provided, show ellipse for moments measurements
  if not chi_moments_params is None:
    size_a,size_b,pos_angle = sl.axes_lengths_from_chi_moments(chi_moments_params[1],chi_moments_params[2],chi_moments_params[3])
    '''
    chi_moments_ellipse = patches.Ellipse(chi_moments_params[0][::-1],
                            2*size_b, 2*size_a, #diameter, not radius
                            #pos_angle,
                            pos_angle + 90,
                            fill=False,edgecolor='r')
    print('moments params', chi_moments_params[0],
                            size_b, size_a,
                            pos_angle+90)
    '''  
    epar = [chi_moments_params[0][::-1],
            2*size_b, 2*size_a, #diameter, not radius
            pos_angle + 90]
  chi_moments_ellipse = make_ellipse(epar,colour='r')
    
                            
  
  # If provided, show ellipse for Sersic fitting measurements
  if not model_morph is None:
    # Using r=sqrt(a*b)=a*sqrt(q), and ellip=epsilon
    q_sersic= (1-model_morph.sersic_ellip)/(1+model_morph.sersic_ellip)
    '''
    sersic_fit_ellipse = patches.Ellipse((model_morph.sersic_xc,
                                            model_morph.sersic_yc),
                            model_morph.sersic_rhalf*np.sqrt(q_sersic)*2, 
                            model_morph.sersic_rhalf/np.sqrt(q_sersic)*2, # diameter, not radius 
                            math.degrees(model_morph.sersic_theta) + 90,
                            fill=False,edgecolor='g')
    print('Sersic  params', (model_morph.sersic_xc,
                                            model_morph.sersic_yc),
                            model_morph.sersic_rhalf*np.sqrt(q_sersic), 
                            model_morph.sersic_rhalf/np.sqrt(q_sersic), 
                            math.degrees(model_morph.sersic_theta)+90)
    '''
    epar=[(model_morph.sersic_xc,model_morph.sersic_yc),
            model_morph.sersic_rhalf*np.sqrt(q_sersic), 
            model_morph.sersic_rhalf/np.sqrt(q_sersic), 
            (model_morph.sersic_theta*180/np.pi)+90]
  sersic_fit_ellipse = make_ellipse(epar,colour='g')
  
  # If provided, show ellipse (circle) for R500 and R200 (in pixels!)
  if r500_05>0:
    epar=[(bcg_im.shape[0]/2,bcg_im.shape[1]/2),
           2*r500_05,2*r500_05,0.]           # uses diameter, so *2
    r500_05_ellipse = make_ellipse(epar,colour='grey',linestyle='--')
    '''patches.Ellipse((bcg_im.shape[0]/2,
                                    bcg_im.shape[1]/2),
                            2*r500,2*r500,0.,           # uses diameter, so *2
                            fill=False,edgecolor='grey',
                            linestyle='--')
    '''
  
  if r500>0:
    epar=[(bcg_im.shape[0]/2,bcg_im.shape[1]/2),
           2*r500,2*r500,0.]           # uses diameter, so *2
    r500_ellipse = make_ellipse(epar,colour='grey',linestyle='--')
    '''patches.Ellipse((bcg_im.shape[0]/2,
                                    bcg_im.shape[1]/2),
                            2*r500,2*r500,0.,           # uses diameter, so *2
                            fill=False,edgecolor='grey',
                            linestyle='--')
    '''
  if r200>0:
    epar=[(bcg_im.shape[0]/2,bcg_im.shape[1]/2),
           2*r200,2*r200,0.]           # uses diameter, so *2
    r200_ellipse = make_ellipse(epar,colour='grey',linestyle='--')
    '''r200_ellipse = patches.Ellipse((bcg_im.shape[0]/2,
                                    bcg_im.shape[1]/2),
                            2*r200,2*r200,0.,
                            fill=False,edgecolor='grey',
                            linestyle='--')   
    '''                        
  # Display the (zoomed in) image and plot all contours found
  s=bcg_im.shape
  fig, axes = plt.subplots(1,len(zoomins), figsize=(15,15))
  for ii in range(len(axes)):
    axes[ii].imshow(bcg_im[zoomins[ii]:s[0]-zoomins[ii],zoomins[ii]:s[0]-zoomins[ii]],
                    interpolation='nearest',origin='lower',
                    norm=LogNorm(vmin=vmin,vmax=vmax))
  axes[0].plot(contour[:, 1]-zoomins[0], contour[:, 0]-zoomins[0], linewidth=2)
  axes[0].add_artist(bestfit_contour_ellipse)
  if not chi_moments_params is None:
    axes[0].add_artist(chi_moments_ellipse)
  if not model_morph is None:
    axes[0].add_artist(sersic_fit_ellipse)
  if r500_05>0:    
    axes[0].add_artist(r500_05_ellipse)
  if r500>0:    
    axes[0].add_artist(r500_ellipse)
  if r200>0:    
    axes[0].add_artist(r200_ellipse)
    
  if not save is None:
      fig.savefig(save,format=save.split('.')[-1])
  if show: 
      fig.show()
  else:
      fig.clf()





def plot_model(image,model_morph):  
  ny, nx = image.shape
  imy, imx = np.mgrid[0:ny, 0:nx]


  sersic_model = models.Sersic2D(
    amplitude=model_morph.sersic_amplitude,
    r_eff=model_morph.sersic_rhalf,
    n=model_morph.sersic_n,
    x_0=model_morph.sersic_xc,
    y_0=model_morph.sersic_yc,
    ellip=model_morph.sersic_ellip,
    theta=model_morph.sersic_theta)
    
  fitted_model = sersic_model(imx, imy)
  plt.subplot(121)
  plt.imshow(image, cmap='rainbow', origin='lower',
           norm=simple_norm(image, stretch='log', log_a=10000))
  plt.subplot(122)
  plt.imshow(fitted_model, cmap='rainbow', origin='lower',
           norm=simple_norm(image, stretch='log', log_a=10000))
  plt.show()

  plt.imshow(image-fitted_model, cmap='rainbow', origin='lower')
  plt.colorbar()
  plt.show()
  plt.hist((image-fitted_model).flatten(),bins=25)
  plt.show()
  
