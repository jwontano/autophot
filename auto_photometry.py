# John Montano - UCI @ montano3@uci.edu
# Photometry routines for reverberation mapping 

# Packages used 
import lmfit as lft
import matplotlib.pyplot as plt
import numpy as np
import numpy.ma as ma
import pandas as pd
import photutils as pu
import random as random

from astropy import wcs
from astropy.io import fits
from astropy.stats import sigma_clipped_stats, signal_to_noise_oir_ccd
from astropy.time import Time
from astropy.utils.exceptions import AstropyWarning
from astropy.wcs import WCS
from datetime import datetime, timedelta
from matplotlib import rc
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from PyAstronomy import pyasl
from scipy import stats

import argparse, csv, datetime, matplotlib, os, requests, shutil, sys, time, warnings, pycali

warnings.simplefilter('ignore', category=AstropyWarning)
def phot(directory=None,
    outdir=None,
    coords=None,
    img_filter=None,
    target_name=None,
    snr_cutoff=0., 
    telname=None, 
    rdnoise=0., 
    pixsc=1., 
    gain=1., 
    rad=[5.,15.,20.,12.], 
    sky_method='DAOphot', 
    sky_hist=False, 
    img_hist=False, 
    renew=True, 
    naver_on=False, 
    centroid_show=False,
    MPC=None,
    c_func='com',
    mjd_cut=None):
    
    if (telname == None or gain == 1. or pixsc == 1.) and telname != 'lco':
        print("Warning: Default values have been enabled. Please check that the proper gain and pixelscale have been inputted, default values pixelscale = 1.0, gain = 1.0, read_noise = 1.0.")
        # set up response to ask if user wishes to continue
        print("Do you wish to continue knowing this? y or n")
        user_response = input()
        if user_response == 'n':
            raise SystemExit()
    # creating the autophot directories if they don't exist
    if os.path.exists(outdir+'/autophot') == False:
        os.mkdir(outdir+'/autophot')
    if img_filter != None:
        if os.path.exists(outdir+'/autophot/'+str(img_filter)) == False:
            os.mkdir(outdir+'/autophot/'+str(img_filter))
        if os.path.exists(outdir+'/autophot/'+str(img_filter)+'/nightcsv') == False:
            os.mkdir(outdir+'/autophot/'+str(img_filter)+'/nightcsv')

    # given a directory will look for and only return a list of filenames that end with fits - .new are something that came in a set of images I was given, so if the fits ends in something else you can add that
    files = [f for f in os.listdir(directory) if os.path.isfile(directory+'/'+f) and (f.endswith(".fits") or f.endswith(".new") or f.endswith(".fit"))]
    # Useful condition if you don't want to rerun all the images in the directory you're feeding - it will check if some data files already exist and skip them thereby only doing new images 
    if renew == False:
        files = [f for f in files if os.path.exists(outdir+'/autophot/'+str(img_filter)+'/nightcsv/phot_'+f+'.csv') == False]
    # big loop that chews through the fits files - does the photometry and saves the results 
    for f,file in enumerate(files):
        # Open your fits file and seperate the data and header.
        hdul = fits.open(directory+"/"+file)
        data = np.asarray(hdul[0].data)
        hdr = hdul[0].header
        # specific for Las Cumbres Observatory since the headers always have these labels
        if telname == 'lco':
            siteid = hdr['SITEID']
            telid = hdr['TELID'] 
            encid = hdr['ENCID']
            tn = siteid+'-'+encid+'-'+telid
            rdnoise = hdr['RDNOISE']
            gain = hdr['GAIN']
            pixsc = hdr['PIXSCALE']
        else:
            tn = telname
        # load up some header information
        exp = hdr['EXPTIME']
        # for images that are averaged, specific for Zowada data 
        if naver_on:
            naver = hdr['NAVER']
        # convert from DN to e-/s
        data = data * gain / exp

        # OPTION
        # Quick look at the histogram of the entire image
        if img_hist == True:
            fig, ax = plt.subplot()
            ax.hist(data.flatten(),bins=5000,log=True)
            ax.set_title('Image Pixel Histogram')
            ax.set_xlabel('Pixel Value')
            ax.set_ylabel('e-/s')
            plt.show()
        # time stamp
        DATE = hdr['DATE-OBS']
        t = Time(DATE, scale='utc')
        MJD = t.mjd
        HJD = pyasl.helio_jd(MJD+.5, coords[0][0], coords[0][1]) + 2.4e6
        if mjd_cut != None:
            if mjd_cut >= MJD:
                continue

        hdul.close()
        # adding half the exposure time to the HJD
        HJD = HJD + (exp*.5)/(3600*24)
        MJD = MJD + (exp*.5)/(3600*24)
        # convert aperture sizes from arcsec to pixel values
        radii = np.asarray(rad)/pixsc
        boxsize = int(rad[3]/pixsc)
        # Some telescopes have bad wcs - or use astrometry to create a wcs. Headers ends up being super messy and astropy's wcs package will get confused if any wcs headers exist from other wcs routines. 
        # Best solution would be for the individual telescopes to fix this themselves or create a good routine to fix the wcs headers and/or add new wcs? 
        w = WCS(hdr)
        # converts the RA, DEC array to pixel coords for the aperture later
        # if this fails then the WCS is most likely bad or cannot be read by astropy
        try:
            pixlist = w.all_world2pix(coords,1)
            pixlist[np.isnan(pixlist)] = 0
        except:
            print('Failed WCS')
            continue

        datshape = np.shape(data)
        xmax, ymax = datshape[0], datshape[1]
        # centering function - 2d gaussian or center of mass
        if c_func == '2dg':
            c_func = pu.centroid_2dg
        elif c_func == 'com':
            c_func = pu.centroid_com
        centers = []
        # go through the list of pixel positions given and find the centroid of a box on pixel position. 
        # there is a way to input just the array of positions, but if one fails then the whole thing stops
        # so I just do pair by pair
        for pix in pixlist:
            try:
                # if the converted pixel positions aren't "real"
                if pix[0] < 0 or pix[1] < 0 or np.isnan(pix[0]):
                    centers.append((xmax/2.,ymax/2.))
                else:
                    sx,sy= pu.centroid_sources(data,xpos=pix[0],ypos=pix[1],box_size=boxsize,centroid_func=c_func)
                    # checking if the centroid is good with simple conditions, doesn't check anything to do with the box
                    if sx[0] < 0 or sx[0] > xmax or sy[0] < 0 or sy[0] > ymax or np.isnan(sx[0]) or np.isnan(sy[0]):
                        centers.append((xmax/2.,ymax/2.))
                        print('1st centering failed.')
                    else:
                        ssx,ssy = pu.centroid_sources(data,xpos=sx[0],ypos=sy[0],box_size=int(boxsize/2),centroid_func=c_func)
                        if np.isnan(ssx[0]) or np.isnan(ssy[0]):
                            centers.append(xmax/2.,ymax/2.)
                            print('2nd centering failed.')
                        else:
                            centers.append((ssx[0],ssy[0]))
                            # print('Starting: %5f, %5f' % (pix[0],pix[1]))
                            # print('First Center: %5f , %5f' % (sx[0],sy[0]))
                            # print('Second Center: %5f , %5f\n' % (ssx[0],ssy[0]))
            except:
                centers.append((xmax/2.,ymax/2.))
                print('Pixel position error')

        # cleaner version of the above chunk however if there is a single bad case it doesn't work at all
        # px, py = np.transpose(pixlist)
        # sx, sy = pu.centroid_sources(data,xpos=px,ypos=py,box_size=boxsize,centroid_func=pu.centroid_com)
        # ssx, ssy = pu.centroid_sources(data,xpos=sx,ypos=sy,box_size=int(boxsize/2),centroid_func=pu.centroid_com)


        # Aperture photometry. Straight forward - run an aperture for the object and an annulus ring to estimate the background 
        # checks = pu.CircularAperture(pixlist, r=radii[0])
        aperture = pu.CircularAperture(centers, r=radii[0])
        annulus_aperture = pu.CircularAnnulus(centers, r_in=radii[1], r_out=radii[2])
        box = pu.RectangularAperture(pixlist,w=radii[3],h=radii[3])
        box2 = pu.RectangularAperture(pixlist,w=radii[3]/2,h=radii[3]/2)

        # little sanity check to see view the centroiding on the images
        if centroid_show == True:
            fig, ax = plt.subplots(figsize=(15,15))
            # hard to get the projection to work on all images uniformly :(
            # if you know you have a safe WCS uncomment the line below and the plots will have a WCS projection
            # plt.subplot(projection=w)
            # can change the inputs of vmin and vmax this is just simple check
            mean, std = np.nanmean(data),np.nanstd(data)
            lower, upper = mean-abs(std), mean+abs(std)
            ax.imshow(data, origin='lower',cmap='Greys_r',vmax=upper, vmin=lower)
            aperture.plot(color='yellow', lw=2, label='Centroid')
            annulus_aperture.plot(color='red', lw=2)
            # checks.plot(color='green',lw=2, label='Initial Pixel Position')
            box.plot(color='red',lw=2)
            box2.plot(color='red', lw=2)
            ax.set_title(str(file)+' centroiding')
            plt.legend()
            plt.show()

        # This will produce a phot table specific to astropy, a lot like a pandas dataframe
        phot = pu.aperture_photometry(data, aperture)
        # sky background estimation
        # create an annulus mask
        # some chunks taken from Astropy examples
        annulus_masks = annulus_aperture.to_mask(method='center')
        bkg_median, bkg_stdev = [], []
        # Allow for modular sky subtraction method
        for m,mask in enumerate(annulus_masks):
            # Do stats on the annulus ring to find a robust sky background
            # grabbing annulus data and flattening it
            annulus_data = mask.multiply(data)
            annulus_data_1d = annulus_data[mask.data > 0]
            # sigma clip default is 3.0 sigma
            mean, median_sigclip, stdev = sigma_clipped_stats(annulus_data_1d)
            # if true will produce a histogram plot of the annulus ring aperture
            # can do this to check if the background estimate is reasonable
            if sky_hist == True:
                x, bins, __ = plt.hist(annulus_data_1d,bins=25)
                plt.title('Sky Annulus Histogram '+str(file)+' '+str(m))
                plt.vlines(mean,0,x.max(),color='red', label='Mean='+str(np.round(mean,4)))
                plt.vlines(median_sigclip,0,x.max(),color='green', label='Median='+str(np.round(median_sigclip,4)))
                plt.vlines((3*median_sigclip - 2*mean), 0, x.max(), color='black', label='DAOphot='+str(np.round(3*median_sigclip - 2*mean,4)))
                plt.xlabel('pixel value')
                plt.legend(frameon=False)
                plt.show()
            # types of background estimation DAOphot and a clipped median, default is DAOphot other methods could be added easily if wanted
            if sky_method == 'DAOphot':
                if mean >= median_sigclip:
                    bkg_median.append(3*median_sigclip - 2*mean)
                else:
                    bkg_median.append(mean)
            elif sky_method == 'clipped median':
                bkg_median.append(median_sigclip)
            bkg_stdev.append(stdev)
        bkg_median, bkg_stdev = np.array(bkg_median), np.array(bkg_stdev)
        # Like pandas dataframes you can create new columns. Perform that for various info we want to save
        phot['annulus_median'] = bkg_median
        # background standard deviation
        phot['bkg_stdev'] = bkg_stdev
        # aperture background value
        phot['aper_bkg'] = bkg_median * aperture.area
        # Can do array operations 
        # aperture sum subtracted the background
        phot['aper_sum_bkgsub'] = phot['aperture_sum'] - phot['aper_bkg']
        # Signal to noise is calculated from the general SNR equation and error is found from going backwards from it.
        SN, error = np.zeros(len(coords)), np.zeros(len(coords))
        # some sanitation
        for i in range(len(coords)):
            if phot['aper_sum_bkgsub'][i] <= 0 or phot['xcenter'][i].value == xmax/2.:
                phot['aper_sum_bkgsub'][i] = 0
                SN[i] = 0
            else:
                # I don't see dark eps ever so it is assumed zero. There is a gain option, but I already did that way in the beginning so we don't need to input it here 
                SN[i] = signal_to_noise_oir_ccd(t=exp, source_eps=phot['aper_sum_bkgsub'][i], sky_eps=phot['annulus_median'][i], dark_eps=0, rd=rdnoise, npix=np.pi*radii[0]**2)
                # can set a signal to noise cutoff value, useful if you have high SNR data
                if SN[i] <= snr_cutoff:
                    phot['aper_sum_bkgsub'][i] = 0 
                    SN[i] = 0
                else:
                    error[i] = 1/SN[i]
        # bunch of tags that are useful to keep track of at this level. Can prune them if needed, but they help track down issues.
        phot['telname'] = tn
        if tn == None:
            phot['MPC'] = tn
        else:
            phot['MPC'] = MPC[tn]
        # specific for zowada data, but can also be made specific to any combined/averaged image by adding naver into a header. naver being the number of images used to average.
        # Error is fixed in these cases so that they are 1/sqrt(naver) along with the photometric error.
        if naver_on:
            phot['error'] = error/np.sqrt(naver)
            phot['naver_error'] = error/np.sqrt(naver)
        else:
            phot['error'] = error
        # signal to noise ratio
        phot['SN'] = SN
        phot['time'] = HJD
        phot['modtime'] = MJD
        # file name
        phot['file'] = file
        # filter band
        phot['filter'] = img_filter
        # exposure time 
        phot['expt'] = exp
        # raw counts calculated from the apeture photometry, for sanity checks
        phot['Raw Counts'] = phot['aper_sum_bkgsub']*exp/gain
        pandaphot = phot.to_pandas()      
        # save data for the single image 
        pandaphot.to_csv(outdir+'/autophot/'+str(img_filter)+'/nightcsv/phot_'+file+'.csv')
        # simple way to track how much its workign through
        print(f+1, ' out of ', len(files), ' done')
    print('Photometry Complete')
    return

def phot_to_one(outdir,img_filter,coords):
    '''
    Go through photometry data tables which produces a single data table of the apertures at given positions and seperate them by telescope, comparison star, AGN etc. 
    Kind of redundant, but having tables of just these values is useful.
    '''
    files = [f for f in os.listdir(outdir+'/autophot/'+str(img_filter)+'/nightcsv') if os.path.isfile(outdir+'/autophot/'+str(img_filter)+'/nightcsv/'+f) and (f.endswith(".csv"))]
    agnpd = pd.DataFrame()
    for file in files:
        pdcsv = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/nightcsv/'+file)
        agnpd = agnpd.append(pdcsv.loc[0])
    agnpd.to_csv(outdir+'/autophot/'+str(img_filter)+'/AGNphots.csv')
    if os.path.exists(outdir+'/autophot/'+str(img_filter)+'/telphot') == False:
        os.mkdir(outdir+'/autophot/'+str(img_filter)+'/telphot')
        print('Creating telphot directory in '+str(img_filter)+' directory')

    for j in range(1,len(coords)):
        starpd = pd.DataFrame()  
        for file in files:
            pdcsv = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/nightcsv/'+file)
            starpd = starpd.append(pdcsv.loc[j])
        starpd.to_csv(outdir+'/autophot/'+str(img_filter)+'/star_photometry_'+str(j)+'.csv')  
        tel_star = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/star_photometry_'+str(j)+'.csv')
        tel_list = tel_star["MPC"]
        telnames = set(tel_list)
        telnames = list(telnames)
        for name in telnames:
            name_star = tel_star.loc[tel_star["MPC"] == name]
            name_star.to_csv(outdir+'/autophot/'+str(img_filter)+'/telphot/'+name+'_star_phot'+str(j)+'.csv')
            print(name+' '+str(j)+' star filed ')
    agnpd = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/AGNphots.csv')
    for name in telnames:
        agn = agnpd.loc[agnpd["MPC"] == name]
        agn.to_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_AGNphots.csv')
        print(name+' AGN csv created')
    return

def scale_AGN(outdir,img_filter,target_name,who=[],clean=False,clean_cutoff=3):
    # takes the nightly scaling factors and applies it to the AGN photometry
    # also has an option to clean large error points
    tel_pd = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/star_photometry_1.csv')
    tel_list = tel_pd["MPC"]
    telnames = set(tel_list)
    for name in telnames:
        if (name in who) == False and len(who) != 0:
            continue
        else:
            agn =  pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/'+str(name)+'_AGNphots.csv')
            flux, error, time, modtime, telname, MPC, files = agn["aper_sum_bkgsub"].to_numpy(), agn["error"].to_numpy(), agn["time"].to_numpy(), agn["modtime"].to_numpy(), agn["telname"].to_numpy(), agn["MPC"].to_numpy(), agn["file"]
            s = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/'+str(name)+'_star_chisq_scale.csv')
            scale = s["scale"].to_numpy()
            scaled_flux = flux * scale
            scaled_flux[scaled_flux <= 0] = np.nan
            scaled_error = flux * error * scale
            serror = error * scale
            serror[serror <= 0] = np.nan
            scaled_error[scaled_error <= 0] = np.nan
            median = np.nanmedian(scaled_flux)
            # this applies a sigma clipping on the UNCERTAINTIES, so this will get rid of data points with unreasonably large errors
            if clean:
                mean, median, std = sigma_clipped_stats(scaled_error)
                sigmas = (scaled_error - median) / std
                for p in range(len(scaled_flux)):
                    if sigmas[p] > clean_cutoff: 
                        scaled_error[p], scaled_flux[p] = np.nan, np.nan
                    else:
                        continue
            agnpd = pd.DataFrame()
            # save the scaled AGN data
            agnpd["scaled AGN flux"], agnpd["scaled AGN error"], agnpd["telname"], agnpd["MPC"], agnpd["time"], agnpd["modtime"], agnpd["file"] = scaled_flux, scaled_error, telname, MPC, time, modtime, files
            agnpd.to_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_scaledAGNphot.csv')
            print(name+' scaled AGN csv created')
    # another cleaning for sanity
    clean_lc(outdir,img_filter)
    return

def plot_scaled(outdir,img_filter,target_name='NNF',cdict={},who=[],filename=False):
    # plot the scaled AGN light curve
    tel_pd = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/star_photometry_1.csv')
    tel_list = tel_pd["MPC"]
    telnames = set(tel_list)
    for name in telnames:
        if (name in who) == False and len(who) != 0:
            continue
        else:
            fig, ax = plt.subplots(figsize=(10,10))
            agn = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_scaledAGNphot.csv')
            scaled_flux = agn['scaled AGN flux']
            scaled_error = agn['scaled AGN error']
            time = agn['time']
            # define a color dictionary before hand for telescopes
            # if not defined it will generate a random color
            if len(cdict) != 0:
                try:
                    color = cdict[name]
                except:
                    r = random.random() 
                    b = random.random() 
                    g = random.random() 
                    color = (r,g,b) 
            else:
                r = random.random() 
                b = random.random() 
                g = random.random() 
                color = (r,g,b) 
            ax.errorbar(time,scaled_flux,yerr=scaled_error,fmt='.',ecolor=color,color=color,ms=2,elinewidth=.45)
            # will print the filename for each associated point, in interactive mode this can be useful to track down bad images
            if filename == True:
                for i, txt in enumerate(agn['file']):
                    ax.annotate(txt, (time[i],scaled_flux[i]),fontsize=8,rotation=45.,fontstretch=.01,alpha=.2)
            ax.set_title(str(target_name)+' '+name+' '+str(img_filter),fontsize=24)
            ax.set_xlabel('HJD',fontsize=16)
            ax.set_ylabel('Flux',fontsize=16)
            ax.minorticks_on()
            ax.tick_params(axis='y',direction='in',which='both')
            ax.tick_params(axis='x',direction='in',which='both')
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            ax.tick_params(which='both',width=1.2)
            ax.tick_params(which='major',length=6,labelsize=12)
            ax.tick_params(which='minor',length=3)
            ax.yaxis.set_ticks_position('both')
            ax.xaxis.set_ticks_position('both')
            # uncomment if you want timestamps
            # s = datetime.datetime.now().strftime('%c')
            # plt.figtext(0.005,.98,s=s)
            # can add the datetime timestamp to the saved jpg file name too if wanted
            plt.savefig(outdir+'/autophot/'+str(img_filter)+'/'+name+'_'+img_filter+'_plot_'+str(img_filter)+'.jpg',dpi=300)
        plt.show()
    return

def plot_align(outdir,img_filter,target_name,reference=None,cdict={},all_black=False):
    # for before intercalibration - plot all the light curves to a reference light curve
    # simple alignment with just scaling the means
    tel_pd = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/star_photometry_1.csv')
    tel_list = tel_pd["MPC"]
    telnames = set(tel_list)
    telnames = list(telnames)
    if reference != None and reference in telnames:
        telnames.remove(reference)
    else:
        print("choosing random reference")
        # set() is kinda random
        reference = telnames[0]
        telnames.remove(telnames[0])
    r_pd = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/'+reference+'_scaledAGNphot.csv')
    ref_flux = r_pd["scaled AGN flux"].to_numpy()
    ref_error = r_pd["scaled AGN error"].to_numpy()
    ref_time = r_pd["time"].to_numpy()
    ref_mean = np.nanmean(ref_flux)
    fig, ax = plt.subplots(figsize=(10,10))
    ms, elw = 2, .45
    if all_black == True:
        color = 'black'
    elif len(cdict) != 0:
        try:
            color = cdict[reference]
        except:
            r = random.random() 
            b = random.random() 
            g = random.random() 
            color = (r,g,b) 
    else:
        r = random.random() 
        b = random.random() 
        g = random.random() 
        color = (r,g,b) 
    ax.errorbar(ref_time,ref_flux,yerr=ref_error,fmt='.',color=color,label=reference,ecolor=color, ms=ms, elinewidth=elw)
    # plot other telescopes light curves
    for name in telnames:
        lcdf = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_scaledAGNphot.csv')
        flux = lcdf["scaled AGN flux"].to_numpy()
        error = lcdf["scaled AGN error"].to_numpy()
        time = lcdf["time"].to_numpy()
        lc_mean = np.nanmean(flux)
        lc_factor = ref_mean/lc_mean
        fflux = flux*lc_factor
        ferror = error*lc_factor
        if all_black == True:
            color = 'black'
        elif len(cdict) != 0:
            try:
                color = cdict[name]
            except:
                r = random.random() 
                b = random.random() 
                g = random.random() 
                color = (r,g,b) 
        else:
            r = random.random() 
            b = random.random() 
            g = random.random() 
            color = (r,g,b) 
        ax.errorbar(time,fflux,yerr=ferror,fmt='.',color=color,label=name,ecolor=color, ms=ms, elinewidth=elw)
    ax.set_ylabel("Scaled Counts")
    ax.set_xlabel("HJD")
    ax.set_title(str(target_name)+' '+str(img_filter)+' Simple Align')
    ax.minorticks_on()
    ax.tick_params(axis='y',direction='in',which='both')
    ax.tick_params(axis='x',direction='in',which='both')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='both',width=1.2)
    ax.tick_params(which='major',length=6,labelsize=12)
    ax.tick_params(which='minor',length=3)
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.legend()
    plt.savefig(outdir+'/autophot/'+img_filter+'/'+target_name+'_'+img_filter+'_simple_align.jpg')
    plt.show()
    return

def clean_lc(outdir,img_filter):
    '''
    Cleans up a light curve if they have any bad values somehow, such as negatives, 0's, nans
    '''
    tel_pd = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/star_photometry_1.csv')
    tel_list = tel_pd["MPC"]
    telnames = set(tel_list)
    for name in telnames:
        lc = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_scaledAGNphot.csv')
        lc = lc.replace(np.nan,0)
        lc = lc.drop(lc[lc["scaled AGN flux"]  <= 0].index)
        lc = lc.sort_values(by=['time'],ascending=True)  
        lc.to_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_scaledAGNphot.csv')
    return 


def condense(outdir,img_filter,who=[],delta=0.25):
    '''
    Time averaging routine from Aaron Barth, converted from IDL to python, default is a delta of .25 days, will take a single light curve and condense it.
    This will take the real flux-value light curves and condense them, some line changes can make it to condense anything else
    '''
    tel_pd = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/star_photometry_1.csv')
    tel_list = tel_pd["MPC"]
    telnames = list(set(tel_list))
    for name in telnames:
        if (name in who) == False and len(who) != 0:
            continue
        lc = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_'+img_filter+'_flamb.csv')
        lc = lc.drop(lc[lc["flux"]  <= 0].index)
        lc = lc.sort_values(by=['time'],ascending=True)
        flux = lc['flux'].to_numpy()
        error = lc['error'].to_numpy()
        time = lc['time'].to_numpy()
        telname = lc['telname']
        mpc = lc['MPC']
        cflux = []
        cerror = []
        ctime = []
        cmodtime = []
        time_points = np.size(time)
        mtime = lc['modtime'].to_numpy()
        i=0
        while i < time_points:
            ntime = time[i]
            w = np.where(abs(time - ntime) < delta)
            wsize = np.size(w)
            wtop = np.max(w)+1
            if wsize == 1:
                cflux.append(flux[i])
                cerror.append(error[i])
                ctime.append(time[i])
                cmodtime.append(mtime[i])
            else:
                # load the slices
                sflux = np.copy(flux[i:wtop])
                serror = np.copy(error[i:wtop])
                stime = np.copy(time[i:wtop])
                smodtime = np.copy(mtime[i:wtop])
                weights = 1.0 / serror**2
                nflux = np.sum(weights * sflux) / np.sum(weights)
                cflux.append( nflux )
                cerror.append( np.sqrt( 1.0/ np.sum(weights)))
                ctime.append(np.mean(stime))
                cmodtime.append(np.mean(smodtime))
            i = i + wsize
        cflux = np.asarray(cflux)
        condensed_AGN = pd.DataFrame(cflux,columns=['flux'])
        condensed_AGN['error'] = np.asarray(cerror)
        condensed_AGN['time'] = np.asarray(ctime)
        condensed_AGN['modtime'] = np.asarray(cmodtime)
        condensed_AGN['MPC'] = mpc
        condensed_AGN['telname'] = telname
        condensed_AGN = condensed_AGN.sort_values(by=['time'],ascending=True)
        condensed_AGN.to_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_avg_lc.csv')
    return


def plot_condense(outdir,img_filter,target_name,cdict={},mpc={},who=[],multi=True,single=False,counts=False):
    # single filter-band plot
    tel_pd = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/star_photometry_1.csv')
    tel_list = tel_pd["MPC"]
    telnames = list(set(tel_list))
    fig, ax = plt.subplots(figsize=(10,10))
    for name in telnames:
        if (name in who) == False and len(who) != 0:
            continue
        if len(cdict) != 0:
            try:
                color = cdict[name]
            except:
                r = random.random() 
                b = random.random() 
                g = random.random() 
                color = (r,g,b)     
        else:
            r = random.random() 
            b = random.random() 
            g = random.random() 
            color = (r,g,b) 
        clc = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_avg_lc.csv')
        time = clc["time"]
        flux = clc["flux"]
        error = clc["error"]
        ax.errorbar(time,flux,yerr=error,fmt='.',color=color,label=name,ecolor=color,ms=2,elinewidth=.45)
    ax.set_xlabel("HJD")
    ax.set_ylabel('f$_\lambda$ ($10^{-14}$ erg cm$^{-2}$ s$^{-1}$ $\mathrm{\AA}^{-1}$)')
    ax.set_title(str(target_name)+" "+str(img_filter))
    ax.minorticks_on()
    ax.tick_params(axis='y',direction='in',which='both')
    ax.tick_params(axis='x',direction='in',which='both')
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which='both',width=1.2)
    ax.tick_params(which='major',length=6,labelsize=12)
    ax.tick_params(which='minor',length=3)
    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.legend()
    plt.savefig(outdir+'/autophot/'+str(img_filter)+'/'+str(target_name)+'_'+img_filter+'_condense_multitel_plot.jpg',dpi=300,bbox_inches='tight')
    plt.show()

    if single == True:
        for name in telnames:
            if (name in who) == False and len(who) != 0:
                continue
            if len(cdict) != 0:
                color = cdict[name]
            else:
                r = random.random()

                b = random.random()
                g = random.random()
                color = (r,g,b)
            fig, ax = plt.subplots(figsize=(10,10))
            if counts:
                clc = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_scaledAGNphot.csv')
            else:
                clc = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_avg_lc.csv')
            time = clc["time"]
            flux = clc["flux"]
            error = clc["error"]
            ax.errorbar(time,flux,yerr=error,fmt='.',color=color,label=name,ecolor=color,ms=2,elinewidth=.45)
            ax.set_xlabel("HJD")
            ax.set_ylabel('f$_\lambda$ ($10^{-14}$ erg cm$^{-2}$ s$^{-1}$ $\mathrm{\AA}^{-1}$)')
            ax.set_title(str(name)+" "+str(target_name)+" "+str(img_filter))
            ax.minorticks_on()
            ax.tick_params(axis='y',direction='in',which='both')
            ax.tick_params(axis='x',direction='in',which='both')
            ax.xaxis.set_minor_locator(AutoMinorLocator())
            ax.yaxis.set_minor_locator(AutoMinorLocator())
            ax.tick_params(which='both',width=1.2)
            ax.tick_params(which='major',length=6,labelsize=12)
            ax.tick_params(which='minor',length=3)
            ax.yaxis.set_ticks_position('both')
            ax.xaxis.set_ticks_position('both')
            ax.legend()
            plt.savefig(outdir+'/autophot/'+str(img_filter)+'/'+str(target_name)+'_'+str(name)+'_'+img_filter+'_condense_multitel_plot.jpg',dpi=300,bbox_inches='tight')
        plt.show()

# excess variance measurement that was used to do uncertainty expansion before intercalibration
# can still be used to check excess variance for data sets in general
def excess_var(outdir,img_filter,starlist,who=[]):
    tel_pd = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/star_photometry_1.csv')
    tel_list = tel_pd["MPC"]
    telnames = set(tel_list)
    for name in telnames:
        if (name in who) == False and len(who) != 0:
            continue
        else:
            star_norm_exvar = []
            star_means = []
            for i in range(1,len(starlist)):
                star = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/telphot/'+name+'_star_phot'+str(i)+'_normed.csv')
                normed_flux = star["scaled flux"].to_numpy()
                mean_flux = np.nanmean(normed_flux)
                star_means.append(mean_flux)
                normed_flux_error = star["scaled flux error"].to_numpy()
                normed_flux_error[normed_flux_error == 0] = np.nan
                normalized_excess_variance = np.nansum(((normed_flux - mean_flux)**2 - normed_flux_error**2)) / (np.size(normed_flux) * mean_flux**2) 
                if normalized_excess_variance < 0.:
                    star_norm_exvar.append(np.nan)
                else:
                    star_norm_exvar.append(normalized_excess_variance)
            sigmas = np.zeros(len(starlist))
            try:
                stdev = np.nanstd(star_norm_exvar)
                sigmas = (star_norm_exvar-np.nanmedian(star_norm_exvar))/stdev
            except:
                sigmas[sigmas == 0] = np.nan
            for sig in range(0,len(sigmas)):
                if abs(sigmas[sig]) >= 2:
                    star_norm_exvar[sig] = np.nan
            if len(set(star_norm_exvar)) == 1:
                avg_ex_var = 0.
            else:
                avg_ex_var = np.nanmean(star_norm_exvar)
            agn = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_scaledAGNphot.csv')
            flux = agn["scaled AGN flux"].to_numpy()
            agn_mean = np.nanmean(flux)
            old_error = np.nan_to_num(agn["scaled AGN error"].to_numpy())
            time, telname, modtime, mpc, files = agn["time"].to_numpy(), agn["telname"], agn["modtime"], agn["MPC"], agn["file"]
            exvar_error = avg_ex_var*(agn_mean**2)
            new_error = np.sqrt( old_error**2 + exvar_error )
            name_AGN = pd.DataFrame()
            name_AGN["scaled AGN flux"], name_AGN["scaled AGN error"], name_AGN["telname"], name_AGN["MPC"], name_AGN["time"], name_AGN["modtime"], name_AGN["file"] = flux, new_error, telname, mpc, time, modtime, files
            name_AGN.to_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_scaledAGNphot.csv')

def cali_format(outdir,flt,target_name,reference=None,condense=True,pycali_datadir=None,raw=False,mjd=False):
    # takes the lightcurves and formats them to be used by pycali
    # user chooses the reference light curve
    for f in flt:
        tel_pd = pd.read_csv(outdir+'/autophot/'+str(f)+'/star_photometry_1.csv')
        tel_list = tel_pd["MPC"]
        names = set(tel_list)
        names = list(names)
        if reference != None and (reference in names) == True:
            r = names.index(reference)
            names[r], names[0] = names[0], names[r]
        cali_av_lc = pd.DataFrame(columns=['time', 'flux', 'error'])
        for n in names:
            if condense == True:
                tf_lc = pd.read_csv(outdir+'/autophot/'+str(f)+'/'+n+'_avg_lc.csv')
            elif raw == True:
                # if this option is chosen then normalize the light curve flux
                tf_lc = pd.read_csv(outdir+'/autophot/'+str(f)+'/'+n+'_scaledAGNphot.csv')
            else:
                tf_lc = pd.read_csv(outdir+'/autophot/'+str(f)+'/'+n+'_'+f+'_flamb.csv')
            if mjd == True:
                tf_time, tf_flux, tf_error, tf_mpc, tf_telnames = tf_lc['modtime'].to_numpy(), tf_lc['flux'].to_numpy()*10**(14), tf_lc['error'].to_numpy()*10**(14), tf_lc['MPC'], tf_lc['telname']
            else:
                tf_time, tf_flux, tf_error, tf_mpc, tf_telnames = tf_lc['time'].to_numpy(), tf_lc['flux'].to_numpy()*10**(14), tf_lc['error'].to_numpy()*10**(14), tf_lc['MPC'], tf_lc['telname']
            enc_time, enc_flux, enc_error = [], [], []
            for y in range(0,len(tf_time)):
                enc_time.append(tf_time[y])
                enc_flux.append(tf_flux[y])
                enc_error.append(tf_error[y])
            to_append = pd.DataFrame({'time':enc_time, 'flux':enc_flux, 'error':enc_error})
            to_append = to_append.replace(np.nan,0)
            to_append = to_append.drop(to_append[to_append['flux']  <= 0.].index)
            to_append['flux'] = to_append['flux'].map('{:,.8e}'.format)
            to_append['error'] = to_append['error'].map('{:,.8e}'.format)
            to_append = to_append.sort_values(by=['time'],ascending=True)
            length = int(len(to_append['flux']))
            try:
                string = '# '+tf_mpc[0]
            except:
                continue
            head = pd.DataFrame({'time': [string], 'flux': [length], 'error': [np.nan]})
            precali = head.append(to_append)
            cali_av_lc = cali_av_lc.append(precali)
        # save light curves to a txt file and format
        if condense == True:
            cali_av_lc.to_csv(outdir+'/autophot/'+target_name+'_'+f+'_cont_avg.txt',sep=' ',index=False,float_format='%.6f')
            with open(outdir+'/autophot/'+target_name+'_'+f+'_cont_avg.txt', 'r') as infile:
                temp = infile.read().replace("\"","")
                temp = temp.replace("time flux error\n","")
            with open(outdir+'/autophot/'+target_name+'_'+f+'_cont_avg.txt','w') as outfile:
                outfile.write(temp)
            print(f+' Averaged Continuum cali file created')
        else:
            cali_av_lc.to_csv(outdir+'/autophot/'+target_name+'_'+f+'_cont.txt',sep=' ',index=False,float_format='%.6f')

            with open(outdir+'/autophot/'+target_name+'_'+f+'_cont.txt', 'r') as infile:
                temp = infile.read().replace("\"","")
                temp = temp.replace("time flux error\n","")
            with open(outdir+'/autophot/'+target_name+'_'+f+'_cont.txt','w') as outfile:
                outfile.write(temp)
            print(f+' Averaged Continuum cali file created')
    return 
################################################################################################################################################################
# comparison star functions

def chisqfunc(params,flux,error,length,mean=None):
    # function that will be inputted into the minimizer - lmfit funcs require params to be passed in, independent variable(s) can also be inputted along with a functions args and kwargs. 
    # residual function to be fed to lmfit
    # need to have the input be as clean as possible or else the minimizer freaks out
    x = params
    y = []
    if length == 1:
        for i in range(1,len(params)+1):
            y.append(x['A'+str(i)]*1)
    else:
        for i in range(length,len(params)):
            y.append(x['A'+str(i)]*1)
    y = np.asarray(y)   
    residual = []
    if length == 1:
        temp = []
        for z in range(np.shape(y)[0]):
            # returning zero immediately to avoid running into any issues with lmfit - lmfit tends to use zero and nan values as markers for ending the fit 
                # standard chi sq. (factor*stars_nights_flux - stars_meanflux)**2 / (factor*stars_nights_error)
            chi = ((y[z]*flux[0][z]-mean)**2)/((y[z]*flux[0][z]*error[0][z])**2)
            temp.append(chi)
        residual.append(temp)
    else: 
        for i in range(0,length):
            temp = []
            # one iteration of this loop should return the chi sq array of a single star, the entire nested loop will return a large array of chi sq values for each star and each night
            for z in range(len(flux[i])):
                # returning zero immediately to avoid running into any issues with lmfit - lmfit tends to use zero and nan values as markers for ending the fit 
                if y[z] == 0 or error[i][z] == 0 or x['star'+str(i)] == 0:
                    chi = 0
                else:
                    # standard chi sq. (factor*stars_nights_flux - stars_meanflux)**2 / (factor*stars_nights_error)
                    chi = ((y[z]*flux[i][z]-x['star'+str(i)])**2)/((y[z]*flux[i][z]*error[i][z])**2)
                temp.append(chi)
            residual.append(temp)
    residual = np.asarray(residual)
    # print(np.nansum(residual)/len(params))    
    # print(len(params))
    return residual 

def fit(outdir,coords,img_filter,who=[]):
    # takes the photometry of all the comparison stars to "flatten"
    # assume that the star has a non-variable flux so the mean is what the fluxes are being fitted to.
    # using all given stars do a chi^2 minimization for the scale factors per night
    tel_star = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/star_photometry_1.csv')
    tel_list = tel_star["MPC"]
    names = set(tel_list)
    names = list(names)
    for name in names:
        if (name in who) == False and len(who) != 0:
            continue
        else:
            all_flux, all_error, all_mean, all_scale = [], [], [], []
            for p in range(1,len(coords)):
                star = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/telphot/'+name+'_star_phot'+str(p)+'.csv')
                flux, err = star["aper_sum_bkgsub"].to_numpy(), star["error"].to_numpy()
                # last second sanitation
                flux[flux <= 0] = np.nan
                if len(flux) == 0 or np.sum(np.isnan(flux)) == np.size(flux):
                    flux_median = 0.
                else:
                    flux_median =  np.nanmedian(flux)
                # first guess at scale factors
                scale = flux_median / flux
                # if the initial guess is huge assume it is a bad measurement (super low count below the median)
                for i in range(len(scale)):
                    if scale[i] > 10:
                        scale[i] = 0
                flux_error = flux * err
                all_flux.append(flux)
                all_error.append(flux_error)
                all_mean.append(flux_median)
                all_scale.append(scale)
            avg_scale = []
            all_flux, all_error, all_mean, all_scale = np.asarray(all_flux), np.asarray(all_error), np.asarray(all_mean), np.asarray(all_scale)
            n_stars, n_pts = np.shape(all_flux)
            # print(sigma_clipped_stats(all_scale.flatten()),np.nanstd(all_scale.flatten()),np.nanmedian(all_scale.flatten()))

            # fig, ax = plt.subplots()
            # ax.hist(all_scale.flatten(),bins=1000)
            # plt.show()
            # if there exists only one point
            if n_pts == 1:
                bf = pd.DataFrame([1.],columns=['scale'])
                bf['labels'] = ['A1']
                bf.to_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_star_chisq_scale.csv')
                print(name+' stars chi sq min fitted')
                continue
            else:
            # should add more restrictions on scaling factors
                # There is probably a way to transpose the arrays to make this work in much fewer lines of code
                for l in range(n_pts):
                    star_scale = np.zeros(n_stars)
                    for k in range(n_stars):
                        star_scale[k] = (all_scale[k][l])
                    # checking if all stars have real values
                    if np.sum(np.isnan(star_scale)) == (n_stars):
                        avg_scale.append(0.)
                    # elif np.nanstd(star_scale) > np.nanmean(star_scale)*.5:
                    #     avg_scale.append(0.)
                    else:
                        # print(star_scale, np.nanmean(star_scale),np.nanmedian(star_scale), np.nanstd(star_scale))
                        avg_scale.append(np.nanmedian(star_scale))
                # fig, ax = plt.subplots()
                # ax.hist(avg_scale,bins=1000)
                # plt.show()

                params = np.concatenate((all_mean,avg_scale),axis=None)
                params = np.nan_to_num(params)
                print('Setting fit parameters')
                all_flux, all_error = np.nan_to_num(all_flux), np.nan_to_num(all_error)
                pars = lft.Parameters()
                for v in range(np.size(params)):
                    if v < n_stars and n_stars != 1:
                        if params[v] == 0:
                            pars.add(name='star'+str(v),value=params[v],vary=False)
                        else:
                            pars.add(name='star'+str(v),value=params[v],min=(params[v]/2.),max=(params[v]*2.))
                    elif n_stars == 1 and v < 1:
                        continue
                    elif params[v] == 0.:
                        pars.add(name='A'+str(v),value=params[v],vary=False)
                    else:
                        # removed upper limit
                        pars.add(name='A'+str(v),value=params[v],min=0)
                # print(pars)
                print('Fitting chi-squared minimization for '+name)
                if n_stars == 1:
                    mini = lft.minimize(fcn=chisqfunc, params=pars, args=(all_flux,all_error,n_stars,all_mean[0]),nan_policy='omit')
                else:
                    mini = lft.minimize(fcn=chisqfunc, params=pars,args=(all_flux,all_error,n_stars),nan_policy='omit')#,ftol=1e-10,xtol=1e-10)
                place = []
                labels = []
                for usls,param in mini.params.items():
                    if 'star' in usls:
                        continue
                    else:
                        place.append(param.value)
                        labels.append(usls)
                print(lft.fit_report(mini,show_correl=False))

                bf = pd.DataFrame(place,columns=["scale"])
                bf['labels'] = labels
                bf.to_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_star_chisq_scale.csv')
                print(name+' stars chi sq min fitted')

def star_norm(outdir,coords,img_filter,who=[],clean=False):
    # apply the nights scale to stars to check for flatness and create files that have the normalized values of the stars
    # then we can use to get for excess variance, error ratios etc. 
    tel_pd = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/star_photometry_1.csv')
    tel_list = tel_pd["MPC"]
    telnames = list(set(tel_list))
    for name in telnames:
        if (name in who) == False and len(who) != 0:
            continue
        else:
            if os.path.exists(outdir+'/autophot/'+str(img_filter)+'/'+name+'_star_chisq_scale.csv') == False:
                print(name+' has no scale factor file - skipping')
                continue
            star_scale = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_star_chisq_scale.csv')
            applied_scale = star_scale["scale"].to_numpy()
            for z in range(1,len(coords)):
                cstar = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/telphot/'+name+'_star_phot'+str(z)+'.csv')
                flux, telname, names, SN = cstar["aper_sum_bkgsub"].to_numpy(), cstar["telname"].to_numpy(), cstar["file"].to_numpy(), cstar["SN"].to_numpy()
                error = cstar["error"].to_numpy()
                error[error <= 0] = np.nan
                scaled_flux = flux * applied_scale
                scaled_flux[scaled_flux <= 0] = np.nan          
                scaled_error = error * applied_scale * flux
                scaled_error[scaled_error <= 0] = np.nan
                time = cstar["time"].to_numpy()
                if clean == True:
                    emean, emedian, estdev = sigma_clipped_stats(scaled_error)
                    fmean, fmedian, fstdev = sigma_clipped_stats(scaled_flux)
                    for p in range(len(scaled_flux)):
                        if np.isnan(scaled_flux[p]):
                            continue
                        elif (np.abs(scaled_error[p]-emedian)/estdev > 5) or (np.abs(scaled_flux[p]-fmedian)/fstdev > 3):
                            scaled_flux[p], scaled_error[p] = np.nan, np.nan
                        else:
                            continue
                scaled_star = pd.DataFrame(scaled_flux,columns=["scaled flux"])
                scaled_star["scaled flux error"], scaled_star["file"], scaled_star["SN"], scaled_star["telname"], scaled_star["time"] = scaled_error, names, SN, telname, time
                scaled_star.to_csv(outdir+'/autophot/'+str(img_filter)+'/telphot/'+name+'_star_phot'+str(z)+'_normed.csv')
            print(name+' comparison star lightcurves normalized')

def plot_normed(outdir,img_filter,coords,cdict={},mpc={},who=[]):
    tel_pd = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/star_photometry_1.csv')
    telnames = tel_pd["MPC"]
    telnames = set(telnames)
    for name in telnames:
        if (name in who) == False and len(who) != 0:
            continue

        fig, ax = plt.subplots(figsize=(10,10),ncols=2,nrows=len(coords)-1,sharex=True)
        plt.subplots_adjust(hspace=0,wspace=0) 
        for i in range(len(coords)-1):
            old = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/telphot/'+name+'_star_phot'+str(i+1)+'.csv')
            oflux, otime = old['aper_sum_bkgsub'].to_numpy(), old['time']
            oerror = old['error'].to_numpy()*oflux
            oflux[oflux <= 0] = np.nan
            star =  pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/telphot/'+name+'_star_phot'+str(i+1)+'_normed.csv')
            flux, error, time = star["scaled flux"].to_numpy(), star["scaled flux error"].to_numpy(), star["time"]
            if len(cdict) != 0:
                try:
                    color = cdict[name]
                except:
                    r = random.random() 
                    b = random.random() 
                    g = random.random() 
                    color = (r,g,b)
            else:
                r = random.random() 
                b = random.random() 
                g = random.random() 
                color = (r,g,b)
            ms, elw = 2, .45
            if len(coords) == 2:
                ax[1].axhline(y=np.nanmedian(oflux),linestyle='-.',color=color,alpha=.5)
                ax[1].errorbar(otime,oflux,yerr=oerror,fmt='.',color=color,ecolor=color,ms=ms,elinewidth=elw)
                ax[0].axhline(y=np.nanmedian(flux),linestyle='-.',color=color,alpha=.5)
                ax[0].errorbar(time,flux,yerr=error,fmt='.',color=color, ecolor=color,ms=ms,elinewidth=elw)
                ax[0].set_ylabel("Flux",fontsize=16)
                ax[1].yaxis.set_label_position("right")
                ax[1].yaxis.tick_right()
                ax[0].minorticks_on()
                ax[1].minorticks_on()
                ax[0].tick_params(axis='y',direction='in',which='both')
                ax[0].tick_params(axis='x',direction='in',which='both')
                ax[1].tick_params(axis='y',direction='in',which='both')
                ax[1].tick_params(axis='x',direction='in',which='both')
                ax[1].yaxis.set_ticks_position('both')
                ax[1].xaxis.set_ticks_position('both')
                ax[0].yaxis.set_ticks_position('both')
                ax[0].xaxis.set_ticks_position('both')
            else:
                ax[i,0].axhline(y=np.nanmedian(oflux),linestyle='-.',color=color,alpha=.5)
                ax[i,0].errorbar(otime,oflux,yerr=oerror,fmt='.',color=color,ecolor=color,ms=ms,elinewidth=elw)
                ax[i,1].axhline(y=np.nanmedian(flux),linestyle='-.',color=color,alpha=.5)
                ax[i,1].errorbar(time,flux,yerr=error,fmt='.',color=color, ecolor=color,ms=ms,elinewidth=elw)
                ax[i,0].set_ylabel("Flux",fontsize=16)
                ax[i,1].yaxis.set_label_position("right")
                ax[i,1].yaxis.tick_right()
                ax[i,0].minorticks_on()
                ax[i,1].minorticks_on() 
                ax[i,0].tick_params(axis='y',direction='in',which='both')
                ax[i,0].tick_params(axis='x',direction='in',which='both')
                ax[i,1].tick_params(axis='y',direction='in',which='both')
                ax[i,1].tick_params(axis='x',direction='in',which='both')
                ax[i,1].yaxis.set_ticks_position('both')
                ax[i,1].xaxis.set_ticks_position('both')
                ax[i,0].yaxis.set_ticks_position('both')
                ax[i,0].xaxis.set_ticks_position('both')
        if len(coords) == 2:
            ax[0].set_xlabel("HJD",fontsize=16)
            ax[1].set_xlabel("HJD",fontsize=16)
            ax[0].set_title('Raw',fontsize=20)
            ax[1].set_title('Normalized',fontsize=20)
        else:
            ax[-1,0].set_xlabel("HJD",fontsize=16)
            ax[-1,1].set_xlabel("HJD",fontsize=16)
            ax[0,0].set_title('Raw',fontsize=20)
            ax[0,1].set_title('Normalized',fontsize=20)
        s=datetime.datetime.now().strftime('%c')
        plt.figtext(0.005,.98,s=s)
        fig.suptitle(name+' Comparison Stars '+img_filter)
        plt.savefig(outdir+'/autophot/'+str(img_filter)+'/'+name+'_'+img_filter+'_stars_normed.jpg',dpi=300,bbox_inches='tight')
    plt.show()




def true_mag_retrieve(outdir,tar=None, mags=None, who=[],feedback=False,f_to=[]):
    # Takes true magnitude values for specific comparison stars used in various filters used to put photometry measurements on a true flux scale (cgs units). 
    # Standard flux to magnitude routine
    # zp = zero point in flux 2.5 * log_10(wavelength zeropoint)
    # fp = flux point of comparison star 10^(-(M-zp))/2.5
    # input wavelength in angstroms
    # make function to input zero points for filters?
    if mags != None:
        tags, B, V, u, g, r, ip, z = mags[0], mags[1], mags[2], mags[3], mags[4], mags[5], mags[6], mags[7] 
    # this shoud be removed at some point
    if tar == 'Mrk817': 
        tags, B, V, u, g, r, ip, z = [1,2,3,4,5],[16.349,16.218,16.274,13.558,15.29],[15.504,15.019,15.264,12.519,13.978],[17.00,17.56,17.10,15.18,16.51],[15.86,15.58,15.72,13.68,14.85],[15.49,14.86,15.25,12.70,14.01],[15.37,14.62,15.09,12.81,13.80],[15.35,14.51,15.04,13.31,13.61]

    # list of zeropoints for BVugriz filters (for STORM2)
    # should add option to input different zeropoints
    zp = [2.5*np.log10(632*10**(-11)),2.5*np.log10(363.1*10**(-11)),2.5*np.log10(859.5*10**(-11)),2.5*np.log10(466.9*10**(-11)),2.5*np.log10(278.0*10**(-11)),2.5*np.log10(185.2*10**(-11)),2.5*np.log10(131.5*10**(-11))]
    zp, B, V, u, g, r, ip, z = np.asarray(zp),np.asarray(B),np.asarray(V),np.asarray(u),np.asarray(g),np.asarray(r),np.asarray(ip),np.asarray(z)
    # using the zero points convert the flux 
    fp = np.asarray([10**(-(B-zp[0])/2.5),10**(-(V-zp[1])/2.5),10**(-(u-zp[2])/2.5),10**(-(g-zp[3])/2.5),10**(-(r-zp[4])/2.5),10**(-(ip-zp[5])/2.5),10**(-(z-zp[6])/2.5)])
    if len(f_to) != 0:
        img_filters = f_to
        temp_fp = []
        # maybe a dict would be useful here
        for filt in f_to:
            if filt == 'B':
                temp_fp.append(fp[0])
            elif filt == 'V':
                temp_fp.append(fp[1])
            elif filt == 'u':
                temp_fp.append(fp[2])
            elif filt == 'g':
                temp_fp.append(fp[3])
            elif filt == 'r':
                temp_fp.append(fp[4])
            elif filt == 'i':
                temp_fp.append(fp[5])
            elif filt == 'z':
                temp_fp.append(fp[6])
        fp = np.asarray(temp_fp)
    # in the case that no filters are given assume the standard BVugriz filters (?) 
    # this should be changed in the future
    else:
        img_filters = ['B','V','u','g','r','i','z']

    # Finding the real flux scale factor by taking the comparison star fluxes and converting to real flux then taking the average of those factors
    # while using a sigma clip to get rid of outlier factors (e.g. one stars scale factor is off compared to the rest that are somewhat uniform with each other)
    for h,img_filter in enumerate(img_filters):
        tel_pd = pd.read_csv(outdir+'autophot/'+str(img_filter)+'/star_photometry_1.csv')
        telnames = set(tel_pd["MPC"])
        for name in telnames:
            if (name in who) == False and len(who) != 0:
                continue
            else:
                mf = []
                for i in tags:
                    s = pd.read_csv(outdir+'autophot/'+str(img_filter)+'/telphot/'+str(name)+'_star_phot'+str(i)+'_normed.csv')
                    mf.append(np.nanmean(s["scaled flux"].to_numpy()))
                mf = np.asarray(mf)
                mf[mf == 0] = np.nan
                factor = fp[h]/mf
                stdev = np.nanstd(factor)
                mean = np.nanmean(factor)
                for k in range(len(tags)):
                    if np.abs(factor[k]-mean) > stdev:
                        factor[k] = np.nan
                mean = np.nanmean(factor)
                if feedback == True:
                    plt.title(name+' '+img_filter)
                    plt.scatter(tags,factor)
                    plt.hlines(mean,1,5,label='mean')
                    plt.hlines(mean+stdev,1,5,label='1 $\sigma$',linestyles='dashed')
                    plt.hlines(mean-stdev,1,5,linestyles='dashed')
                    plt.ylabel('Scale Factor')
                    plt.legend()
                    plt.show()
                agn = pd.read_csv(outdir+'autophot/'+str(img_filter)+'/'+str(name)+'_scaledAGNphot.csv')
                agnflux = agn["scaled AGN flux"].to_numpy()
                true_agn = agnflux*mean
                agnerror = agn["scaled AGN error"].to_numpy()
                time, modtime, mpc, telname = agn["time"].to_numpy(), agn["modtime"].to_numpy(), agn["MPC"], agn["telname"]
                true_error = agnerror*mean
                true_agn[true_agn == 0] = np.nan
                # save the agn f_lambda into a csv 
                tagn = pd.DataFrame(true_agn,columns=["flux"])
                tagn["error"], tagn["time"], tagn["modtime"], tagn["MPC"], tagn["telname"], tagn["file"] = true_error, time, modtime, mpc, telname, agn["file"]
                tagn = tagn.replace(np.nan,0)
                tagn = tagn.drop(tagn[tagn["flux"]  <= 0].index)
                dropping = np.where(tagn["error"]/tagn["flux"]*100 > 10.0)
                tagn = tagn.drop(index = dropping[0])
                tagn = tagn.sort_values(by=['time'],ascending=True)
                tagn.to_csv(outdir+'/autophot/'+str(img_filter)+'/'+str(name)+'_'+img_filter+'_flamb.csv')
    return

def plot_cgs(cdict={},whoo=None,outdir=None,who=[],f_to=[],multi=True,single=False,all_black=False,filenames=False,exclude=[]):
    if len(f_to) != 0:
        img_filters = f_to
    else:
        img_filters = ['B','V','u','g','r','i','z']
    ndict = {'u': 'u$^\prime$' , 'B': 'B', 'V': 'V', 'g': 'g$^\prime$', 'r': 'r$^\prime$', 'i': 'i$^\prime$', 'z': 'z$_s$'}
    
    for img_filter in img_filters:
        tel_pd = pd.read_csv(outdir+'/autophot/'+img_filter+'/star_photometry_1.csv')
        names = set(tel_pd["MPC"])
        fig, ax = plt.subplots(figsize=(10,10))
        ax.set_title(whoo+' '+ndict[img_filter],fontsize=24)
        ms, elw = 2, .45
        for name in names:
            if (name in who) == False and len(who) != 0:
                continue
            if (name in exclude) == True and len(exclude) != 0:
                continue
            tagn = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/'+str(name)+'_'+img_filter+'_flamb.csv')
            true_agn = tagn["flux"]
            true_error = tagn["error"]
            time = tagn["time"]
            files = tagn["file"]
            if all_black == True:
                color = 'black'
            elif len(cdict) != 0:
                try:
                    color = cdict[name]
                except:
                    r = random.random() 
                    b = random.random() 
                    g = random.random() 
                    color = (r,g,b) 
            else:
                r = random.random() 
                b = random.random() 
                g = random.random() 
                color = (r,g,b)
            ax.errorbar(time,true_agn,yerr=true_error,fmt='.',color=color, ecolor=color,ms=ms,elinewidth=elw,label=name)
            if filenames:
                for i, txt in enumerate(files):
                    ax.annotate(txt, (time[i],true_agn[i]),fontsize=8,rotation=45.,fontstretch=.01,alpha=.1)
        ax.minorticks_on()
        ax.tick_params(axis='y',direction='in',which='both')
        ax.tick_params(axis='x',direction='in',which='both')
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(which='both',width=1.2)
        ax.tick_params(which='major',length=6,labelsize=12)
        ax.tick_params(which='minor',length=3)
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        ax.set_xlabel('HJD',fontsize=16)
        ax.set_ylabel('f$_\lambda$ (ergs cm$^{-2}$ s$^{-1}$ $\AA^{-1}$)',fontsize=16)
        ax.legend(markerscale=2,borderaxespad=0,handletextpad=0,prop={'size': 14})
        s = datetime.datetime.now().strftime('%c')
        ax.yaxis.set_label_coords(-.08,.5)
        plt.figtext(0.001,.98,s=s)
        plt.savefig(outdir+'/autophot/'+str(img_filter)+'/uncal_true_flux_'+whoo+'_'+str(img_filter)+'.pdf', dpi=900,bbox_inches='tight')
        plt.savefig(outdir+'/autophot/'+str(img_filter)+'/uncal_true_flux_'+whoo+'_'+str(img_filter)+'.jpg', dpi=300,bbox_inches='tight')
        plt.show()
        if single == True:
            for name in names:
                if (name in who) == False and len(who) != 0:
                    continue
                tsagn = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/'+str(name)+'_'+img_filter+'_flamb.csv')
                true_agn = tsagn["flux"]
                true_error = tsagn["error"]
                time = tsagn["time"]
                color = 'black'
                fig, ax = plt.subplots(figsize=(10,10))
                ax.set_title(whoo+' '+str(name)+' '+ndict[img_filter],fontsize=24)
                ax.errorbar(time,true_agn,yerr=true_error,fmt='.',color=color, ecolor=color, ms=ms, elinewidth=elw, label=name)
                ax.minorticks_on()
                ax.tick_params(axis='y',direction='in',which='both')
                ax.tick_params(axis='x',direction='in',which='both')
                ax.xaxis.set_minor_locator(AutoMinorLocator())
                ax.yaxis.set_minor_locator(AutoMinorLocator())
                ax.tick_params(which='both',width=1.2)
                ax.tick_params(which='major',length=6,labelsize=12)
                ax.tick_params(which='minor',length=3)
                ax.yaxis.set_ticks_position('both')
                ax.xaxis.set_ticks_position('both')
                ax.set_xlabel('HJD',fontsize=16)
                ax.set_ylabel('f$_\lambda$ (ergs cm$^{-2}$ s$^{-1}$ $\AA^{-1}$)',fontsize=16)
                ax.legend(markerscale=2,borderaxespad=0,handletextpad=0,prop={'size': 14})
                s = datetime.datetime.now().strftime('%c')
                ax.yaxis.set_label_coords(-.08,.5)
                plt.figtext(0.001,.98,s=s)
                plt.savefig(outdir+'/autophot/'+str(img_filter)+'/uncal_true_flux_'+whoo+'_'+str(name)+'_'+str(img_filter)+'.pdf', dpi=900,bbox_inches='tight')
                plt.savefig(outdir+'/autophot/'+str(img_filter)+'/uncal_true_flux_'+whoo+'_'+str(name)+'_'+str(img_filter)+'.jpg', dpi=900,bbox_inches='tight')
            plt.show()

def error_ratio(outdir,coords,img_filter):
    tel_pd = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/star_photometry_1.csv')
    tel_list = tel_pd["MPC"]
    telnames = set(tel_list)
    telnames = list(telnames)
    for name in telnames:
        star_norm_exvar = []
        star_means = []
        for i in range(0,len(coords)-1):
            star = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/telphot/'+name+'_star_phot'+str(i+1)+'_normed.csv')
            normed_flux = star["scaled flux"].to_numpy()
            mean_flux = np.nanmean(normed_flux)
            star_means.append(mean_flux)
            normed_flux_error = star["scaled flux error"].to_numpy()
            normed_flux_error[normed_flux_error == 0] = np.nan
            normalized_excess_variance = np.nansum(((normed_flux - mean_flux)**2 - normed_flux_error**2)) / (np.size(normed_flux) * mean_flux**2) 
            if normalized_excess_variance < 0.:
                star_norm_exvar.append(np.nan)
            else:
                star_norm_exvar.append(normalized_excess_variance)
        sigmas = np.zeros(len(coords)-1)
        try:
            stdev = np.nanstd(star_norm_exvar)
            sigmas = (star_norm_exvar-np.nanmedian(star_norm_exvar))/stdev
        except:
            sigmas[sigmas == 0] = np.nan
        for sig in range(0,len(sigmas)):
            if abs(sigmas[sig]) >= 2:
                star_norm_exvar[sig] = np.nan
        if len(set(star_norm_exvar)) == 1:
            avg_ex_var = 0.
        else:
            avg_ex_var = np.nanmean(star_norm_exvar)
        # print(star_norm_exvar)
        # take average excess variance and add in quadrature to the agn
        agn = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_scaledAGNphot.csv')
        flux = agn["scaled AGN flux"].to_numpy()
        error = agn["scaled AGN error"].to_numpy()
        agn["scaled AGN error"] = np.sqrt(error**2 + (flux*avg_ex_var)**2)
        agn.to_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_scaledAGNphot.csv')

        amax = np.nanmax(flux)
        amin = np.nanmin(flux)
        agn_mean = np.nanmean(flux)
        fig, ax = plt.subplots(figsize=(10,10))
        ax.axvline(x=agn_mean,linestyle='--',color='red',label='$<f_{AGN}>$')
        ax.scatter(star_means,star_norm_exvar,color='black',marker='x',label='$<f_{star}>$')
        ax.set_title('Fractional Excess Scatter for '+str(name))
        ax.set_ylabel('Fractional Excess Scatter $\sigma_{nx}$')
        ax.set_xlabel('Flux')
        ax.fill_betweenx(y=[np.nanmin(star_norm_exvar)*.98,np.nanmax(star_norm_exvar)*1.02],x1=[amin,amin],x2=[amax,amax],color='red',alpha=0.3,label='Flux Range')
        plt.legend(frameon=False)
        plt.savefig(outdir+'/autophot/'+str(img_filter)+'/'+str(name)+'_excess_scatter.jpg',dpi=600,bbox_inches='tight')
    plt.show()

########################################################################################################################################################################
# pycali plots
# ask where to find the cali files or give the cali files?
def plot_multicali(tn,cdict,avg=False,cali_file):
    rc('font', **{'family':'Dejavu Sans'})
    rc('text', usetex=True)
    filters = ['u','B','g','V','r','i','z']
    fig, ax = plt.subplots(nrows=len(filters),ncols=1,sharex=True,figsize=(15,10))
    plt.subplots_adjust(hspace=0)
    ax[0].set_title(tn,fontsize=28)
    for i,f in enumerate(filters):

        fn = cali_file+'/'+tn+'_'+str(f)+'_cont.txt_cali'
        cali = pd.read_csv(fn,header=None,usecols=[0,1,3,5],delimiter=' ',names=['time','flux','error','names'])
        names = set(cali['names'])
        names = list(names)
        # maybe add plot last 
        try:
            xb = names.index('liverpool')
            names[xb], names[-1] = names[-1], names[xb]
        except:
            pass
        ms, elw, alpha = 1, .7, 1
        for n in names:
            sub = cali.loc[cali['names'] == n]
            flux = sub['flux']
            error = sub['error']
            time = sub['time']

            # if n == 'asiago':
            #     as_v = pd.read_csv('/home/korbinite5646/AGN_home/MRK817/autophot/'+f+'/asiago_avg_lc.csv')
            #     lc_size = np.size(flux)
            #     error = (as_v['error']*10**(14))[0:lc_size]

            if len(cdict) != 0 and (n in cdict) == True:
                color = cdict[n]
            else:
                r = random.random() 
                b = random.random() 
                g = random.random() 
                color = 'black'
            if n == 'V39' or n == 'V37' or n == 'Z24' or n == 'Z31':
                ax[i].errorbar(time,flux,yerr=error,fmt='.',ms=ms,ecolor=color,elinewidth=elw,color=color,label='LCO-'+n,alpha=alpha)
            elif n == 'wise':
                ax[i].errorbar(time,flux,yerr=error,fmt='.',ms=ms,ecolor=color,elinewidth=elw,color=color,label='WISE',alpha=alpha)
            elif n == 'zowada':
                ax[i].errorbar(time,flux,yerr=error,fmt='.',ms=ms,ecolor=color,elinewidth=elw,color=color,label='ZOW',alpha=alpha)
            elif n == 'liverpool':
                ax[i].errorbar(time,flux,yerr=error,fmt='.',ms=ms,ecolor=color,elinewidth=elw,color=color,label='LT',alpha=alpha)
            elif n == 'lijiang':
                ax[i].errorbar(time,flux,yerr=error,fmt='.',ms=ms,ecolor=color,elinewidth=elw,color=color,label='LJ',alpha=alpha)
            elif n == 'ratir':
                ax[i].errorbar(time,flux,yerr=error,fmt='.',ms=ms,ecolor=color,elinewidth=elw,color=color,label='SPM',alpha=alpha)
            elif n == 'F65':
                ax[i].errorbar(time,flux,yerr=error,fmt='.',ms=ms,ecolor=color,elinewidth=elw,color=color,label='FTN-M3',alpha=alpha)
            else:
                ax[i].errorbar(time,flux,yerr=error,fmt='.',ms=ms,ecolor=color,elinewidth=elw,color=color,label=n)
            ax[i].text(0.02,0.1,s=f,transform=ax[i].transAxes,fontsize= 16)
        ax[i].minorticks_on()
        ax[i].tick_params(axis='y',direction='in',which='both')
        ax[i].tick_params(axis='x',direction='in',which='both')
        ax[i].xaxis.set_minor_locator(AutoMinorLocator())
        ax[i].yaxis.set_minor_locator(AutoMinorLocator())
        ax[i].tick_params(which='both',width=1.2)
        ax[i].tick_params(which='major',length=6,labelsize=12)
        ax[i].tick_params(which='minor',length=3)
        ax[i].yaxis.set_ticks_position('both')
        ax[i].xaxis.set_ticks_position('both')
        if i == 0:
            handles, labels = ax[0].get_legend_handles_labels()
            handles = [h[0] for h in handles]
        else:
            temph, templ = ax[i].get_legend_handles_labels()
            temph = [t[0] for t in temph]
            for p in range(0,len(templ)):
                if templ[p] not in labels:
                    handles.append(temph[p])
                    labels.append(templ[p]) 
    # add ability to change ncol and where to place legend
    fig.legend(handles,labels,markerscale=1.,borderaxespad=0,handletextpad=0.0,prop={'size': 10},loc='upper left',bbox_to_anchor=(.12,.88),bbox_transform=plt.gcf().transFigure, ncol=10, frameon=False)
    ax[len(filters)-1].set_xlabel('HJD',fontsize=20)
    ax[3].set_ylabel( 'f$_\lambda$ ($10^{-14}$ erg cm$^{-2}$ s$^{-1}$ $\mathrm{\AA}^{-1}$)',fontsize=20)
    ax[3].yaxis.set_label_coords(-.05,.5)
    # s = datetime.datetime.now().strftime('%c')
    # plt.figtext(0.005,.98,s=s)
    # plt.savefig(wheretosave+tn+'cali_avg.jpg',dpi=300)
    plt.show()
    
def plot_cali(tn,cdict,calidir,filters=[]):
    rc('font', **{'family':'Dejavu Sans'})
    rc('text', usetex=True)
    ms = 2
    elw = .45
    for i,f in enumerate(filters):
        cali = pd.read_csv(calidir,header=None,usecols=[0,1,3,5],delimiter=' ',names=['time','flux','error','names'])
        names = set(cali['names'])
        fig, ax = plt.subplots(figsize=(20,20))
        ax.set_title(tn+' '+f,fontsize=24)
        alpha = 1
        for n in names:
            sub = cali.loc[cali['names'] == n]
            flux = sub['flux']
            error = sub['error']
            time = sub['time']
            if len(cdict) != 0 and (n in cdict) == True:
                color = cdict[n]
            else:
                r = random.random() 
                b = random.random() 
                g = random.random() 
                color = (r,g,b)
            ax.errorbar(time,flux,yerr=error,fmt='.',ms=ms,ecolor=color,elinewidth=elw,color=color,label=n)
        ax.minorticks_on()
        ax.tick_params(axis='y',direction='in',which='both')
        ax.tick_params(axis='x',direction='in',which='both')
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
        ax.tick_params(which='both',width=1.2)
        ax.tick_params(which='major',length=6,labelsize=12)
        ax.tick_params(which='minor',length=3)
        ax.yaxis.set_ticks_position('both')
        ax.xaxis.set_ticks_position('both')
        handles, labels = ax.get_legend_handles_labels()
        handles = [h[0] for h in handles]
        fig.legend(handles,labels,borderaxespad=0,handletextpad=0.0,prop={'size': 12},loc='center',bbox_to_anchor=(.5,.862),bbox_transform=plt.gcf().transFigure, framealpha=0.8,markerscale=2,ncol=len(names),frameon=False)
        
        ax.set_xlabel('HJD',fontsize=16)
        ax.set_ylabel( 'f$_\lambda$ ($10^{-14}$ erg cm$^{-2}$ s$^{-1}$ $\mathrm{\AA}^{-1}$)',fontsize=16)
        ax.yaxis.set_label_coords(-.07,.5)
        # s = datetime.datetime.now().strftime('%c')
        # plt.figtext(0.005,.98,s=s)
        # plt.savefig(tn+'_'+f+'_cali_avg.pdf',dpi=300)
    plt.show()


'''
LCO API program to see which proposals I am a part of and I have options to download images to my station. The code is from the LCO API documentation and github with some changes for myself. 
see https://github.com/LCOGT/
'''
def get_token(plist=False):
    response = requests.post(
        'https://archive-api.lco.global/api-token-auth/',
        data = {
            'username': 'yourusername',
            'password': 'yourpassword'
        }
    ).json()
    token = response.get('token')
    headers = {'Authorization': 'Token ' + token}
    if plist:
        profile_response = requests.get('https://archive-api.lco.global/profile/',headers=headers)
        try: 
            profile_response.raise_for_status()
        except requests.exceptions.HTTPError as exc:
            print('Request failed: {}'.format(profile_response.content))
            raise exc
        profile_dict = profile_response.json()
        print('Username: {} proposals'.format(profile_dict['username']))
        print('You are a member of {} proposals'.format(len(profile_dict['profile']['proposals'])))
        for p in profile_dict['profile']['proposals']:
            print(p)
    return token

def download(start_date,end_date,proposal_id,data_directory,target=None,public=False):
    token = get_token()
    headers = {'Authorization': 'Token ' + token}
    if target == None:
        archive_response = requests.get('https://archive-api.lco.global/frames/?reduction_level=91&limit=50&proposal_id='+proposal_id+'&'+'start='+start_date+'&'+'end='+end_date+' 23%3A59&&',headers=headers).json()
    elif public:
        archive_response = requests.get('https://archive-api.lco.global/frames/?reduction_level=91&limit=100&'+'proposal_id='+proposal_id+'&'+'start='+start_date+'&'+'end='+end_date+' 23%3A59&'+'target_name='+target+'&public=true',headers=headers).json()
    else:
        archive_response = requests.get('https://archive-api.lco.global/frames/?reduction_level=91&limit=20&'+'proposal_id='+proposal_id+'&'+'start='+start_date+'&'+'end='+end_date+' 23%3A59&'+'target_name='+target+'&public=false',headers=headers).json()
    frames = archive_response['results']
    files = []
    for i,j,k in os.walk('/home/korbinite5646/AGN_home/FITS/'):
        if k != []:
            files = files + k
        else:
            continue
    while True:
        for frame in frames:
            if (str(frame['filename'][:-3]) in files) or (str(frame['filename']) in files):
                print('{} file exists already skipping...'.format(frame['filename']))
                continue
            else:
                with open(os.path.join(data_directory, frame['filename']), 'wb') as f:
                    print('Downloading {}...\n'.format(frame['filename']))
                    f.write(requests.get(frame['url']).content)
        if archive_response.get('next'):
            archive_response = requests.get(archive_response['next'], headers=headers).json()
            frames = archive_response['results']
        else:
            break

def lcoarchive_download(start_date,end_date,propID,target,prop_list=False,public=False,datadir=None):


    if propID == None and prop_list == False:
        print('Please input a proposal identification to choose from. If you would like to see the proposal list use the -list option.')
        exit()
    if prop_list:
        token = get_token(True)
        exit()
    # Where to download the files into
    # need to add the size of and number of images that'll be added 
    # and after that a confirmation to download
    data_d = '/home/korbinite5646/AGN_home/FITS/API/'

    if end_date is None and start_date is None:
        end_date =  datetime.today()
        start_date = end_date - timedelta(days=7)
        end_date = end_date.strftime("%Y-%m-%d")
        start_date = start_date.strftime("%Y-%m-%d")
        print('Searching last 7 days from {} to {}...\n'.format(start_date,end_date))
    elif end_date is None:
        end_date = datetime.today()
        end_date = end_date.strftime("%Y-%m-%d")
        print('Searching date range from {} to {}...\n'.format(start_date,end_date))
        
    download(start_date,end_date,propID,data_d,target,public=public)

# will take a folder of fits images and MOVES them to organize into a given directory
def filter_n_sort_images(directory,outdir=None,lco=False):
    # would have to change the f_dict based onthe names of the filters in your own headers
    f_dict = {'gp': 'g', 'g': 'g', 'SDSS-G': 'g', 'ip': 'i', 'i': 'i', 'SDSS-I': 'i', 'rp': 'r', 'r': 'r', 'SDSS-R': 'r', 'R': 'r', 'up': 'u', 'u': 'u', 'SDSS-U': 'u', 'zs': 'z', 'z': 'z', 'SDSS-Z': 'z', 'V': 'V', 'Bessell-V': 'V', 'B': 'B', 'Bessell-B': 'B', "Y": 'Y'}
    # look in the directory given and get a list of the fits files - only does surface directory does NOT go into subdirectories
    files = [f for f in os.listdir(directory) if os.path.isfile(directory+'/'+f) and (f.endswith(".fits") or f.endswith(".new"))]
    for file in files:
        # Las Cumbres Observatory data has a number of useful header tags
        if lco:
            imagedir = directory+'/'+file
            hdul = fits.open(image)
            hdr = hdul[0].header
            target = hdr['OBJECT']
            prop = hdr['PROPID']
            fltr = hdr['FILTER']
            fltr = f_dict[fltr]
            hdul.close()
            if os.path.exists(outdir+prop) == False:
                print('Creating proposal storage directory for '+prop)
                os.mkdir(outdir+prop)
            if os.path.exists(outdir+prop+'/'+target)==False:
                print('Creating directory for '+target+' in '+prop)
                os.mkdir(outdir+prop+'/'+target)
            if os.path.exists(outdir+prop+'/'+target+'/'+fltr) == False:
                print('Creating filter storage directory for '+fltr)
                os.mkdir(outdir+prop+'/'+target+'/'+fltr)
            try:
                shutil.move(imagedir,outdir+prop+'/'+target+'/'+fltr+'/'+file)
                print(file+' moved to '+outdir+prop+'/'+target+'/'+fltr)
            except:
                print('Failed to move file '+f)
                continue
        else:
            imagedir = directory+'/'+ file
            hdul = fits.open(imagedir)
            hdr = hdul[0].header
            # hoping that the images given have a header filter
            try:
                fltr = hdr['FILTER']
            except:
                fltr = hdr['FILTER1']
            # small string dictionary for various names that header filters come in
            fltr = f_dict[fltr]
        try:
            shutil.move(imagedir,outdir+'/'+fltr+'/'+file)
            print(file+' moved to '+outdir+'/'+fltr+'/'+file)
        except:
            print('Failed to move file :-(')
            continue
    return

# YOU NEED PYCALI INSTALLED FOR THIS TO RUN
# ENSURE PYCALI IS PROPERLY INSTALLED AND CAN BE CALLED
def pycali_run(filt='g', run=1, outdir=''):
    # need to add wrapper inputs
    cfg = pycali.Config()
    # fcont at this point should take from the autophot directory
    # need to add directory creation for /cali and /cali/save
    if run == 1:
        fcont = outdir+"/autophot/"+TARGET+"_"+filt+"_cont_avg.txt"
    elif run == 2:
        fcont = outdir+"/autophot/cali/save/"+TARGET+"_"+filt+"_cont_avg.txt"
    cfg.setup(fcont=fcont,
              nmcmc=15000, ptol=0.1,
              scale_range_low=0.5, scale_range_up=2.0,
              shift_range_low=-1.0, shift_range_up=1.0,
              syserr_range_low=0.0, syserr_range_up=0.2,
              errscale_range_low=0.5, errscale_range_up=2.0,
              sigma_range_low=1.0e-4, sigma_range_up=1.0,
              tau_range_low=1.0, tau_range_up=1.0e4,
              fixed_scale=False, fixed_shift=False,
              fixed_syserr=True, fixed_error_scale=True)
    cali = pycali.Cali(cfg)
    cali.mcmc()
    cali.get_best_params()
    cali.output()
    cali.recon()

# for when pycali is run the first time 
def produce_clean(filters,outdir,obj_name,ref):
    for filt in filters:
        # load reconstruction and continuum
        r_fn = outdir+'/autophot/'+obj_name+'_'+filt+'_cont_avg.txt_recon'
        c_fn = outdir+'/autophot/'+obj_name+'_'+filt+'_cont_avg.txt_cali'
        recon_lc = pd.read_csv(r_fn,header=None,usecols=[0,3,5],delimiter=' ',names=['time','flux','error'])
        cali_lc = pd.read_csv(c_fn,header=None,usecols=[0,1,3,5],delimiter=' ',names=['time','flux','error','tel'])
        # calculate residuals by interpolation with the reconstruction
        residual = cali_lc['flux'] - np.interp(cali_lc['time'], recon_lc['time'], recon_lc['flux'])
        cali_lc['residual'] = residual
        std_res = residual/cali_lc['error']
        cali_lc['std_res'] = std_res
        res_std = np.std(cali_lc['residual'])
        cali_lc['std normed res'] = cali_lc['residual']/res_std
        clean_idx = np.where((abs(cali_lc['std normed res']) < 5))
        tel_pd = pd.read_csv(outdir+'/autophot/'+str(filt)+'/star_photometry_1.csv')
        tel_list = tel_pd['MPC']
        names = set(tel_list)
        names = list(names)
        names.remove('zowada')
        names.append('zowada_storm')
        names.append('zowada_x')
        if (ref in names) == True:
                r = names.index(ref)
                names[r], names[0] = names[0], names[r]
        cali_av_lc = pd.DataFrame(columns=['time', 'flux', 'error'])
        for n in names:    
            tf_lc = pd.read_csv(outdir+'/autophot/'+str(filt)+'/'+n+'_avg_lc.csv')
            tf_time, tf_flux, tf_error, tf_mpc, tf_telnames = tf_lc['time'].to_numpy(), tf_lc['flux'].to_numpy()*10**(14), tf_lc['error'].to_numpy()*10**(14), tf_lc['MPC'], tf_lc['telname']
            enc_time, enc_flux, enc_error = [], [], []
            for y in range(len(tf_time)):
                if np.sum(np.round(tf_time[y],6) == np.round(cali_lc['time'][clean_idx[0]],6)) > 0:
                    enc_time.append(tf_time[y])
                    enc_flux.append(tf_flux[y])
                    enc_error.append(tf_error[y])
                else:
                    continue
            length = int(len(enc_time))
            try:
                string = '# '+tf_mpc[0]
            except:
                continue
            head = pd.DataFrame({'time': [string], 'flux': [length], 'error': [np.nan]})
            to_append = pd.DataFrame({'time':enc_time, 'flux':enc_flux, 'error':enc_error})
            to_append = to_append.sort_values(by=['time'],ascending=True)
            precali = head.append(to_append)
            cali_av_lc = cali_av_lc.append(precali)
        cali_av_lc.to_csv(outdir+'/autophot/cali/save/'+TARGET+'_'+filt+'_cont_avg.txt',sep=' ',index=False)
        with open(outdir+'/autophot/cali/save/'+TARGET+'_'+filt+'_cont_avg.txt', 'r') as infile:
            temp = infile.read().replace("\"","")
            temp = temp.replace("time flux error\n","")
        with open(outdir+'/autophot/cali/save/'+TARGET+'_'+filt+'_cont_avg.txt','w') as outfile:
            outfile.write(temp)
        print(filt+' Averaged Continuum cali file created')
    return

# def check_clean(filters):
#   for filt in filters:
#       r_fn, c_fn = '/home/korbinite5646/AGN_home/MRK817/autophot/cali/error/'+filt+'_recon.txt','/home/korbinite5646/AGN_home/MRK817/autophot/cali/error/mrk817_'+filt+'_cont_avg.txt_cali'
#       recon_lc = pd.read_csv(r_fn,header=None,usecols=[0,3,5],delimiter=' ',names=['time','flux','error'])
#       cali_lc = pd.read_csv(c_fn,header=None,usecols=[0,1,3,5],delimiter=' ',names=['time','flux','error','tel'])
#       names = list(set(cali_lc['tel']))
#       residual = cali_lc['flux'] - np.interp(cali_lc['time'], recon_lc['time'], recon_lc['flux'])
#       cali_lc['residual'] = residual
#       std_res = residual/cali_lc['error']
#       cali_lc['std_res'] = std_res
#       res_std = np.std(cali_lc['residual'])
#       cali_lc['std normed res'] = cali_lc['residual']/res_std
#   # full histogram all points
#       fig, ax = plt.subplots(figsize=(15,5))
#       ax.hist(residual, bins=25, color='black')
#       ax.set_title('Mrk 817 Residual Histogram '+filt)
#       ax.set_xlabel('Flux')
#       fig, ax = plt.subplots(figsize=(15,5))
#       for n in names:
#           tel_indx = np.where(cali_lc['tel'] == n)
#           ax.errorbar(cali_lc['time'][tel_indx[0]],residual[tel_indx[0]],yerr=cali_lc['error'][tel_indx[0]],fmt='x',label=n,elinewidth=.7, color=cdict[n])
#       ax.set_title('Mrk 817 Residual over Time '+filt)
#       ax.set_ylabel('Flux')
#       ax.legend()
#   # residual over time plot all tels together
#       fig, ax = plt.subplots(figsize=(15,5))
#       ax.hist(std_res, bins=25, color='gray')
#       ax.set_title('Mrk 817 Residual/Error Histogram '+ filt)
#       ax.set_xlabel('$\sigma$')
#       fig, ax = plt.subplots(figsize=(15,5))
#       for n in names:
#           tel_indx = np.where(cali_lc['tel'] == n)
#           ax.scatter(cali_lc['time'][tel_indx[0]],std_res[tel_indx[0]],marker='x',label=n, color=cdict[n])
#       ax.set_title('Mrk 817 Residual/Error '+filt)
#       ax.set_ylabel('$\sigma$')
#       ax.legend()
#       fig, ax = plt.subplots(figsize=(15,5))
#       ax.hist(cali_lc['std normed res'], bins=25, color='gray')
#       ax.set_title('Mrk 817 Residual/Std Histogram '+filt)
#       ax.set_xlabel('$\sigma$')  
#       fig, ax = plt.subplots(figsize=(15,5))
#       for n in names:
#           tel_indx = np.where(cali_lc['tel'] == n)
#           ax.scatter(cali_lc['time'][tel_indx[0]],cali_lc['std normed res'][tel_indx[0]],marker='x',label=n, color=cdict[n])
#       ax.set_title('Mrk 817 Residual/Std '+filt)
#       ax.set_ylabel('$\sigma$')
#       ax.legend()
#       plt.show()