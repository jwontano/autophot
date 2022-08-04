

from astropy.io import fits
import photutils as pu
import numpy as np
import matplotlib.pyplot as plt
from astropy import wcs
from astropy.stats import sigma_clipped_stats, signal_to_noise_oir_ccd
import pandas as pd
import numpy.ma as ma
from astropy.time import Time
from PyAstronomy import pyasl
from astropy.wcs import WCS
import random as random
from matplotlib import rc
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from astropy.utils.exceptions import AstropyWarning
import datetime, os, star_norm, shutil, warnings
# rc('font', **{'family': 'serif', 'serif': ['Computer Modern Roman'], 'size': 14})
# rc('text', usetex=True)

warnings.simplefilter('ignore', category=AstropyWarning)

def phot(directory,
    outdir,
    coords,
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
    c_func='com'):
    
    # creating the autophot directories
    # FIX THIS TO BE BETTER/EFFICIENT
    if os.path.exists(outdir+'/autophot') == False:
        os.mkdir(outdir+'/autophot')
    if img_filter != None:
        if os.path.exists(outdir+'/autophot/'+str(img_filter)) == False:
            os.mkdir(outdir+'/autophot/'+str(img_filter))
        if os.path.exists(outdir+'/autophot/'+str(img_filter)+'/nightcsv') == False:
            os.mkdir(outdir+'/autophot/'+str(img_filter)+'/nightcsv')

    # given a directory will look for and only return a list of filenames that end with fits - .new are something that came in a set of images I was given, so if the fits ends in something else you can add that
    files = [f for f in os.listdir(directory) if os.path.isfile(directory+'/'+f) and (f.endswith(".fits") or f.endswith(".new") or f.endswith(".fit"))]
    # Useful condition if you don't want to rerun all the images in the directory you're feeding - it will check if some data files already exist and skip them  thereby only doing new images 
    if renew == False:
        files = [f for f in files if os.path.exists(outdir+'/autophot/'+str(img_filter)+'/nightcsv/phot_'+f+'.csv') == False]

    # Giant gross main loop to do the photometry on images
    for f,file in enumerate(files):
        # Open your fits file and seperate the data and header.
        hdul = fits.open(directory+"/"+file)
        data = np.asarray(hdul[0].data)
        hdr = hdul[0].header
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
        # load up some header information that we'll be using
        exp = hdr['EXPTIME']
        if naver_on:
            naver = hdr['NAVER']
        # convert from DN to e-/s

        data = data * gain / exp


        # OPTION
        # Quick look at the histogram of the entire image
        # Should probably clean this up at some
        if img_hist == True:
            fig, ax = plt.subplot()
            ax.hist(data.flatten(),bins=5000,log=True)
            ax.set_title('Image Pixel Histogram')
            ax.set_xlabel('Pixel Value')
            ax.set_ylabel('e-/s')
            plt.show()

        # D:< ugly
        try:
            if telname == 'wise':
                HJD = hdr['JD-HELIO']
                MJD = hdr['JD'] - 2.4e6
            else:
                HJD = hdr['HJD']
                DATE = hdr['DATE-OBS']
                t = Time(DATE, scale='utc')
                MJD = t.mjd
        except:
            try:
                MJD = hdr['MJD']
                HJD = pyasl.helio_jd(MJD, coords[0][0], coords[0][1]) + 2.4e6
            except: 
                try:
                    MJD = hdr['MJD-OBS']
                    HJD = pyasl.helio_jd(MJD,coords[0][0], coords[0][1]) + 2.4e6       
                except:
                    try:
                        DATE = hdr['DATE-OBS']
                        t = Time(DATE, scale='utc')
                        MJD = t.mjd
                        HJD = pyasl.helio_jd(MJD, coords[0][0], coords[0][1]) + 2.4e6 
                    except:
                        try:
                            JD = hdr['JD']
                            HJD = pyasl.helio_jd(JD-2.4e6, coords[0][0], coords[0][1]) +2.4e6
                            MJD = JD - 2.4e6
                        except:
                            print('Failed to find date in header')
                            break
        hdul.close()
        HJD = HJD + (exp*.5)/(3600*24)
        # convert from arcsec to pixel values
        radii = np.asarray(rad)/pixsc
        boxsize = int(rad[3]/pixsc)

        # Need to find a clean solution for this - headers are not neat at all. WCS headers are even WORSE. If old WCS headers are not completely commented out astropy's WCS package messes up.
        # should add optional wcs to fix this? Find header cleaning solution. 
        # if it isn't found could run astrometry ? how would other people run it 
        # Some telescopes have bad wcs - or use astrometry to create a wcs. Headers ends up being super messy and astropy's wcs package will get confused if any wcs headers exist from other wcs routines. Best solution would be for the individual telescopes to fix this themselves or create a good routine to fix the wcs headers and add new wcs? 
        if tn == 'ratir':
            if img_filter == 'z':
                w = WCS(hdr[80:])
            else:
                w = WCS(hdr[45:])
        elif tn == 'caha':
            w = WCS(hdr[60:])
        # elif tn == 'lijiang':
        #     w = WCS(hdr[:60])
        else: 
            w = WCS(hdr)

        # better version of the above
        try:
            pixlist = w.all_world2pix(coords,1)
            pixlist[np.isnan(pixlist)] = 0
            # print(pixlist)
        except:
            print('Failed WCS')
            continue

# FIX HERE
        # data's shape
        datshape = np.shape(data)
        xmax, ymax = datshape[0], datshape[1]
        if c_func == '2dg':
            c_func = pu.centroid_2dg
        elif c_func == 'com':
            c_func = pu.centroid_com
        centers = []
        for pix in pixlist:
            try:
                # 
                if pix[0] < 0 or pix[1] < 0 or np.isnan(pix[0]):
                    centers.append((xmax/2.,ymax/2.))
                else:
                    sx,sy= pu.centroid_sources(data,xpos=pix[0],ypos=pix[1],box_size=boxsize,centroid_func=c_func)
                    # print(sx,sy)
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

        # cleaner version of the above chunk
        # px, py = np.transpose(pixlist)
        # sx, sy = pu.centroid_sources(data,xpos=px,ypos=py,box_size=boxsize,centroid_func=pu.centroid_com)
        # ssx, ssy = pu.centroid_sources(data,xpos=sx,ypos=sy,box_size=int(boxsize/2),centroid_func=pu.centroid_com)


        # Aperture photometry chunk, nothing fancy. Straight forward run an aperture for the object and an annulus ring 
        # checks = pu.CircularAperture(pixlist, r=radii[0])
        # agnap = pu.CircularAperture(pixlist[0], r=radii[0])
        aperture = pu.CircularAperture(centers, r=radii[0])
        annulus_aperture = pu.CircularAnnulus(centers, r_in=radii[1], r_out=radii[2])
        box = pu.RectangularAperture(pixlist,w=radii[3],h=radii[3])
        box2 = pu.RectangularAperture(pixlist,w=radii[3]/2,h=radii[3]/2)

        if centroid_show == True:
            fig, ax = plt.subplots(figsize=(15,15))
            # plt.subplot(projection=w)
            mean, std = np.nanmean(data),np.nanstd(data)
            lower, upper = mean-abs(std), mean+abs(std)
            ax.imshow(data, origin='lower',cmap='Greys_r',vmax=upper, vmin=lower)
            aperture.plot(color='yellow', lw=2, label='Centroid')
            # annulus_aperture.plot(color='red', lw=2)
            # agnap.plot(color='blue',lw=2)
            # checks.plot(color='green',lw=2, label='Initial Pixel Position')
            box.plot(color='red',lw=2)
            box2.plot(color='red', lw=2)
            ax.set_title(str(file))
            plt.show()

        # This will produce a phot table special to the astropy code, a lot like a pandas df
        phot = pu.aperture_photometry(data, aperture)
        # sky background estimation
        annulus_masks = annulus_aperture.to_mask(method='center')
        bkg_median, bkg_stdev = [], []
        # Allow for modular sky subtraction method
        for m,mask in enumerate(annulus_masks):
            # Do stats on the annulus ring to find a robust sky background
            # grabbing annulus data and flattening it
            annulus_data = mask.multiply(data)
            annulus_data_1d = annulus_data[mask.data > 0]
            # sigma clip default = 3.0
            mean, median_sigclip, stdev = sigma_clipped_stats(annulus_data_1d)
            if sky_hist == True:
                x, bins, __ = plt.hist(annulus_data_1d,bins=25)
                plt.title('Sky Annulus Histogram '+str(file)+' '+str(m))
                plt.vlines(mean,0,x.max(),color='red', label='Mean='+str(np.round(mean,4)))
                plt.vlines(median_sigclip,0,x.max(),color='green', label='Median='+str(np.round(median_sigclip,4)))
                plt.vlines((3*median_sigclip - 2*mean), 0, x.max(), color='black', label='DAOphot='+str(np.round(3*median_sigclip - 2*mean,4)))
                plt.xlabel('pixel value')
                plt.legend(frameon=False)
                plt.show()
            if sky_method == 'DAOphot':
                if mean >= median_sigclip:
                    bkg_median.append(3*median_sigclip - 2*mean)
                else:
                    bkg_median.append(mean)
            elif sky_method == 'clipped median':
                bkg_median.append(median_sigclip)
            bkg_stdev.append(stdev)
        bkg_median, bkg_stdev = np.array(bkg_median), np.array(bkg_stdev)
        # Like pandas df you can create new columns, neat! Perform that for various info we need
        phot['annulus_median'] = bkg_median
        phot['bkg_stdev'] = bkg_stdev
        phot['aper_bkg'] = bkg_median * aperture.area
        # Can do array operations so neatly
        phot['aper_sum_bkgsub'] = phot['aperture_sum'] - phot['aper_bkg']
        

        # Signal to noise is calculated from the general SNR equation and error is found from going backwards from it.
        SN, error = np.zeros(len(coords)), np.zeros(len(coords))
        for i in range(len(coords)):
            if phot['aper_sum_bkgsub'][i] <= 0 or phot['xcenter'][i].value == xmax/2.:
                phot['aper_sum_bkgsub'][i] = 0
                SN[i] = 0
            else:
                # I don't use dark eps I just assumed it was zero, easily possible to add it in though. There is a gain option, but I already did that way in the beginning so we don't need to input it here :)
                SN[i] = signal_to_noise_oir_ccd(t=exp, source_eps=phot['aper_sum_bkgsub'][i], sky_eps=phot['annulus_median'][i], dark_eps=0, rd=rdnoise, npix=np.pi*radii[0]**2)
                # SN[i] = signal_to_noise_oir_ccd(t=exp, source_eps=phot['aperture_sum'][i], sky_eps=0., dark_eps=0, rd=rdnoise, npix=np.pi*radii[0]**2)

                if SN[i] <= snr_cutoff:
                    phot['aper_sum_bkgsub'][i] = 0 
                    SN[i] = 0
                else:
                    error[i] = 1/SN[i]

        # bunch of tags that are useful to keep track of at this level. As the program moves forward we can prune them if needed, but they help us track down weird things.
        phot['telname'] = tn
        if tn == None:
            phot['MPC'] = tn
        else:
            phot['MPC'] = MPC[tn]
        

        
        if naver_on:
            phot['error'] = error/np.sqrt(naver)
            phot['naver_error'] = error/np.sqrt(naver)
        else:
            phot['error'] = error
        phot['SN'] = SN
        phot['time'] = HJD
        phot['modtime'] = MJD
        phot['file'] = file
        phot['filter'] = img_filter
        phot['expt'] = exp
        phot['Raw Counts'] = phot['aper_sum_bkgsub']*exp/gain
        pandaphot = phot.to_pandas()      
        pandaphot.to_csv(outdir+'/autophot/'+str(img_filter)+'/nightcsv/phot_'+file+'.csv')
        print(f+1, ' out of ', len(files), ' done')
    print('Photometry Complete')
    return


def phot_to_one(outdir,img_filter,coords):
    '''
    Neat function to go through the slew of single photometry files and seperate them by specifics like telescope, made easy by pandas! This can definitely be made sleeker!
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
    agnpd = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/AGNphots.csv')
    for name in telnames:
        agn = agnpd.loc[agnpd["MPC"] == name]
        agn.to_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_AGNphots.csv')
    return



def scale_AGN(outdir,img_filter,target_name,who=[],clean=False):
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
            #####################################################
            # look at this some more later
            std = np.nanstd(scaled_flux)
            print(std,median)
            nsig = abs((scaled_flux-median)/median)
            # plt.hist(nsig)
            # plt.show()
            
            # print(name)
            # print(nsig)
            # print('#####################')
            # if clean:
            # for p in range(len(scaled_flux)):
            #     if scaled_error[p] > scaled_flux[p]*.1 or (nsig[p] > 5): 
            #         scaled_error[p], scaled_flux[p] = np.nan, np.nan
            #     else:
            #         continue
            agnpd = pd.DataFrame()
            agnpd["scaled AGN flux"], agnpd["scaled AGN error"], agnpd["telname"], agnpd["MPC"], agnpd["time"], agnpd["modtime"], agnpd["file"] = scaled_flux, scaled_error, telname, MPC, time, modtime, files
            agnpd.to_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_scaledAGNphot.csv')
    return


def scale_plot(outdir,img_filter,target_name='NNF',cdict={},who=[]):
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
            s = datetime.datetime.now().strftime('%c')
            plt.figtext(0.005,.98,s=s)
            plt.savefig(outdir+'/autophot/'+str(img_filter)+'/'+name+'_'+img_filter+'_plot_'+str(img_filter)+'.jpg',dpi=300)
        plt.show()
    return

def dumb_align(outdir,img_filter,target_name,reference,cdict={},all_black=False):
    tel_pd = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/star_photometry_1.csv')
    tel_list = tel_pd["MPC"]
    telnames = set(tel_list)
    telnames = list(telnames)
    if reference in telnames:
        telnames.remove(reference)
    else:
        print("Invalid reference given - choosing random reference")
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

    plt.show()

    return




# def curve_elementary_align(outdir,img_filter,coords,target_name='NNF',reference='NA',cdict={}):
#     # Light curves will be aligned to the reference provided - if no reference is provided then a random name will be picked from the list
#     tel_pd = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/star_photometry_1.csv')
#     tel_list = tel_pd["MPC"]
#     telnames = set(tel_list)
#     telnames = list(telnames)
#     telnames.remove(reference)
#     r_pd = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/'+reference+'_scaledAGNphot.csv')
#     ref_flux = r_pd["scaled AGN flux"].to_numpy()
#     ref_error = r_pd["scaled AGN error"].to_numpy()
#     ref_time = r_pd["time"].to_numpy()
#     ref_name = r_pd["telname"].to_numpy()
#     ref_file = r_pd["file"].to_numpy()
#     ref_mpc = r_pd["MPC"]
#     ref_modtime = r_pd["modtime"].to_numpy()
#     fig, ax = plt.subplots(figsize=(10,10))
#     mpd = pd.DataFrame({"flux":ref_flux,"error":ref_error,"time":ref_time,"modtime":ref_modtime,"telname":ref_name,"MPC":ref_mpc,"file":ref_file})
#     if len(cdict) != 0:
#         color = cdict[ref_mpc[0]]
#     else:
#         r = random.random() 
#         b = random.random() 
#         g = random.random() 
#         color = (r,g,b) 
#     ax.errorbar(ref_time,ref_flux,yerr=ref_error,fmt='.',color=color,label=ref_mpc[0],ecolor=color,ms=2,elinewidth=.45)
#     ref_means = np.zeros(len(coords)-1)
#     for i in range(1,len(coords)):
#         prim = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/telphot/'+reference+'_star_phot'+str(i)+'_normed.csv')
#         prim_flux = prim["scaled flux"].to_numpy()
#         try:
#             ref_means[i-1] = np.nanmedian(prim_flux)
#         except:
#             ref_means[i-1] = np.nan
#     for name in telnames:
#         star_means = np.zeros(len(coords)-1)
#         for k in range(1,len(coords)):
#             # load star
#             star = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/telphot/'+name+'_star_phot'+str(k)+'_normed.csv')
#             # get star mean
#             star_flux = star["scaled flux"].to_numpy()
#             try:
#                 star_means[k-1] = np.nanmedian(star_flux)
#             except:
#                 star_means[k-1] = np.nan
#         star_means[star_means == 0] = np.nan
#         ratio_means = ref_means / star_means
#         ratio = np.nanmean(ratio_means)
#         # open scaled lightcurves 
#         tel_lightcurve = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_scaledAGNphot.csv')
#         flux = tel_lightcurve["scaled AGN flux"].to_numpy()
#         error = tel_lightcurve["scaled AGN error"].to_numpy()
#         time = tel_lightcurve["time"].to_numpy()
#         modtime = tel_lightcurve["modtime"].to_numpy()
#         tn = tel_lightcurve["telname"]
#         mpc = tel_lightcurve["MPC"]
#         fl = tel_lightcurve["file"]
#         # apply scale to flux and error
#         scdown_flux = flux * ratio
#         scdown_error = error * ratio
#         tpd = pd.DataFrame({"flux":scdown_flux,"error":scdown_error,"time":time,"modtime":modtime,"telname":tn,"MPC":mpc,"file":fl})
#         mpd = mpd.append(tpd)
#         if len(cdict) != 0:
#             color = cdict[mpc[0]]
#         else:
#             r = random.random() 
#             b = random.random() 
#             g = random.random() 
#             color = (r,g,b) 
#         ax.errorbar(time,scdown_flux,yerr=scdown_error,fmt='.',color=color,label=mpc[0],ecolor=color,ms=2,elinewidth=.45)
#     mpd.to_csv(outdir+'/autophot/'+str(img_filter)+'/'+target_name+'_aligned_lc.csv')
#     ax.set_ylabel("Flux")
#     ax.set_xlabel("HJD")
#     ax.set_title(str(target_name)+' '+str(img_filter))
#     ax.minorticks_on()
#     ax.tick_params(axis='y',direction='in',which='both')
#     ax.tick_params(axis='x',direction='in',which='both')
#     ax.xaxis.set_minor_locator(AutoMinorLocator())
#     ax.yaxis.set_minor_locator(AutoMinorLocator())
#     ax.tick_params(which='both',width=1.2)
#     ax.tick_params(which='major',length=6,labelsize=12)
#     ax.tick_params(which='minor',length=3)
#     ax.yaxis.set_ticks_position('both')
#     ax.xaxis.set_ticks_position('both')
#     ax.legend()
#     plt.savefig(outdir+'/autophot/'+str(img_filter)+'/'+str(target_name)+'_'+str(img_filter)+'_multitel_plot.jpg',dpi=900)
#     plt.show()


def clean_lc(outdir,img_filter):
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


def condense(outdir,img_filter,who=[],delta=0.25,counts=False):
    tel_pd = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/star_photometry_1.csv')
    tel_list = tel_pd["MPC"]
    telnames = set(tel_list)
    for name in telnames:
        if (name in who) == False and len(who) != 0:
            continue
        if counts:
            lc = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_scaledAGNphot.csv')
        else:
            lc = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_'+img_filter+'_f_lambda_AGNphot.csv')
        flux = lc['flux'].to_numpy()
        error = lc['error'].to_numpy()
        time = lc['time'].to_numpy()

        # lc = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_scaledAGNphot.csv')
        # flux = lc['scaled AGN flux'].to_numpy()
        # error = lc['scaled AGN error'].to_numpy()
        # time = lc['time'].to_numpy()


        telname = lc['telname']
        mpc = lc['MPC']
        cflux = []
        cerror = []
        ctime = []
        time_points = np.size(time)
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
            else:
                # load the slices
                sflux = np.copy(flux[i:wtop])
                serror = np.copy(error[i:wtop])
                stime = np.copy(time[i:wtop])
                weights = 1.0 / serror**2
                nflux = np.sum(weights * sflux) / np.sum(weights)
                cflux.append( nflux )
                cerror.append( np.sqrt( 1.0/ np.sum(weights)))
                ctime.append(np.mean(stime))
            i = i + wsize
        cflux = np.asarray(cflux)
        condensed_AGN = pd.DataFrame(cflux,columns=['flux'])
        condensed_AGN['error'] = np.asarray(cerror)
        condensed_AGN['time'] = np.asarray(ctime)
        condensed_AGN['MPC'] = mpc
        condensed_AGN['telname'] = telname
        condensed_AGN.to_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_uncal_cgs_condensed_lc.csv')
    return

def plot_condense(outdir,img_filter,target_name,cdict={},mpc={},who=[],multi=True,single=False,counts=False):
    tel_pd = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/star_photometry_1.csv')
    tel_list = tel_pd["MPC"]
    telnames = set(tel_list)
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
        if counts:
            clc = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_scaledAGNphot.csv')
        else:
            clc = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_uncal_cgs_condensed_lc.csv')
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
    plt.savefig(outdir+'/autophot/'+str(img_filter)+'/'+str(target_name)+'_'+img_filter+'_condense_multitel_plot.jpg',dpi=300)
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
                clc = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_uncal_cgs_condensed_lc.csv')
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
            plt.savefig(outdir+'/autophot/'+str(img_filter)+'/'+str(target_name)+'_'+str(name)+'_'+img_filter+'_condense_multitel_plot.jpg',dpi=300)
        plt.show()


def excess_var(outdir,img_filter,starlist,who=[]):
    tel_pd = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/star_photometry_1.csv')
    tel_list = tel_pd["telname"]
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
            print(old_error)
            print(new_error)
            name_AGN = pd.DataFrame()
            name_AGN["scaled AGN flux"], name_AGN["scaled AGN error"], name_AGN["telname"], name_AGN["MPC"], name_AGN["time"], name_AGN["modtime"], name_AGN["file"] = flux, new_error, telname, mpc, time, modtime, files
            name_AGN.to_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_scaledAGNphot.csv')
# def plot_exvar(outdir,AGN_RADEC,starlist,img_filter,target_name='NNF'):
#     tel_pd = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/star_photometry_0.csv')
#     tel_list = tel_pd["telname"]
#     telnames = set(tel_list)
#     for name in telnames:
#         fig, ax = plt.subplots(figsize=(10,10))
#         agn = pd.read_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_exvar_scaled_AGNphot.csv')
#         flux = agn["scaled AGN flux"].to_numpy()
#         error = agn["scaled AGN error"].to_numpy() 
#         time = agn["time"].to_numpy()
#         telname = agn["telname"].to_numpy()
#         flux = ma.masked_less_equal(flux,0)
#         if len(cdict) != 0:
#             color = cdict[]
#         else:
#             r = random.random() 
#             b = random.random() 
#             g = random.random() 
#             color = (r,g,b) 
#         ax.errorbar(time,flux,yerr=error,fmt='.',ecolor='grey',color=color,ms=2,elinewidth=.45)
#         ax.set_title(str(target_name)+' '+name+' '+str(img_filter), fontsize=24)
#         ax.minorticks_on()
#         ax.tick_params(axis='y',direction='in',which='both')
#         ax.tick_params(axis='x',direction='in',which='both')
#         ax.xaxis.set_minor_locator(AutoMinorLocator())
#         ax.yaxis.set_minor_locator(AutoMinorLocator())
#         ax.tick_params(which='both',width=1.2)
#         ax.tick_params(which='major',length=6,labelsize=12)
#         ax.tick_params(which='minor',length=3)
#         ax.set_xlabel('HJD',fontsize=16)
#         ax.set_ylabel('Flux',fontsize=16)
#         plt.savefig(outdir+'/autophot/'+str(img_filter)+'/'+str(target_name)+'_'+name+'_exvar_AGNplot_'+str(img_filter)+'.jpg',dpi=900)
#         plt.show()  

# def sn_plot():

#     for name in names:

#         for i in range(1,len(coords)):


#             time
#             sn
#             error = *100 #%error
#




# add at the end that the files get added to the pycali data directory!
def cali_format(outdir,flt,target_name,reference=None,condense=True):
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
                tf_lc = pd.read_csv(outdir+'/autophot/'+str(f)+'/'+n+'_uncal_cgs_condensed_lc.csv')
            else:
                tf_lc = pd.read_csv(outdir+'/autophot/'+str(f)+'/'+n+'_'+f+'_f_lambda_AGNphot.csv')
            tf_time, tf_flux, tf_error, tf_mpc, tf_telnames = tf_lc['time'].to_numpy(), tf_lc['flux'].to_numpy()*10**(14), tf_lc['error'].to_numpy()*10**(14), tf_lc['MPC'], tf_lc['telname']
            # tf_flux = tf_flux.map('%.6f')
            enc_time, enc_flux, enc_error = [], [], []
            for y in range(0,len(tf_time)):
                enc_time.append(tf_time[y])
                enc_flux.append(tf_flux[y])
                enc_error.append(tf_error[y])
            length = int(len(tf_time))
            try:
                string = '# '+tf_mpc[0]
            except:
                continue
            head = pd.DataFrame({'time': [string], 'flux': [length], 'error': [np.nan]})
            to_append = pd.DataFrame({'time':enc_time, 'flux':enc_flux, 'error':enc_error})
            to_append['flux'] = to_append['flux'].map('{:,.8e}'.format)
            to_append['error'] = to_append['error'].map('{:,.8e}'.format)
            to_append = to_append.sort_values(by=['time'],ascending=True)
            precali = head.append(to_append)
            cali_av_lc = cali_av_lc.append(precali)
        cali_av_lc.to_csv(outdir+'/autophot/'+target_name+'_'+f+'_cont_avg.txt',sep=' ',index=False,float_format='%.6f')
        with open(outdir+'/autophot/'+target_name+'_'+f+'_cont_avg.txt', 'r') as infile:
            temp = infile.read().replace("\"","")
            temp = temp.replace("time flux error\n","")
        with open(outdir+'/autophot/'+target_name+'_'+f+'_cont_avg.txt','w') as outfile:
            outfile.write(temp)
        print(f+' Averaged Continuum cali file created')
    return 
