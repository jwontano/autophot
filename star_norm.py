import lmfit as lft
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib
from astropy.stats import sigma_clipped_stats
import random as random
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from matplotlib import rc
import os, pandas, datetime, time
from matplotlib.ticker import FormatStrFormatter

# function that will be inputted into the minimizer - lmfit funcs require params to be passed in, independent variable(s) can also be inputted along with a functions args and kwargs. 
# residual function to be fed to lmfit
def chisqfunc(params,flux,error,length,mean=None):
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
	tel_star = pandas.read_csv(outdir+'/autophot/'+str(img_filter)+'/star_photometry_1.csv')
	tel_list = tel_star["MPC"]
	names = set(tel_list)
	names = list(names)
	for name in names:
		if (name in who) == False and len(who) != 0:
			continue
		else:
			all_flux, all_error, all_mean, all_scale = [], [], [], []
			for p in range(1,len(coords)):
				star = pandas.read_csv(outdir+'/autophot/'+str(img_filter)+'/telphot/'+name+'_star_phot'+str(p)+'.csv')
				flux, err = star["aper_sum_bkgsub"].to_numpy(), star["error"].to_numpy()
				flux[flux <= 0] = np.nan
				if len(flux) == 0 or np.sum(np.isnan(flux)) == np.size(flux):
					flux_median = 0.
				else:
					flux_median =  np.nanmedian(flux)
				scale = flux_median / flux
				flux_error = flux * err
				all_flux.append(flux)
				all_error.append(flux_error)
				all_mean.append(flux_median)
				all_scale.append(scale)
			avg_scale = []
			all_flux, all_error, all_mean, all_scale = np.asarray(all_flux), np.asarray(all_error), np.asarray(all_mean), np.asarray(all_scale)
			n_stars, n_pts = np.shape(all_flux)
			if n_pts == 1:
				bf = pandas.DataFrame([1.],columns=['scale'])
				bf['labels'] = ['A1']
				bf.to_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_star_chisq_scale.csv')
				print(name+' stars chi sq min fitted')
				continue
			else:
			# should add more restrictions on scaling factors - if only two scale factors are left or reasonable throw out the array?
				# There should be a way to transpose the arrays to make this work in much fewer lines of code
				for l in range(n_pts):
					star_scale = np.zeros(n_stars)
					for k in range(n_stars):
						star_scale[k] = (all_scale[k][l])
					if np.sum(np.isnan(star_scale)) == (n_stars):
						avg_scale.append(0.)
					elif np.nanstd(star_scale) > np.nanmean(star_scale)*.5:
						avg_scale.append(0.)
					else:
						# print(star_scale, np.nanmean(star_scale),np.nanmedian(star_scale), np.nanstd(star_scale))
						avg_scale.append(np.nanmedian(star_scale))
				params = np.concatenate((all_mean,avg_scale),axis=None)
				params = np.nan_to_num(params)
				
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
				if n_stars == 1:
					mini = lft.minimize(fcn=chisqfunc, params=pars, args=(all_flux,all_error,n_stars,all_mean[0]),nan_policy='omit')
				else:
					mini = lft.minimize(fcn=chisqfunc, params=pars,args=(all_flux,all_error,n_stars),nan_policy='omit')
				place = []
				labels = []
				for usls,param in mini.params.items():
					if 'star' in usls:
						continue
					else:
						place.append(param.value)
						labels.append(usls)
				bf = pandas.DataFrame(place,columns=["scale"])
				bf['labels'] = labels
				bf.to_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_star_chisq_scale.csv')
				print(name+' stars chi sq min fitted')

def star_norm(outdir,coords,img_filter,who=[],clean=False):
	# apply the nights scale to stars to check for flatness and create files that have the normalized values of the stars
	# then we can use to get for excess variance, error ratios etc. 
	tel_pd = pandas.read_csv(outdir+'/autophot/'+str(img_filter)+'/star_photometry_1.csv')
	tel_list = tel_pd["MPC"]
	telnames = list(set(tel_list))
	for name in telnames:
		if (name in who) == False and len(who) != 0:
			continue
		else:
			if os.path.exists(outdir+'/autophot/'+str(img_filter)+'/'+name+'_star_chisq_scale.csv') == False:
				print(name+' has no scale factor file - skipping')
				continue
			star_scale = pandas.read_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_star_chisq_scale.csv')
			applied_scale = star_scale["scale"].to_numpy()
			for z in range(1,len(coords)):
				cstar = pandas.read_csv(outdir+'/autophot/'+str(img_filter)+'/telphot/'+name+'_star_phot'+str(z)+'.csv')
				flux, telname, names, SN = cstar["aper_sum_bkgsub"].to_numpy(), cstar["telname"].to_numpy(), cstar["file"].to_numpy(), cstar["SN"].to_numpy()
				# this does not seem as clean as it should be 
				error = cstar["error"].to_numpy()
				error[error <= 0] = np.nan 
				# error[error > 0.2] = np.nan
				scaled_flux = flux * applied_scale
				scaled_flux[scaled_flux <= 0] = np.nan			
				scaled_error = error * applied_scale * flux
				scaled_error[scaled_error <= 0] = np.nan
				time = cstar["time"].to_numpy()
				if clean == True:
					mean, median, stdev = sigma_clipped_stats(scaled_flux)
					for p in range(len(scaled_flux)):
						if np.isnan(scaled_flux[p]):
							continue
						elif (scaled_error[p] > scaled_flux[p]*.1) or (np.abs(scaled_flux[p]-median)/stdev > 5): 
							scaled_flux[p], scaled_error[p] = np.nan, np.nan
						else:
							continue
				scaled_star = pandas.DataFrame(scaled_flux,columns=["scaled flux"])
				scaled_star["scaled flux error"], scaled_star["file"], scaled_star["SN"], scaled_star["telname"], scaled_star["time"] = scaled_error, names, SN, telname, time
				scaled_star.to_csv(outdir+'/autophot/'+str(img_filter)+'/telphot/'+name+'_star_phot'+str(z)+'_normed.csv')
			print(name+' stars normalized')

def normed_plots(outdir,img_filter,coords,cdict={},mpc={},who=[]):
	tel_pd = pandas.read_csv(outdir+'/autophot/'+str(img_filter)+'/star_photometry_1.csv')
	telnames = tel_pd["MPC"]
	telnames = set(telnames)
	for name in telnames:
		if (name in who) == False and len(who) != 0:
			continue

		fig, ax = plt.subplots(figsize=(10,10),ncols=2,nrows=len(coords)-1,sharex=True)
		plt.subplots_adjust(hspace=0,wspace=0) 
		for i in range(len(coords)-1):
			old = pandas.read_csv(outdir+'/autophot/'+str(img_filter)+'/telphot/'+name+'_star_phot'+str(i+1)+'.csv')
			oflux, otime = old['aper_sum_bkgsub'].to_numpy(), old['time']
			oerror = old['error'].to_numpy()*oflux
			oflux[oflux <= 0] = np.nan
			star =  pandas.read_csv(outdir+'/autophot/'+str(img_filter)+'/telphot/'+name+'_star_phot'+str(i+1)+'_normed.csv')
			flux, error, time = star["scaled flux"].to_numpy(), star["scaled flux error"].to_numpy(), star["time"]
			fig.patch.set_facecolor('gray')
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
				ax[1].set_ylabel("Flux",fontsize=16)
				ax[1].yaxis.set_label_position("right")
				ax[1].yaxis.tick_right()
				ax[0].minorticks_on()
				ax[1].minorticks_on()
				ax[0].xaxis.set_minor_locator(AutoMinorLocator())
				ax[0].yaxis.set_minor_locator(AutoMinorLocator())
				ax[1].xaxis.set_minor_locator(AutoMinorLocator())
				ax[1].yaxis.set_minor_locator(AutoMinorLocator())
				ax[0].tick_params(axis='y',direction='in',which='both')
				ax[0].tick_params(axis='x',direction='in',which='both')
				ax[1].tick_params(axis='y',direction='in',which='both')
				ax[1].tick_params(axis='x',direction='in',which='both')
				ax[1].yaxis.set_ticks_position('both')
				ax[1].xaxis.set_ticks_position('both')
				ax[0].yaxis.set_ticks_position('both')
				ax[0].xaxis.set_ticks_position('both')
			else:
				ax[i,1].axhline(y=np.nanmedian(oflux),linestyle='-.',color=color,alpha=.5)
				ax[i,1].errorbar(otime,oflux,yerr=oerror,fmt='.',color=color,ecolor=color,ms=ms,elinewidth=elw)
				ax[i,0].axhline(y=np.nanmedian(flux),linestyle='-.',color=color,alpha=.5)
				ax[i,0].errorbar(time,flux,yerr=error,fmt='.',color=color, ecolor=color,ms=ms,elinewidth=elw)
				ax[i,0].set_ylabel("Flux",fontsize=16)
				ax[i,1].set_ylabel("Flux",fontsize=16)
				ax[i,1].yaxis.set_label_position("right")
				ax[i,1].yaxis.tick_right()
				ax[i,0].minorticks_on()
				ax[i,1].minorticks_on()
				ax[i,0].xaxis.set_minor_locator(AutoMinorLocator())
				ax[i,0].yaxis.set_minor_locator(AutoMinorLocator())
				ax[i,1].xaxis.set_minor_locator(AutoMinorLocator())
				ax[i,1].yaxis.set_minor_locator(AutoMinorLocator())
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
			ax[1].set_title('Raw',fontsize=20)
			ax[0].set_title('Normed',fontsize=20)
		else:
			ax[-1,0].set_xlabel("HJD",fontsize=16)
			ax[-1,1].set_xlabel("HJD",fontsize=16)
			ax[0,1].set_title('Raw',fontsize=20)
			ax[0,0].set_title('Normed',fontsize=20)
		s=datetime.datetime.now().strftime('%c')
		plt.figtext(0.005,.98,s=s)
		fig.suptitle(name+' Comparison Stars '+img_filter)
		plt.savefig(outdir+'/autophot/'+str(img_filter)+'/'+name+'_'+img_filter+'_stars_normed.jpg',dpi=300)
	plt.show()


# Takes true magnitude values for specific comparison stars used in various filters used to put photometry measurements on a true flux scale (cgs units). 
# Standard flux to magnitude routine
# zp = zero point in flux 2.5 * log_10(wavelenght zeropoitn)
# fp = flux point of comparison star 10^(-(M-zp))/2.5

# I need a better way to feed these mags into the function
def true_mag_retrieve(outdir,who=None, mags=None, whom=[],feedback=False,f_to=[]):

	if mags != None:
		# tags = mags[0]
		
		tags, B, V, u, g, r, ip, z = mags[0], mags[1], mags[2], mags[3], mags[4], mags[5], mags[6], mags[7] 
	elif who == 'Mrk817': 
		tags, B, V, u, g, r, ip, z = [1,2,3,4,5],[16.349,16.218,16.274,13.558,15.29],[15.504,15.019,15.264,12.519,13.978],[17.00,17.56,17.10,15.18,16.51],[15.86,15.58,15.72,13.68,14.85],[15.49,14.86,15.25,12.70,14.01],[15.37,14.62,15.09,12.81,13.80],[15.35,14.51,15.04,13.31,13.61]
	elif who == 'Mrk335':
		tags = [1,2,4,5]
		B = [15.547,14.134,14.867,15.365]
		V = [14.849,13.509,14.141,14.23]
		u = [16.54,15.26,15.77,17.17]
		g = [15.167,13.83,14.487,14.88]
		r = [14.738,13.42,14.036,13.99]
		ip = [14.577,13.92,13.874,13.87]
		z = [14.55,13.22,13.86,13.35]
	elif who == 'NGC4395':
		tags = [1]
		B = [1]
		V = [1]
		g = [16.27]
		r = [15.65]
		ip =[15.43]
		z =[15.31]
		u = [1]
	elif who == 'IZ1':
		tags, B, V, u, g, r, ip, z = [1,2,3,4,5],[14.45,15.705,16.305,15.933,15.791],[13.73,14.766,15.336,15.007,15.151],[16.32,17.27,17.66,17.50,16.62],[14.30,15.51,15.84,15.48,14.41],[13.59,14.56,15.08,14.61,14.90],[13.39,14.32,15.22,14.61,15.16],[14.14,14.09,14.70,14.12,14.66]
	# elif who == 'NGC5548':
	# 	tags, B, V, u, g, r, ip, z = [2,3,8],[16.006,14.435,10.996],[15.381,13.774,10.464],[16.99,15.94,14.34],[15.63,15.79,10.72],[15.14,13.73,10.37],[14.98,15.07,10.23],[14.94,14.06,10.80]		
	elif who == 'Mrk1239':
		tags = [1,2,3,4]
		B = [16.28,15.551,15.179,15.321]
		V = [15.076,14.473,14.383,14.626]
		g = [15.462,15.018,14.748,14.951]
		r = [14.662,14.069,14.165,14.385]
		ip =[14.349,13.656,13.953,14.057]
		z =[13.91,13.565,13.529,13.443]
		u = [15.,15.,14.5,14.5]
	elif who=='IC4329A':
		tags = [1,2,3,4]
		B = [1,1,1,1]
		V = [1,1,1,1]
		# g = [15.999,15.673,16.11,16.055]
		g = [13.609,13.781,14.806,13.274]
		r = [1,1,1,1]
		ip =[1,1,1,1]
		z =[12.954,12.532,12.916,12.819]
		u = [1,1,1,1]		
	elif who == 'Ark120':
		tags = [1,2,3,4,5,6]
		B = [1,1,1,1,1,1]
		V = [11.921,11.846,11.786,12.319,12.639,12.481]
		u = [1,1,1,1,1,1]
		g = [12.132,12.196,12.472,12.479,12.707,12.719]
		r = [11.751,11.553,11.198,12.180,12.578,12.272]
		ip = [11.585,11.276,10.455,12.052,12.536,12.028]
		z = [11.591,11.157,10.104,12.061,12.571,11.939]

	zp = [2.5*np.log10(632*10**(-11)),2.5*np.log10(363.1*10**(-11)),2.5*np.log10(859.5*10**(-11)),2.5*np.log10(466.9*10**(-11)),2.5*np.log10(278.0*10**(-11)),2.5*np.log10(185.2*10**(-11)),2.5*np.log10(131.5*10**(-11))]
	zp, B, V, u, g, r, ip, z = np.asarray(zp),np.asarray(B),np.asarray(V),np.asarray(u),np.asarray(g),np.asarray(r),np.asarray(ip),np.asarray(z)
	fp = np.asarray([10**(-(B-zp[0])/2.5),10**(-(V-zp[1])/2.5),10**(-(u-zp[2])/2.5),10**(-(g-zp[3])/2.5),10**(-(r-zp[4])/2.5),10**(-(ip-zp[5])/2.5),10**(-(z-zp[6])/2.5)])
	if len(f_to) != 0:
		img_filters = f_to
		temp_fp = []
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
	else:
		img_filters = ['B','V','u','g','r','i','z']

	for h,img_filter in enumerate(img_filters):
		tel_pd = pandas.read_csv(outdir+'/autophot/'+str(img_filter)+'/star_photometry_1.csv')
		telnames = set(tel_pd["MPC"])
		for name in telnames:
			if (name in whom) == False and len(whom) != 0:
				continue
			else:
				mf = []
				for i in tags:
					s = pandas.read_csv(outdir+'/autophot/'+str(img_filter)+'/telphot/'+str(name)+'_star_phot'+str(i)+'_normed.csv')
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
					# plt.fill_between([1,5],mean-stdev,mean+stdev,alpha=.5)
					# plt.fill_between([1,5],mean-2*stdev,mean+2*stdev,alpha=.2)
					# plt.ylim(np.min(factor),np.max(factor))
					plt.ylabel('Scale Factor')
					plt.legend()
					plt.show()

				agn = pandas.read_csv(outdir+'/autophot/'+str(img_filter)+'/'+str(name)+'_scaledAGNphot.csv')
				agnflux = agn["scaled AGN flux"].to_numpy()
				true_agn = agnflux*mean
				agnerror = agn["scaled AGN error"].to_numpy()
				time, modtime, mpc, telname = agn["time"].to_numpy(), agn["modtime"].to_numpy(), agn["MPC"], agn["telname"]
				true_error = agnerror*mean
				true_agn[true_agn == 0] = np.nan
				# save the agn f_lambda into a csv 
				tagn = pandas.DataFrame(true_agn,columns=["flux"])
				tagn["error"], tagn["time"], tagn["modtime"], tagn["MPC"], tagn["telname"], tagn["file"] = true_error, time, modtime, mpc, telname, agn["file"]
				tagn.to_csv(outdir+'/autophot/'+str(img_filter)+'/'+str(name)+'_'+img_filter+'_f_lambda_AGNphot.csv')

def cgs_plot(cdict={},whoo=None,outdir=None,who=[],f_to=[],multi=True,single=False,all_black=False):
	if len(f_to) != 0:
		img_filters = f_to
	else:
		img_filters = ['B','V','u','g','r','i','z']
	ndict = {'u': 'u$^\prime$' , 'B': 'B', 'V': 'V', 'g': 'g$^\prime$', 'r': 'r$^\prime$', 'i': 'i$^\prime$', 'z': 'z$_s$'}
	for img_filter in img_filters:
		tel_pd = pandas.read_csv(outdir+'/autophot/'+str(img_filter)+'/star_photometry_1.csv')
		names = set(tel_pd["MPC"])
		fig, ax = plt.subplots(figsize=(10,10))
		ax.set_title(whoo+' '+ndict[img_filter],fontsize=24)
		ms, elw = 2, .45
		for name in names:
			if (name in who) == False and len(who) != 0:
				continue
			tagn = pandas.read_csv(outdir+'/autophot/'+str(img_filter)+'/'+str(name)+'_'+img_filter+'_f_lambda_AGNphot.csv')
			true_agn = tagn["flux"]
			true_error = tagn["error"]
			time = tagn["time"]
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
		plt.savefig(outdir+'/autophot/'+str(img_filter)+'/uncal_true_flux_'+whoo+'_'+str(img_filter)+'.pdf', dpi=900)
		plt.savefig(outdir+'/autophot/'+str(img_filter)+'/uncal_true_flux_'+whoo+'_'+str(img_filter)+'.jpg', dpi=300)
		plt.show()
		if single == True:
			for name in names:
				if (name in who) == False and len(who) != 0:
					continue
				tsagn = pandas.read_csv(outdir+'/autophot/'+str(img_filter)+'/'+str(name)+'_'+img_filter+'_f_lambda_AGNphot.csv')
				true_agn = tsagn["flux"]
				true_error = tsagn["error"]
				time = tsagn["time"]
				color = 'black'
				# if len(cdict) != 0:
				# 	color = cdict[name]
				# else:
				# 	r = random.random() 
				# 	b = random.random() 
				# 	g = random.random() 
				# 	color = (r,g,b)
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
				plt.savefig(outdir+'/autophot/'+str(img_filter)+'/uncal_true_flux_'+whoo+'_'+str(name)+'_'+str(img_filter)+'.pdf', dpi=900)
				plt.savefig(outdir+'/autophot/'+str(img_filter)+'/uncal_true_flux_'+whoo+'_'+str(name)+'_'+str(img_filter)+'.jpg', dpi=900)
			plt.show()

def error_ratio(outdir,coords,img_filter):
    tel_pd = pandas.read_csv(outdir+'/autophot/'+str(img_filter)+'/star_photometry_1.csv')
    tel_list = tel_pd["MPC"]
    telnames = set(tel_list)
    telnames = list(telnames)
    for name in telnames:
        star_norm_exvar = []
        star_means = []
        for i in range(0,len(coords)-1):
            star = pandas.read_csv(outdir+'/autophot/'+str(img_filter)+'/telphot/'+name+'_star_phot'+str(i+1)+'_normed.csv')
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





        #
        agn = pandas.read_csv(outdir+'/autophot/'+str(img_filter)+'/'+name+'_scaledAGNphot.csv')
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
        # lims = ax.margins()
        ax.fill_betweenx(y=[np.nanmin(star_norm_exvar)*.98,np.nanmax(star_norm_exvar)*1.02],x1=[amin,amin],x2=[amax,amax],color='red',alpha=0.3,label='Flux Range')
        plt.legend(frameon=False)
        # plt.savefig(outdir+'/autophot/'+str(img_filter)+'/'+str(name)+'_excess_scatter.jpg',dpi=600)

        # fig, ax = plt.subplots(figsize=(10,10))
        # ax.hist(star_norm_exvar)
        # ax.set_xlabel('Fractional Excess Scatter')

    plt.show()

