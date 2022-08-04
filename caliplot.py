
import matplotlib.pyplot as plt
from astropy.io import fits
import numpy as np
import numpy.ma as ma
import pandas as pd
from matplotlib import rc
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoMinorLocator)
import random as random
import sys, datetime, csv, time
from scipy import stats

# ask where to find the cali files or give the cali files?
def multiplot(tn,cdict,avg=False):
	rc('font', **{'family':'Dejavu Sans'})
	rc('text', usetex=True)
	filters = ['u','B','g','V','r','i','z']
	if tn == 'mrk817' or tn == 'Mrk817':
		ttn = 'Mrk 817'
	elif tn == 'mrk335':
		ttn = 'Mrk 335'
	else:
		ttn = tn
	fig, ax = plt.subplots(nrows=len(filters),ncols=1,sharex=True,figsize=(15,10))
	plt.subplots_adjust(hspace=0)
	ax[0].set_title(ttn,fontsize=28)
	for i,f in enumerate(filters):
		if ttn == 'Mrk 817':
			fn = '/home/korbinite5646/AGN_home/MRK817/autophot/cali/error/mrk817_'+str(f)+'_cont_avg.txt_cali'
		elif ttn == 'Mrk 335':
			fn = 'E:/AGN_central/LightCurveAutomation/mrk335/dev/autophot/cali/error/mrk335_'+str(f)+'_cont_avg.txt_cali'
		elif ttn == 'IZ1':
			fn = 'E:/AGN_central/LightCurveAutomation/IZ1/autophot/cali/IZ1_'+str(f)+'_cont_avg.txt_cali'
		else:
			fn = '/home/korbinite5646/AGN_home/'+tn+'_'+str(f)+'_cont.txt_cali'
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

			if n == 'asiago':
				as_v = pd.read_csv('/home/korbinite5646/AGN_home/MRK817/autophot/'+f+'/asiago_uncal_cgs_condensed_lc.csv')
				lc_size = np.size(flux)
				error = (as_v['error']*10**(14))[0:lc_size]

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
	if ttn == 'Mrk 817':
		fig.legend(handles,labels,markerscale=3.5,borderaxespad=0,handletextpad=0.0,prop={'size': 10},bbox_transform=plt.gcf().transFigure,bbox_to_anchor=(.445,.952),frameon=False, ncol=5)
	elif ttn == 'Mrk 335':
		fig.legend(handles,labels,borderaxespad=0,handletextpad=0.0,prop={'size': 8},loc='upper left',bbox_to_anchor=(.12,.8795),bbox_transform=plt.gcf().transFigure, framealpha=0.7,ncol=10,markerscale=2,frameon=False)
	elif tn == 'IZ1':
		fig.legend(handles,labels,markerscale=1.5,borderaxespad=0,handletextpad=0.0,prop={'size': 10},loc='upper left',bbox_to_anchor=(.12,.88),bbox_transform=plt.gcf().transFigure, ncol=6, frameon=False)
	else:
		fig.legend(handles,labels,markerscale=1.5,borderaxespad=0,handletextpad=0.0,prop={'size': 10},loc='upper left',bbox_to_anchor=(.12,.88),bbox_transform=plt.gcf().transFigure, ncol=10, frameon=False)
	ax[len(filters)-1].set_xlabel('HJD',fontsize=20)
	ax[3].set_ylabel( 'f$_\lambda$ ($10^{-14}$ erg cm$^{-2}$ s$^{-1}$ $\mathrm{\AA}^{-1}$)',fontsize=20)
	ax[3].yaxis.set_label_coords(-.05,.5)
	s = datetime.datetime.now().strftime('%c')
	plt.figtext(0.005,.98,s=s)
	# if avg == True:
	# 	# plt.savefig('/home/korbinite5646/AGN_home/'+tn+'/newb/autophot/caliplot/'+tn+'cali_avg.pdf',dpi=600)
	plt.savefig('/home/korbinite5646/AGN_home/MRK817/autophot/caliplot/'+tn+'cali_avg.jpg',dpi=300)
	# else:
	# 	# plt.savefig('/media/korbinite5646/backup/research_backup/AGN_central/LightCurveAutomation/'+tn+'/collab/autophot/caliplot/'+tn+'cali.pdf',dpi=300)
		# plt.savefig('/media/korbinite5646/backup/research_backup/AGN_central/LightCurveAutomation/'+tn+'/collab/autophot/caliplot/'+tn+'cali.jpg',dpi=300)
	plt.show()

def plot(tn,cdict,avg=False):
	rc('font', **{'family':'Dejavu Sans'})
	rc('text', usetex=True)
	filters = ['u','B','g','V','r','i','z']
	ms = 2
	elw = .45
	if tn == 'mrk817':
		ttn = 'Mrk 817'
	elif tn == 'mrk335':
		ttn = 'Mrk 335'
	else: 
		ttn = tn
	for i,f in enumerate(filters):
		if ttn == 'Mrk 817':
			fn =  '/home/korbinite5646/AGN_home/MRK817/autophot/cali/error/mrk817_'+str(f)+'_cont_avg.txt_cali'
		elif ttn == 'Mrk 335':
			fn = 'E:/AGN_central/LightCurveAutomation/mrk335/diff_stars/autophot/cali/mrk335_'+str(f)+'_cont_avg.txt_cali'
		elif ttn == 'IZ1':
			fn = 'E:/AGN_central/LightCurveAutomation/IZ1/autophot/cali/IZ1_'+str(f)+'_cont.txt_cali'
		cali = pd.read_csv(fn,header=None,usecols=[0,1,3,5],delimiter=' ',names=['time','flux','error','names'])
		names = set(cali['names'])
		fig, ax = plt.subplots(figsize=(10,10))
		ax.set_title(ttn+' '+f,fontsize=24)
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
			if n == 'V39' or n == 'V37':
				ax.errorbar(time,flux,yerr=error,fmt='.',ms=ms,ecolor=color,elinewidth=elw,color=color,alpha=alpha,label='LCO-'+n)
			elif n == 'zowada':
				ax.errorbar(time,flux,yerr=error,fmt='.',ms=ms,ecolor=color,elinewidth=elw,color=color,alpha=alpha,label='Zowada')
			elif n == 'ratir':#
				ax.errorbar(time,flux,yerr=error,fmt='.',ms=ms,ecolor=color,elinewidth=elw,color=color,label='SPM',alpha=alpha)
			elif n == 'wise' or n == 'WISE':
				ax.errorbar(time,flux,yerr=error,fmt='.',ms=ms,ecolor=color,elinewidth=elw,color=color,label='Wise',alpha=alpha)
			elif n == 'lijiang':
				ax.errorbar(time,flux,yerr=error,fmt='.',ms=ms,ecolor=color,elinewidth=elw,color=color,label='Lijiang',alpha=alpha)
			elif n == 'liverpool':
				ax.errorbar(time,flux,yerr=error,fmt='.',ms=ms,ecolor=color,elinewidth=elw,color=color,label='Liverpool',alpha=alpha)
			elif n == 'F65':
				ax.errorbar(time,flux,yerr=error,fmt='.',ms=ms,ecolor=color,elinewidth=elw,color=color,alpha=alpha,label='FTN-M3')
			elif n == 'CAHA':
				ax.errorbar(time,flux,yerr=error,fmt='.',ms=ms,ecolor=color,elinewidth=elw,color=color,label='CAHA',alpha=alpha)
			elif n == 'ztf':
				ax.errorbar(time,flux,yerr=error,fmt='.',ms=ms,ecolor='b',elinewidth=elw,color='b',label='ZTF')
			else:
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
		if ttn == 'Mrk 335':
			fig.legend(handles,labels,borderaxespad=0,handletextpad=0.0,prop={'size': 12},loc='center',bbox_to_anchor=(.5,.862),bbox_transform=plt.gcf().transFigure, framealpha=0.8,ncol=len(names),markerscale=2,frameon=False)
		elif ttn == 'Mrk 817':
			fig.legend(handles,labels,borderaxespad=0,handletextpad=0.0,prop={'size': 12},loc='center',bbox_to_anchor=(.5,.862),bbox_transform=plt.gcf().transFigure, framealpha=0.8,markerscale=2,ncol=len(names),frameon=False)
		elif ttn == 'IZ1':
			fig.legend(handles,labels,borderaxespad=0,handletextpad=0.0,prop={'size': 12},loc='center',bbox_to_anchor=(.5,.85),bbox_transform=plt.gcf().transFigure,ncol=5,markerscale=2.5,frameon=False)
		ax.set_xlabel('HJD',fontsize=16)
		ax.set_ylabel( 'f$_\lambda$ ($10^{-14}$ erg cm$^{-2}$ s$^{-1}$ $\mathrm{\AA}^{-1}$)',fontsize=16)
		ax.yaxis.set_label_coords(-.07,.5)
		s = datetime.datetime.now().strftime('%c')
		plt.figtext(0.005,.98,s=s)
		# plt.savefig(tn+'_'+f+'_cali_avg.pdf',dpi=300)
	plt.show()

def ccf_plot(tn,ref):
	rc('font', **{'family':'Dejavu Sans'})
	rc('text', usetex=True)
	cf = {'u': 'violet','B': 'purple','g': 'blue', 'V': 'green', 'r': 'yellow', 'i': 'orange', 'z': 'red'}

	if ref == 'u':
		filters = ['B','g','V','r','i','z']
	elif ref == 'B':
		filters = ['u','g','V','r','i','z']
	elif ref == 'g':
		filters = ['u','B','V','r','i','z']
	perclim = 84.1344746
	clag = []
	cerror = [[],[]]
	fig,ax = plt.subplots(figsize=(10,10))
	for f in filters:
		lag,re = np.loadtxt('E:/AGN_central/LightCurveAutomation/'+tn+'/autophot/pyccf/'+tn+'_'+ref+'-'+f+'_ccf_yap.dat', unpack=True, usecols=[0,1])
		cent = np.loadtxt('E:/AGN_central/LightCurveAutomation/'+tn+'/autophot/pyccf/'+tn+'_'+ref+'-'+f+'_centtab_yap.dat', unpack=True, usecols=[0])
		peak = np.loadtxt('E:/AGN_central/LightCurveAutomation/'+tn+'/autophot/pyccf/'+tn+'_'+ref+'-'+f+'_peaktab_yap.dat', unpack=True, usecols=[0])
		ax.plot(lag,re,label=f,color=cf[f])
		centau = stats.scoreatpercentile(cent,50)
		centau_uperr = (stats.scoreatpercentile(cent, perclim))-centau
		centau_loerr = centau-(stats.scoreatpercentile(cent, (100.-perclim)))
		clag.append(centau)
		cerror[0].append(centau_loerr)
		cerror[1].append(centau_uperr)
	if tn == 'mrk817':
		tn = 'Mrk 817'
	elif tn == 'mrk335':
		tn = 'Mrk 335'
	ax.legend(frameon=False,loc='upper left')
	ax.set_xlabel('Lag (days)',fontsize=16)
	ax.set_ylim(0,1)
	ax.set_xlim(-20,20)
	ax.set_ylabel('CCF',fontsize=16)
	ax.set_title(tn+' '+ref+' CCF',fontsize=24)
	ax.minorticks_on()
	ax.tick_params(axis='y',direction='in',which='both')
	ax.tick_params(axis='x',direction='in',which='both')
	ax.xaxis.set_minor_locator(AutoMinorLocator())
	ax.yaxis.set_minor_locator(AutoMinorLocator())
	ax.tick_params(which='both',width=1.2)
	ax.tick_params(which='major',length=6,labelsize=12)
	ax.tick_params(which='minor',length=3)
	# ax.yaxis.tick_right()
	ax.yaxis.set_ticks_position('both')
	ax.xaxis.set_ticks_position('both')
	s = datetime.datetime.now().strftime('%c')
	plt.figtext(0.005,.98,s=s)
	if tn == 'Mrk 817':
		ntn = 'mrk817'
	elif tn == 'Mrk 335':
		ntn = 'mrk335'
	else:
		ntn = tn
	# plt.savefig('E:/AGN_central/LightCurveAutomation/'+ntn+'/autophot/ccfplot/'+tn+'_'+ref+'_ccf.pdf',dpi=900)
	plt.savefig('E:/AGN_central/LightCurveAutomation/'+ntn+'/autophot/ccfplot/'+tn+'_'+ref+'_ccf.jpg',dpi=300)
	plt.show()
	if ref == 'u':
		angs = [4380,4830,5450,6260,7670,9100]
	elif ref == 'B':
		angs = [3560,4830,5450,6260,7670,9100]
	elif ref == 'g':
		angs = [3560,4380,5450,6260,7670,9100]
	fig, ax = plt.subplots(figsize=(10,10))
	ax.scatter(angs,clag,marker='D',color='black',s=50)
	ax.set_xlabel('$\lambda(\mathrm{\AA})$',fontsize=16)
	ax.set_ylabel('Lag (days)',fontsize=16)
	ax.set_title(tn+' '+ref+' Lags',fontsize=24)
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
	plt.figtext(0.005,.98,s=s)
	# plt.savefig('E:/AGN_central/LightCurveAutomation/'+ntn+'/autophot/ccfplot/'+tn+'_'+ref+'_lambda_lag.pdf',dpi=900)
	plt.savefig('E:/AGN_central/LightCurveAutomation/'+ntn+'/autophot/ccfplot/'+tn+'_'+ref+'_lambda_lag.jpg',dpi=300)
	plt.show()
	return


def main(name=None):
	fdict = {'u': 'u$^\prime$' , 'B': 'B', 'V': 'V', 'g': 'g$^\prime$', 'r': 'r$^\prime$', 'i': 'i$^\prime$', 'z': 'z$^\prime$'}
	cdict = {'liverpool': 'orange', 'zowada': 'blue', 'wise': 'red', 'lijiang': 'y', 'ratir': 'seagreen', 'CAHA': 'cadetblue','F65': 'saddlebrown', 'V39': 'indigo', 'V37': 'green', 'W85': 'orchid', 'W86': 'indigo', 'W87': 'indianred', 'K93': 'olive', 'K91': 'darkorange', 'K92': 'y', 'Q63': 'sienna', 'Q64': 'chocolate','WISE': 'red','asassn': 'lime','ztf': 'k', 'Z00': 'purple', 'Z01': 'magenta', 'asiago': 'c'}
	# multiplot(name,cdict)
	# plot(name,cdict)
	multiplot(name,cdict,avg=True)
	plot(name,cdict)
name = sys.argv[1]
main(name)