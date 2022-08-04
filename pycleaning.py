import pandas as pd
import numpy as np
import sys
import matplotlib.pyplot as plt

cdict = {'liverpool': 'orange', 'zowada': 'blue', 'wise': 'red', 'lijiang': 'y', 'ratir': 'seagreen', 'CAHA': 'cadetblue','F65': 'saddlebrown', 'V39': 'indigo', 'V37': 'green', 'W85': 'orchid', 'W86': 'indigo', 'W87': 'indianred', 'K93': 'olive', 'K91': 'darkorange', 'K92': 'y', 'Q63': 'sienna', 'Q64': 'chocolate','WISE': 'red','asassn': 'lime','ztf': 'k', 'Z24': 'purple', 'Z31': 'magenta','WMO': 'black', 'asiago': 'c'}
filters = ['u','B','V','g','r','i','z']

def produce_clean(filters):
	for filt in filters:
	    # load reconstruction and continuum
	    r_fn, c_fn = '/home/korbinite5646/AGN_home/MRK817/autophot/cali/save/'+filt+'_recon.txt','/home/korbinite5646/AGN_home/MRK817/autophot/cali/save/mrk817_'+filt+'_cont_avg.txt_cali'
	    recon_lc = pd.read_csv(r_fn,header=None,usecols=[0,3,5],delimiter=' ',names=['time','flux','error'])
	    cali_lc = pd.read_csv(c_fn,header=None,usecols=[0,1,3,5],delimiter=' ',names=['time','flux','error','tel'])
	    # calculate residuals by interpolation with the reconstruction
	    residual = cali_lc['flux'] - np.interp(cali_lc['time'], recon_lc['time'], recon_lc['flux'])
	    cali_lc['residual'] = residual
	    std_res = residual/cali_lc['error']
	    cali_lc['std_res'] = std_res
	    res_std = np.std(cali_lc['residual'])
	    cali_lc['std normed res'] = cali_lc['residual']/res_std
	    
	    clean_idx = np.where((abs(cali_lc['std normed res']) < 5) & (abs(cali_lc['std_res']) < 5))
	#     error_normed_idx = np.where(abs(cali_lc['std_res']) < 5)
	        
	    
	    tel_pd = pd.read_csv('/home/korbinite5646/AGN_home/MRK817/autophot/'+str(filt)+'/star_photometry_1.csv')
	    tel_list = tel_pd['MPC']
	    names = set(tel_list)
	    names = list(names)
	    if ('V37' in names) == True:
	            r = names.index('V37')
	            names[r], names[0] = names[0], names[r]
	    if filt == 'z':
	        r = names.index('zowada')
	        names[r], names[0] = names[0], names[r]
	        
	    cali_av_lc = pd.DataFrame(columns=['time', 'flux', 'error'])
	    for n in names:
	        tf_lc = pd.read_csv('/home/korbinite5646/AGN_home/MRK817/autophot/'+str(filt)+'/'+n+'_uncal_cgs_condensed_lc.csv')
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
	    cali_av_lc.to_csv('/home/korbinite5646/repository/PyCALI/data/mrk817_'+filt+'_cont_avg.txt',sep=' ',index=False)
	    with open('/home/korbinite5646/repository/PyCALI/data/mrk817_'+filt+'_cont_avg.txt', 'r') as infile:
	        temp = infile.read().replace("\"","")
	        temp = temp.replace("time flux error\n","")
	    with open('/home/korbinite5646/repository/PyCALI/data/mrk817_'+filt+'_cont_avg.txt','w') as outfile:
	        outfile.write(temp)
	    print(filt+' Averaged Continuum cali file created')



# def before_after(filters):
	


def check_clean(filters):
	for filt in filters:
	    r_fn, c_fn = '/home/korbinite5646/AGN_home/MRK817/autophot/cali/error/'+filt+'_recon.txt','/home/korbinite5646/AGN_home/MRK817/autophot/cali/error/mrk817_'+filt+'_cont_avg.txt_cali'
	    recon_lc = pd.read_csv(r_fn,header=None,usecols=[0,3,5],delimiter=' ',names=['time','flux','error'])
	    cali_lc = pd.read_csv(c_fn,header=None,usecols=[0,1,3,5],delimiter=' ',names=['time','flux','error','tel'])
	    names = list(set(cali_lc['tel']))
	    residual = cali_lc['flux'] - np.interp(cali_lc['time'], recon_lc['time'], recon_lc['flux'])
	    cali_lc['residual'] = residual
	    std_res = residual/cali_lc['error']
	    cali_lc['std_res'] = std_res
	    res_std = np.std(cali_lc['residual'])
	    cali_lc['std normed res'] = cali_lc['residual']/res_std
	# full histogram all points
	    
	    fig, ax = plt.subplots(figsize=(15,5))
	    ax.hist(residual, bins=25, color='black')
	    ax.set_title('Mrk 817 Residual Histogram '+filt)
	    ax.set_xlabel('Flux')
	    
	    fig, ax = plt.subplots(figsize=(15,5))
	    for n in names:
	        tel_indx = np.where(cali_lc['tel'] == n)
	        ax.errorbar(cali_lc['time'][tel_indx[0]],residual[tel_indx[0]],yerr=cali_lc['error'][tel_indx[0]],fmt='x',label=n,elinewidth=.7, color=cdict[n])
	    ax.set_title('Mrk 817 Residual over Time '+filt)
	    ax.set_ylabel('Flux')
	    ax.legend()
	# residual over time plot all tels together
	    fig, ax = plt.subplots(figsize=(15,5))
	    ax.hist(std_res, bins=25, color='gray')
	    ax.set_title('Mrk 817 Residual/Error Histogram '+ filt)
	    ax.set_xlabel('$\sigma$')
	    
	    fig, ax = plt.subplots(figsize=(15,5))
	    for n in names:
	        tel_indx = np.where(cali_lc['tel'] == n)
	        ax.scatter(cali_lc['time'][tel_indx[0]],std_res[tel_indx[0]],marker='x',label=n, color=cdict[n])
	    ax.set_title('Mrk 817 Residual/Error '+filt)
	    ax.set_ylabel('$\sigma$')
	    ax.legend()
	    
	    fig, ax = plt.subplots(figsize=(15,5))
	    ax.hist(cali_lc['std normed res'], bins=25, color='gray')
	    ax.set_title('Mrk 817 Residual/Std Histogram '+filt)
	    ax.set_xlabel('$\sigma$')
	    
	    fig, ax = plt.subplots(figsize=(15,5))
	    for n in names:
	        tel_indx = np.where(cali_lc['tel'] == n)
	        ax.scatter(cali_lc['time'][tel_indx[0]],cali_lc['std normed res'][tel_indx[0]],marker='x',label=n, color=cdict[n])
	    ax.set_title('Mrk 817 Residual/Std '+filt)
	    ax.set_ylabel('$\sigma$')
	    ax.legend()
	    plt.show()



###########################################################################################################################
wtd = sys.argv[1]
filters = sys.argv[2]

if wtd == '1':
	produce_clean(filters)
elif wtd == '2':
	check_clean(filters)
else:
	print('give correct input')