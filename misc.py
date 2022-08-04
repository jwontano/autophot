import os 
from astropy.io import fits
import sys
import shutil
from tkinter import filedialog
from tkinter import *

def filter_sort_images(directory,lco=False,outdir=None):
	f_dict = {'gp': 'g', 'g': 'g', 'SDSS-G': 'g', 'ip': 'i', 'i': 'i', 'SDSS-I': 'i', 'rp': 'r', 'r': 'r', 'SDSS-R': 'r', 'R': 'r', 'up': 'u', 'u': 'u', 'SDSS-U': 'u', 'zs': 'z', 'z': 'z', 'SDSS-Z': 'z', 'V': 'V', 'Bessell-V': 'V', 'B': 'B', 'Bessell-B': 'B', "Y": 'Y'}
	files = [f for f in os.listdir(directory) if os.path.isfile(directory+'/'+f) and f.endswith(".fits")]
	for file in files:
		if lco:
			image = directory+'/'+file
			hdul = fits.open(image)
			hdr = hdul[0].header
			target = hdr['OBJECT']
			prop = hdr['PROPID']
			fltr = hdr['FILTER']
			fltr = f_dict[fltr]
			hdul.close()
			if os.path.exists('/home/korbinite5646/AGN_home/FITS/'+prop) == False:
				print('Creating proposal storage directory for '+prop)
				os.mkdir('/home/korbinite5646/AGN_home/FITS/'+prop)
			if os.path.exists('/home/korbinite5646/AGN_home/FITS/'+prop+'/'+target)==False:
				print('Creating directory for '+target+' in '+prop)
				os.mkdir('/home/korbinite5646/AGN_home/FITS/'+prop+'/'+target)
			if os.path.exists('/home/korbinite5646/AGN_home/FITS/'+prop+'/'+target+'/'+fltr) == False:
				print('Creating filter storage directory for '+fltr)
				os.mkdir('/home/korbinite5646/AGN_home/FITS/'+prop+'/'+target+'/'+fltr)
			try:
				shutil.move(image,'/home/korbinite5646/AGN_home/FITS/'+prop+'/'+target+'/'+fltr+'/'+file)
				print(file+' moved to '+'/home/korbinite5646/AGN_home/FITS/'+prop+'/'+target+'/'+fltr)
			except:
				print('Failed to move file - dunno why')
				continue
		else:
			image = directory+'/'+ file
			hdul = fits.open(image)
			hdr = hdul[0].header
			try:
				fltr = hdr['FILTER']
			except:
				fltr = hdr['FILTER1']
			fltr = f_dict[fltr]
		try:
			shutil.move(image,outdir+'/'+fltr+'/'+file)
			print(file+' moved to '+outdir+'/'+fltr+'/'+file)
		except:
			print('Failed to move file :-(')
			continue
	return


def ask():
	root = Tk()
	root.withdraw()
	folder_selected = filedialog.askdirectory()
	return folder_selected



def main():
	d = '/home/korbinite5646/AGN_home/FITS/KEY2020B-006/Mrk 817/other_tels/liverpool/'
	filter_sort_images(d,False,outdir=d)
	# direcs = ['20220103','20220104','20220107','20220108','20220109','20220110','20220111','20220114','20220115','20220117','20220118','20220120','20220121','20220122','20220124','20220125','20220127','20220128','20220129','20220130','20220201','20220202','20220203','20220205','20220206','20220207','20220208','20220209','20220210','20220211','20220212','20220223','20220224','20220228','20220301','20220304','20220306','20220310','20220312','20220317','20220318','20220321','20220322','20220323','20220324','20220326','20220327','20220328','20220329','20220330']

	# for d in direcs:
	# 	filter_sort_images('/home/korbinite5646/AGN_home/FITS/KEY2020B-006/Mrk 817/other_tels/mrk817_wise/'+d,False,'/home/korbinite5646/AGN_home/FITS/KEY2020B-006/Mrk 817/other_tels/mrk817_wise/')
main()