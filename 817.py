
'''
Main python code to run all of autophot for mrk817, includes all telescopes/observatories. 
'''
from star_norm import *
from AGN_phot import *

def main():
	
	# img_filters = ['B','V','u','g','r','i','z'] 
	# directory = '/home/korbinite5646/AGN_home/FITS/KEY2020B-006/Mrk 817/'
	# for img_filter in img_filters:
	# 	phot( directory+img_filter,
	# 	outdir,
	# 	coords=coords,
	# 	img_filter=img_filter,
	# 	target_name=target_name,
	# 	telname='lco',
	# 	snr_cutoff=10.,
	# 	renew=True,
	# 	MPC=MPC)


	# img_filters = ['g','r','i','z']
	# directory = '/home/korbinite5646/AGN_home/FITS/KEY2020B-006/Mrk 817/other_tels/mrk817_wise/'
	# for img_filter in img_filters:
	# 	phot( directory+img_filter,
	# 	outdir,
	# 	coords=coords,
	# 	img_filter=img_filter,
	# 	target_name=target_name,
	# 	telname='wise',
	# 	rdnoise=8.0,
	# 	pixsc=0.882,
	# 	gain=0.476,
	# 	centroid_show=False)
# z band test star ratios with FTN and WISE
	img_filters = ['u']#['g','r','i']
	# directory = '/home/korbinite5646/AGN_home/FITS/KEY2020B-006/Mrk 817/other_tels/zowada/'
	# for img_filter in img_filters:
	# 	phot(directory+img_filter,
	# 		outdir,
	# 		coords=coords, 
	# 		img_filter=img_filter,
	# 		target_name=target_name,
	# 		telname='zowada',
	# 		rdnoise=8.0,
	# 		pixsc=0.882,
	# 		gain=0.476,
	# 		naver_on=True,
	# 		rad=[7.,15.,20.,15.],
	# 		centroid_show=False,
	# 		MPC=MPC,
	# 		c_func='2dg',
	# 		sky_method='clipped median')

# # # 	# Do centroid show for liverpool u band?

	# img_filters = ['g','r','i','z']
	# directory = '/home/korbinite5646/AGN_home/FITS/KEY2020B-006/Mrk 817/other_tels/liverpool/'
	# for img_filter in img_filters:
	# 	phot(directory+img_filter+'/ast/',
	# 		outdir,
	# 		coords=coords,
	# 		img_filter=img_filter,
	# 		target_name=target_name,
	# 		telname='liverpool',
	# 		rdnoise=8.00,
	# 		pixsc=0.3037,
	# 		gain=1.62,
	# 		MPC=MPC,
	# 		centroid_show=False)

	# img_filters = ['B','V']
	# directory = '/home/korbinite5646/AGN_home/FITS/KEY2020B-006/Mrk 817/other_tels/WMO/'
	# for img_filter in img_filters:
	# 	phot(directory+img_filter,
	# 		outdir,
	# 		coords=coords,
	# 		img_filter=img_filter,
	# 		target_name=target_name,
	# 		telname='WMO',
	# 		rdnoise=12.0,
	# 		pixsc=0.49,
	# 		gain=1.4142135623731,
	# 		centroid_show=False,
	# 		MPC=MPC)

	# img_filters = ['u']#['V','u','r','i']
	# directory = '/home/korbinite5646/AGN_home/FITS/KEY2020B-006/Mrk 817/other_tels/asiago/'
	# for img_filter in img_filters:
	# 	phot(directory+img_filter,
	# 		outdir,
	# 		coords=coords,
	# 		img_filter=img_filter,
	# 		telname='asiago',
	# 		rdnoise=10.5,
	# 		pixsc=0.87,
	# 		gain=1.55,
	# 		centroid_show=False) 

	# img_filters = ['V']
	# directory = '/home/korbinite5646/AGN_home/FITS/KEY2020B-006/Mrk 817/other_tels/lj/images/ast/'
	# phot(directory,
	# 	outdir,
	# 	coords=coords,
	# 	img_filter='V',
	# 	telname='lijiang',
	# 	rdnoise=13.5,
	# 	pixsc=0.286,
	# 	gain=3.2,
	# 	MPC=MPC,
	# 	centroid_show=False)
	# directory = '/home/korbinite5646/AGN_home/FITS/KEY2020B-006/Mrk 817/other_tels/caha/images/ast/'
	# phot(directory,
	# 	outdir,
	# 	coords=coords,
	# 	img_filter='V',
	# 	telname='caha',
	# 	rdnoise=7.4,
	# 	pixsc=0.53,
	# 	gain=1.4,
	# 	MPC=MPC,
	# 	centroid_show=False)
	
	# img_filters = ['V','u','g','r','i','z'] 
	# img_filters = ['V']
	who = ['zowada']
	# who = ['lijiang','CAHA']
	# for img_filter in img_filters:
		# phot_to_one(outdir,img_filter=img_filter,coords=coords)
		# fit(outdir,coords,img_filter,who=who)
		# star_norm.star_norm(outdir,coords,img_filter,who=who,clean=False)
		# star_norm.normed_plots(outdir,img_filter,coords=coords,cdict=cdict,mpc=MPC,who=who)
		# scale_AGN(outdir, img_filter=img_filter, target_name=target_name,clean=False,who=who)
		# scale_plot(outdir,img_filter=img_filter,target_name=target_name,cdict=cdict,who=who)
		# clean_lc(outdir=outdir,img_filter=img_filter)
		# dumb_align(outdir,img_filter=img_filter,target_name=target_name,reference='zowada',cdict=cdict)	



# 



	# true_mag_retrieve(outdir,
		# who='Mrk817',
		# f_to=img_filters,
		# whom='zowada',
		# mags=[[1,2,3,4],[],[],[15.1,15.24,15.22,15.27],[],[],[],[]])
	# cgs_plot(cdict=cdict,whoo=target_name,outdir=outdir,f_to=['u'],single=False)
	for img_filter in img_filters:
		condense(outdir,img_filter,who=who)
		# plot_condense(outdir,img_filter,target_name,cdict,MPC,single=False)
	# img_filters = ['B','V']
	# cali_format(outdir,img_filters,target_name,reference='V37')
	# img_filters = ['g','r','i','z']
	# cali_format(outdir,img_filters,target_name,reference='V37')

	cali_format(outdir,['u'],target_name,reference='zowada')

# asiago?
# coords = [ [219.09200836,58.7942765], [219.077094,58.677452], [219.05994,58.666196]]

# ZOWADA 4
# coords = [[219.09200836,58.7942765],[219.156599,58.71755],[219.077094,58.677452]]

# LT 1 2 
# coords = [[219.09200836,58.7942765],[219.1357,58.83312], [219.17719,58.84183] ]

# ASIAGO



# original 
# coords = [ [219.09200836,58.7942765], [219.1357,58.83312], [219.17719,58.84183] , [219.150785,58.766423] , [219.077094,58.677452], [219.05994,58.666196]] 


# jakes 
# 15.10 15.24 15.22 15.27
coords = [[219.09200836,58.7942765], [218.7154167,58.6825000],[219.1550000,58.7172222],[219.4291667,58.7075000],[219.3066667,58.6644444]]

target_name = 'Mrk817'
outdir = '/home/korbinite5646/AGN_home/MRK817/'
cdict={'liverpool': 'black',
'zowada': 'blue',
'WISE': 'red',
'wise': 'red',
'lijiang': 'tan',
'ratir': 'darkseagreen',
'CAHA': 'cadetblue',
'F65': 'saddlebrown',
'V39': 'darkorange',
'V37': 'darkgreen',
'W85': 'orchid',
'WMO': 'black',
'Z24': 'orange',
'Z31': 'pink',
'asiago': 'c'}

MPC ={'ogg-clma-2m0a': 'F65', 
 'ogg-clma-0m4b': 'T04', 
 'ogg-clma-0m4c': 'T03', 
 'elp-doma-1m0a': 'V37', 
 'elp-domb-1m0a': 'V39', 
 'elp-aqwa-0m4a': 'V38', 
 'lsc-doma-1m0a': 'W85', 
 'lsc-domb-1m0a': 'W86', 
 'lsc-domc-1m0a': 'W87', 
 'lsc-aqwa-0m4a': 'W89', 
 'lsc-aqwb-0m4a': 'W79', 
 'cpt-doma-1m0a': 'K91', 
 'cpt-domb-1m0a': 'K92', 
 'cpt-domc-1m0a': 'K93', 
 'cpt-aqwa-0m4a': 'L09', 
 'coj-clma-0m4a': 'Q58', 
 'coj-clma-0m4b': 'Q59', 
 'coj-doma-1m0a': 'Q63', 
 'coj-domb-1m0a': 'Q64', 
 'coj-clma-2m0a': 'E10', 
 'tfn-aqwa-0m4a': 'Z21', 
 'tfn-aqwa-0m4b': 'Z17',
 'tfn-domb-1m0a': 'Z24', 
 'tfn-doma-1m0a': 'Z31',
 'ratir': 'ratir', 
 'liverpool': 'liverpool', 
 'zowada': 'zowada', 
 'WISE': 'WISE',
 'wise':'WISE', 
 'lijiang': 'lijiang', 
 'caha': 'CAHA', 
 'WMO': 'WMO', 
 'asiago': 'asiago', 
 'test': 'test',
 'CAHA': 'CAHA'}


main()
