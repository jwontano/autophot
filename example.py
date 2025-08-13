# John W. Montano - UCI
# mrk 817 extended campaign example
from auto_photometry import *
def photometry_update():
    # these are things you need to gather before you start using the program, and only needs to be done once
    # Right ascension and declination for your AGN that is the first pair, and the rest are comparison stars within the FOV of your images that have recorded magnitudes
    # not all need to have magnitudes but that will be needed for the flux conversion part
    coords = [ [219.09200836,58.7942765], [219.1357,58.83312], [219.17719,58.84183] , [219.150785,58.766423] , [219.077094,58.677452], [219.05994,58.666196]] 
    # what you call your target, this is what will show up in filenames and the like
    target_name = 'Mrk817'
    # the output directory, I make directories for every object I work and you can make subdirectories if you want to try out different things :]
    outdir = '/home/korbinite5646/AGN_home/MRK817/extended/'
    # you have to know what telescopes you are using and the label name you will assign to them
    # this is for plots generated so that we know what data points belong to what telescopes
    # in this example the STORM 2 project has 12 telescopes
    cdict={'liverpool': 'black',
    'zowada': 'blue',
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
    # observatories have minor planet center codes the program uses it for labels in plots and filenames
    # you can make up your own labels
    # https://minorplanetcenter.net/iau/lists/ObsCodesF.html
    MPC ={'ogg-clma-2m0a': 'F65', 
     'elp-doma-1m0a': 'V37', 
     'elp-domb-1m0a': 'V39', 
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
     'CAHA': 'CAHA'}
   # where your images are
   # directory = '/home/korbinite5646/AGN_home/FITS/KEY2023B-001/Mrk 817/'
  # the filters you are using 
  img_filters = ['g','r','i','z']
    
    
    directory = '/media/korbinite5646/backup/STORM2/FITS/zowada/z/defringed/'
    for img_filter in img_filters:
      # note the inputs for this function, this is the photometry if your pixelscale or gain is off your data will be too
      phot( directory,
    	outdir,
    	coords=coords,
    	img_filter=img_filter,
    	target_name=target_name,
    	telname='zowada',
    	snr_cutoff=0.,
    	renew=True,
    	MPC=MPC,
        rdnoise=8.0,
        pixsc=0.882,
        gain=0.476)
      # this function creates separate tables/files based on telescopes and object you'll see the output in your directory
    	phot_to_one(outdir,img_filter=img_filter,coords=coords)
      # the normalization process involving the comparison stars
    	fit(outdir,coords,img_filter)
      # normalizes the comparison star light curves
    	star_norm(outdir,coords,img_filter,clean=True)
      
    	plot_normed(outdir,img_filter,coords=coords,cdict=cdict,mpc=MPC)
    	# applies the normalization factors to the AGN light curve
      scale_AGN(outdir, img_filter=img_filter, target_name=target_name,clean=True)
    	plot_scaled(outdir,img_filter=img_filter,target_name=target_name,cdict=cdict)
    	# sanity check to make sure no bad data points exist
      clean_lc(outdir,img_filter)
      # simple alignment plot by way of normalizign to the mean of a reference light curve, chose the reference or it'll be random
    	plot_align(outdir=outdir,img_filter=img_filter,target_name=target_name,cdict=cdict,reference='V37')
    # flux conversion from data units/ electrons/s to f_lambda cgs units
    true_mag_retrieve(outdir,tar='Mrk817',f_to=img_filters)
    plot_cgs(cdict=cdict,whoo=target_name,outdir=outdir,f_to=img_filters)
    # will average points taken on the same night
    for img_filter in img_filters: condense(outdir,img_filter)
    for img_filter in img_filters: plot_condense(outdir,img_filter,target_name,cdict,MPC,single=False)
    # will format the light curves made from the previous functions to the format that pycali requires
    cali_format(outdir,img_filters,target_name,reference='V37')

photometry_update()


outdir = '/home/korbinite5646/AGN_home/MRK817/extended/'
# this is a very basic run of pycali, again need to have pycali installed and properly in your PATH
pycali_run('B',outdir=outdir)
pycali_run('V',outdir=outdir)
pycali_run('u',outdir=outdir)
pycali_run('g',outdir=outdir)
pycali_run('r',outdir=outdir)
pycali_run('i',outdir=outdir)
pycali_run('z',outdir=outdir)
# produce_clean(['g','r','i','z'],'/home/korbinite5646/AGN_home/MRK817/extended/','Mrk817')
# pycali_run('B',outdir=outdir,run=2)
# pycali_run('V',outdir=outdir,run=2)
# pycali_run('u',outdir=outdir,run=2)
# pycali_run('g',outdir=outdir,run=2)
# pycali_run('r',outdir=outdir,run=2)
# pycali_run('i',outdir=outdir,run=2)
# pycali_run('z',outdir=outdir,run=2)

