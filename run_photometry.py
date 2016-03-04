from match_stars import *

file_1 = 'h_e_'
date ='20120905'
file_2 = '_61_1_1_1.fits'
SN_RA = 96.7124 
SN_dec = 59.0840

fwhmpsf=7

make_clean()
write_file_list(date)
fix_headers(date)
iraf.daophot()
iraf.daophot.daofind(image='@imagelist', output='default', fwhmpsf=fwhmpsf, sigma=15,
                     datamin=5, datamax=500000.,  interactive='no', verify='no')

iraf.daophot.phot(image='@imagelist', output='image.mag', coords='*coo*',
                  cbox=6, annulus=4*fwhmpsf, dannulus=10, apertures=16, sigma=5,
                  datamin=50., datamax=50000, verify='no')

mkshortnames(date)
iraf.mknobsfile(photfiles='image.mag', idfilters='G,R,I,Z,U',
                imsets='shortnames', observations='obsout')

file_name = file_1+date+file_2
match_Stars(fits_file=file_name, save_location='./', nobsfile='./obsout',
            standard_stars_file='standard_stars.dat',
            out_nobsfile='final_obsout_'+file_name, SN_coord=[SN_RA, SN_dec], dist_treshold=0.0003)

iraf.fitparams(observations='final_obsout_'+file_name,
               catalogs='standard_stars.dat', config='conf2',
               parameters='fitparams_output', interactive='no')

iraf.invertfit(observations='final_obsout_'+file_name, config='conf2',
               parameters='fitparams_output', calib='final_phot_'+date+'_'+str(fwhmpsf))


