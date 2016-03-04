from match_stars import *

date ='20120905'

SN_RA = 96.7124 
SN_dec = 59.0840

fwhmpsf=7

make_clean()
write_file_list(date)
fix_headers(date)
iraf.daophot()
run_sextractor_daofind('imagelist')
#iraf.daophot.daofind(image='@imagelist', output='default', fwhmpsf=fwhmpsf, sigma=15,
#                     datamin=5, datamax=500000.,  interactive='no', verify='no')

iraf.daophot.phot(image='@imagelist', output='image.mag', coords='*coo*',
                  cbox=6, annulus=4*fwhmpsf, dannulus=10, apertures=16, sigma=5,
                  datamin=20., datamax=50000, verify='no')

mkshortnames(date)
iraf.mknobsfile(photfiles='image.mag', idfilters='G,R,I,Z,U',
                imsets='shortnames', observations='obsout')
file_name = find_band(date, band='G')

match_Stars(fits_file=file_name, out_nobsfile='final_obsout_'+file_name,
            SN_coord=[SN_RA, SN_dec], dist_treshold=0.0003)

iraf.fitparams(observations='final_obsout_'+file_name,
               catalogs='standard_stars.dat', config='conf2',
               parameters='fitparams_output', interactive='no')

iraf.invertfit(observations='final_obsout_'+file_name, config='conf2',
               parameters='fitparams_output', calib='final_phot_'+date+'_'+str(fwhmpsf))


