import json
from match_stars import *

# Clean old temporary files.
make_clean()

# Load input dictionary.
input_dict = json.load(open("input_dict.txt"))
date = input_dict['date']
SN_RA = input_dict['SN_RA']
SN_dec = input_dict['SN_DEC']
imagelist = input_dict['images']


# Write imagelist for daofind or sextractor
write_file_list(date)

# Changes the name of the filter in the header of the fits files. Necessary for LT fits files. 
fix_headers(date)

iraf.daophot()

# Run sextractor to find the sources in the images.
run_sextractor_daofind('imagelist')

##iraf.daophot.daofind(image='@imagelist', output='default', fwhmpsf=fwhmpsf, sigma=15,
##                     datamin=5, datamax=500000.,  interactive='no', verify='no')

##iraf.daophot.phot(image='@imagelist', output='image.mag', coords='*coo*',
##                  cbox=6, annulus=4*fwhmpsf, dannulus=10, apertures=16, sigma=5,
##                  datamin=20., datamax=50000, verify='no')

# Run the daophot phot rutine on the images, creating the obs files with the magnitudes.
for imagename in imagelist.keys():
    image = imagelist[imagename]
    iraf.daophot.phot(image=imagename, output=imagename+'.mag.1',
                      coords=imagename+'.coo.1', cbox=6,
                      annulus=4*image['fwhmpsf'], dannulus=image['dannulus'],
                      apertures=16, sigma=image['sigma'], datamin=20.,
                      datamax=50000, verify='no')

# Create the "shortnames" file with shortened filenames for mknobsfile
mkshortnames(date)
##iraf.mknobsfile(photfiles='image.mag', idfilters='G,R,I,Z,U',
##                imsets='shortnames', observations='obsout')
iraf.mknobsfile(photfiles='*mag.1*', idfilters='G,R,I,Z,U',
                imsets='shortnames', observations='obsout')
file_name = find_band(date, band='G')

# Match the stars in the fits files with the catalog stars.
match_Stars(fits_file=file_name, out_nobsfile='final_obsout_'+file_name,
            SN_coord=[SN_RA, SN_dec], dist_treshold=0.0003, header_line=1)

# Fit the photometries.
iraf.fitparams(observations='final_obsout_'+file_name,
               catalogs='standard_stars.dat', config='conf2',
               parameters='fitparams_output', interactive='no')
iraf.invertfit(observations='final_obsout_'+file_name, config='conf2',
               parameters='fitparams_output', calib='final_phot_'+date ) #+'_'+str(fwhmpsf))



