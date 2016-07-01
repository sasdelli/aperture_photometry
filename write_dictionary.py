from match_stars import *
import json

# Date within the fits filenames.
date ='20120905'

# RA and DEC of the SN.
SN_RA = 96.7124 
SN_dec = 59.0840

# This needs to be the typical FWHM of the PSF of the stars in the image. An appropriate value can be estimated using the bultin functions of ds9. These are sensible initial values.
fwhmpsf=3
sigma=5
dannulus=10

# input fits file format
file_format = 'h_e_'+date+'*.fits'

#############

make_clean()
imagelist={}
for i in glob.glob(file_format):
    imagelist[i]={}
    imagelist[i]['sigma']=sigma
    imagelist[i]['fwhmpsf']=fwhmpsf
    imagelist[i]['dannulus']=dannulus
write_file_list(date)

input_dict={}
input_dict['images']=imagelist
input_dict['SN_RA']=SN_RA
input_dict['SN_DEC']=SN_dec
input_dict['date']=date

json.dump(input_dict, open("input_dict.txt",'w'), indent=2)

