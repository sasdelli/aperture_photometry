from match_stars import *
import json

date ='20120905'

SN_RA = 96.7124 
SN_dec = 59.0840


fwhmpsf=3

make_clean()
imagelist={}
for i in glob.glob('h_e_'+date+'*.fits'):
    imagelist[i]={}
    imagelist[i]['sigma']=5
    imagelist[i]['fwhmpsf']=3
    imagelist[i]['dannulus']=10
write_file_list(date)

input_dict={}
input_dict['images']=imagelist
input_dict['SN_RA']=SN_RA
input_dict['SN_DEC']=SN_dec
input_dict['date']=date

json.dump(input_dict, open("input_dict.txt",'w'), indent=2)

