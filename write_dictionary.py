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


json.dump(imagelist, open(date+"_input_dict.txt",'w'), indent=2)

