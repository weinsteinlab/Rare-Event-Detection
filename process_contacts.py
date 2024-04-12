import numpy as np
import sys
import gzip

infile = sys.argv[1]
outfile = sys.argv[2]

f = open(infile, "r")
savefile = gzip.GzipFile("%s.gz" %outfile,'w')

rez_check = np.arange(24,224+1)
np.save('check/rez_check.npy',rez_check)
rez_count = len(rez_check)

temp=[]
for frame in f:
    framemap = np.zeros((rez_count,rez_count),dtype="bool")
    y = frame.split(',')[:-1]
    for i in range(rez_count):
        rezs_contact = y[i].split(' ')
        if "" in rezs_contact: rezs_contact.remove("") 
        for residue in rezs_contact:
            residue = np.where(rez_check==int(residue))[0][0]
            framemap[i][residue] = 1
    temp.append(framemap)
f.close()

temp = np.asarray(temp)
print(temp.shape)
temp = temp.reshape(temp.shape[0],temp.shape[1]*temp.shape[2])
print(temp.shape)
np.savetxt(savefile, temp, fmt="%d")
savefile.close()
