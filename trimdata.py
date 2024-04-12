import os,sys,glob
import numpy as np
import gzip

system = sys.argv[1]
i = int(sys.argv[2])
furtherstride = int(sys.argv[3])

rez_check = np.load('check/rez_check.npy')
rez_num = len(rez_check)

variable_map = np.load('check/variable_map.npy')
valiable_id = np.where(variable_map.flatten()==1)[0]
readf = gzip.GzipFile("data/%s_traj%04d.dat.gz" %(system,i),'r')
data = np.loadtxt(readf,dtype="bool")
readf.close()
length = data.shape[0]
stride = np.arange(0,length,furtherstride)
datatrimed = data[:,valiable_id]
datatrimed2 = datatrimed[stride,:]
savef = gzip.GzipFile("data/trimed_%s_traj%04d.dat.gz" %(system,i),'w')
np.savetxt(savef, datatrimed2, fmt="%d")
savef.close()
 
