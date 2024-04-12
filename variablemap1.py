import os,sys,glob
import numpy as np
import gzip

system = sys.argv[1]
traj = int(sys.argv[2])

rez_check = np.load('check/rez_check.npy')
rez_num = len(rez_check)

useful = np.zeros((rez_num,rez_num)).astype("bool")
variating = np.zeros((rez_num,rez_num)).astype("bool")

f = gzip.GzipFile("data/%s_traj%04d.dat.gz" %(system,traj),'r')
data = np.loadtxt(f,dtype="bool")
f.close()
variating = variating + (np.std(data,axis=0)>0).reshape(rez_num,rez_num)

for i in range(rez_num):
  for j in range(rez_num):
    if j>i+1:
      useful[i][j] = variating[i][j]

positive_rate = np.sum(useful)/rez_num/rez_num
variating_rate = np.sum(variating)/rez_num/rez_num
print("ratio of useful data: %0.4f%%" %(positive_rate*100))
print("ratio of variating data: %0.4f%%" %(variating_rate*100))
print(useful.shape)

np.save('check/variable_map_%s%04d.npy' %(system,traj), useful)

