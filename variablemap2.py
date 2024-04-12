import os,sys,glob
import numpy as np
import gzip

systems = sys.argv[1]
systems = systems.split(",")
startfile = int(sys.argv[2])
endfile = int(sys.argv[3])

rez_check = np.load('check/rez_check.npy')
rez_num = len(rez_check)

useful = np.zeros((rez_num,rez_num)).astype("bool")

filerange = range(startfile,endfile+1)
for system in systems:
  for i in filerange:
    data = np.load("check/variable_map_%s%04d.npy"%(system,i)).astype("bool")
    useful = useful + data

valiable_id = np.where(useful==1)
f = open('check/variable_name.txt', "w")
for i in range(len(valiable_id[0])):
  name1 = rez_check[valiable_id[0][i]]
  name2 = rez_check[valiable_id[1][i]]
  f.write("%d-%d\n" %(name1,name2))
f.close()

positive_rate = np.sum(useful)/rez_num/rez_num
print("ratio of variating data: %0.4f%%" %(positive_rate*100))
print(useful.shape)

np.save('check/variable_map.npy', useful)


