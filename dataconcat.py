import os,sys,glob
import numpy as np
import gzip

startfile = int(sys.argv[1])
endfile = int(sys.argv[2])
system = sys.argv[3]
trajid = int(sys.argv[4])
filelength = 1000
stride = 5

rez_check = np.load('check/rez_check.npy')
rez_num = len(rez_check)

if os.path.exists("data/%s_traj%04d.dat.gz" %(system,trajid)):
    savefile = gzip.GzipFile("data/%s_traj%04d.dat.gz" %(system,trajid),'r')
    previous_data = np.loadtxt(savefile,dtype="bool")
    savefile.close()
    previouslen = previous_data.shape[0]
    fulldata = np.zeros((int(previouslen+(endfile-startfile+1)*filelength/stride),rez_num*rez_num),dtype="bool")
    fulldata[:previouslen,:] = previous_data
else:
    previouslen = 0
    fulldata = np.zeros((int(previouslen+(endfile-startfile+1)*filelength/stride),rez_num*rez_num),dtype="bool")

filerange = range(startfile,endfile+1)
for i in filerange:
    f = gzip.GzipFile("data/%s_traj%04d/%s_traj%04d_%06d.gz" %(system,trajid,system,trajid,(i*filelength)),'r')
    data = np.loadtxt(f,dtype="bool")
    f.close()
    print(data.shape)
    fulldata[int(previouslen+i*filelength/stride):int(previouslen+(i+1)*filelength/stride),:] = data
    positive_rate = np.sum(data)/data.flatten().shape
    print("%s_traj%04d_%06d positive rate: %0.4f%%" %(system,trajid,(i*filelength),(positive_rate*100)))

print(fulldata.shape)
positive_rate = np.sum(fulldata)/fulldata.flatten().shape
print("average positive rate: %0.4f%%" %(positive_rate*100))
print(fulldata.shape)
savefile = gzip.GzipFile("data/%s_traj%04d.dat.gz" %(system,trajid),'w')
np.savetxt(savefile, fulldata, fmt="%d")
savefile.close()

