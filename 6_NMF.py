# -*- coding: utf-8 -*-
"""
Created on Thu Mar 17 00:30:41 2022

@author: henry
"""

import pickle
import numpy as np
from sklearn.decomposition import NMF
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize
import re
import gzip
import sys

def smooth_contacts(contacts, frame_step, target_window_size):
    # frame step in nanoseconds
    # target window size in nanoseconds (fastest process you want to observe)
    halfstride = int((target_window_size/frame_step)/2)
    actual_window_size = (halfstride*2)*frame_step
    print("using stride: ", stride)
    print("which is a window of: ", actual_window_size, " nanoseconds")
    contact_averages = []
#     for j in range(int(contacts.shape[1]/stride)):
#         contact_averages.append(contacts[:,j*stride:(j+1)*stride,:].mean(axis=1))
    for j in range(int(contacts.shape[0])):
        windows_start = np.max([0,j-halfstride])
        windows_end = np.min([contacts.shape[0],j+halfstride])
        contact_averages.append(contacts[windows_start:windows_end,:].mean(axis=0))
    return np.asarray(contact_averages)

def find_longest_phase(W,threshold = 0.5):
    W_up = W>(np.max(W)*threshold)
    phaselength = []
    phasestart = []
    in_phase = 0
    phase_l = 0
    for i in range(len(W)):
        if in_phase:
            if W_up[i]:
                phase_l = phase_l+1
            else:
                phaselength.append(phase_l)
                in_phase = 0
                phase_l = 0
        else:
            if W_up[i]:
                phasestart.append(i)
                in_phase = 1
    if in_phase:
        phaselength.append(phase_l)
    phaselength = np.asarray(phaselength)
    phasestart = phasestart[np.argmax(phaselength)]
    phaselength = np.max(phaselength)
    return phaselength,phasestart

def load_data_in_range(parmname,range_to_load,strided,smooth):
    data = []
    data_cov = []
    smooth = int((smooth-1)/strided)+1
    cov1 = np.ones(smooth)/smooth
    for i in range(len(range_to_load)):
        traj = range_to_load[i][0]
        startf = int(range_to_load[i][1]/strided)
        endf = int(range_to_load[i][2]/strided)
        d = np.loadtxt('../analysis/%02d/analysis/parameters/%s.txt' %(traj,parmname))[::stride]
        data.append(d[startf:endf])
        d2 = np.convolve(d[startf:endf],cov1,mode='valid')
        data_cov.append(d2)
        
    return data,data_cov

def find_constructive_contact(H, search_window=0.05):
    N_comp = H.shape[0]
    cluster_list = []
    for i in range(N_comp):
        cluster_id,cluster_height = find_constructive_contact_single(H[i], search_window)
        cluster_list.append(cluster_id)
    cluster_list = np.array(cluster_list)
    return cluster_list

def find_constructive_contact_single(H_comp, search_window=0.05):
    top_H = np.max(H_comp)
    windows = int(top_H/(search_window/10))-10
    count_list = np.zeros(windows)
    for i in range(windows):
        height_i = top_H - i*search_window/10
        top_i = height_i+search_window/2
        down_i = height_i-search_window/2
        count_list[i] = np.sum((H_comp<top_i) & (H_comp>down_i))
    max_i = top_H - np.argmax(count_list)*search_window/10
    cluster_id = (H_comp<(max_i+search_window/2)) & (H_comp>(max_i-search_window/2))
    cluster_height = np.median(H_comp[cluster_id])
    return cluster_id,cluster_height

def find_real_x(x_to_find,x_start_in_traj,real_x_start=0,stride=5):
    real_x_array = ((np.array([x_to_find])-x_start_in_traj)/stride+real_x_start).astype(int)
    return real_x_array[0]


def stack_temporal_comp(W_comp,order,traj_id,rename,iflegend=0):
    i=traj_id
    stacked_W = W_comp*0
    stacked_W[:,0] = W_comp[:,order[0]]
    for j in range(1,len(order)):
        stacked_W[:,j] = stacked_W[:,j-1] + W_comp[:,order[j]]
    gs = gridspec.GridSpec(3, 3,figure=plt.figure(figsize=(6,4.5)))
    ax1 = plt.subplot(gs[:1, :])
    startf = np.max([0,int((range_list[i][0]-to_load[i][1]-int((251-1)/2))/stride2)])
    endf = np.min([int((range_list[i][1]-to_load[i][1]-int((251-1)/2))/stride2),len(Ser_Dist_cov[0])])
    xrange2 = np.arange(range_list[i][0]+to_load[i][1]-np.min([0,int((range_list[i][0]-to_load[i][1]-int((251-1)/2))/stride2)])*stride,range_list[i][1]+to_load[i][1]-np.max([int((range_list[i][1]-to_load[i][1]-int((251-1)/2))/stride2)-len(Ser_Dist_cov[0]),0])*stride,stride2)*0.08
    ax1.plot(xrange2,Ser_Dist_cov[i][startf:endf],alpha=0.7,color="g")
    ax1.plot(xrange2,B1H4_Dist_cov[i][startf:endf],alpha=0.7,color="b")
    #ax1.set_ylabel("CHL-Site\nDist (A)",fontsize=16)
    ax1.set_ylim([4.5,12])
    #ax1.set_title("Cholesterol binding mode in traj%02d"%(to_load[i][0]),fontsize=16)
    ax1.set_title("Cholesterol binding mode in RT%d"%(i+1),fontsize=16)
    ax1.set_xticks([])
    ax1.legend(("Ser136/147","β1-H4"),loc='upper left',fontsize=14)
    ax2 = plt.subplot(gs[1:, :])
    linewidth = 2 
    xrange = np.arange(range_list[i][0]+to_load[i][1],range_list[i][1]+to_load[i][1],stride)*0.08
    for j in range(len(order)):
        k = len(order) - 1 - j
        component = order[k]
        ax2.fill_between(xrange,stacked_W[:,k], label="comp " + rename[k],linewidth=linewidth,alpha=1,color = colorlist[k])
    #ax2.set_title("Temporal Weights of Components over Traj %d " %(to_load[i][0]),fontsize=16)
    ax2.set_title("Temporal Weights of Components over Traj %d " %(i+1),fontsize=16)
    ax2.set_ylabel("component weight",fontsize=16)
    ax2.set_xlabel("time (ns)",fontsize=16)
    #ax2.plot([22900*0.08,22900*0.08],[0,1],'k--')
    ax2.plot([23000*0.08,23000*0.08],[0,1],'k--')
    #ax2.plot([26000*0.08,26000*0.08],[0,1],'k--')
    #ax2.plot([27800*0.08,27800*0.08],[0,1],'k--')
    ax2.plot([28600*0.08,28600*0.08],[0,1],'k--')
    ax2.plot([30900*0.08,30900*0.08],[0,1],'k--')
    handles, labels = plt.gca().get_legend_handles_labels()
    #ax2.legend([handles[idx] for idx in np.arange(len(rename)-1,-1,-1)],[labels[idx] for idx in np.arange(len(rename)-1,-1,-1)],loc='upper left',bbox_to_anchor=(1,-0.5),fontsize=14,ncol = 2)
    if iflegend: 
        rename = rename[:len(order)]
        ax2.legend([handles[idx] for idx in len(rename)-1-np.argsort(rename)],[labels[idx] for idx in len(rename)-1-np.argsort(rename)],loc='upper left',bbox_to_anchor=(1,-0.5),fontsize=14,ncol = 2)
    plt.tight_layout()
    plt.show()
    plt.close()

def plot_temporal_comp(W_comp,order,traj_id,rename,iflegend=0):
    i=traj_id
    gs = gridspec.GridSpec(3, 3,figure=plt.figure(figsize=(6,4.5)))
    ax1 = plt.subplot(gs[:1, :])
    startf = np.max([0,int((range_list[i][0]-to_load[i][1]-int((251-1)/2))/stride2)])
    endf = np.min([int((range_list[i][1]-to_load[i][1]-int((251-1)/2))/stride2),len(Ser_Dist_cov[0])])
    xrange2 = np.arange(range_list[i][0]+to_load[i][1]-np.min([0,int((range_list[i][0]-to_load[i][1]-int((251-1)/2))/stride2)])*stride,range_list[i][1]+to_load[i][1]-np.max([int((range_list[i][1]-to_load[i][1]-int((251-1)/2))/stride2)-len(Ser_Dist_cov[0]),0])*stride,stride2)*0.08
    ax1.plot(xrange2,Ser_Dist_cov[i][startf:endf],alpha=0.7,color="g")
    ax1.plot(xrange2,B1H4_Dist_cov[i][startf:endf],alpha=0.7,color="b")
    #ax1.set_ylabel("CHL-Site\nDist (A)",fontsize=16)
    ax1.set_ylim([4.5,12])
    ax1.set_title("Cholesterol binding mode in traj%02d"%(to_load[i][0]),fontsize=16)
    ax1.set_xticks([])
    ax1.legend(("Ser136/147","β1-H4"),loc='upper left',fontsize=14)
    ax2 = plt.subplot(gs[1:, :])
    linewidth = 3 
    xrange = np.arange(range_list[i][0]+to_load[i][1],range_list[i][1]+to_load[i][1],stride)*0.08
    for j in range(len(order)):
        k = len(order) - 1 - j
        component = order[k]
        ax2.plot(xrange,W_comp[:,component], label="comp " + rename[k],linewidth=linewidth,alpha=1,color = colorlist[k])
    ax2.set_title("Temporal Weights of Components over Traj %d " %(to_load[i][0]),fontsize=16)
    ax2.set_ylabel("component weight",fontsize=16)
    ax2.set_xlabel("time (ns)",fontsize=16)
    #ax2.plot([22900*0.08,22900*0.08],[0,1],'k--')
    ax2.plot([23000*0.08,23000*0.08],[0,1],'k--')
    #ax2.plot([26000*0.08,26000*0.08],[0,1],'k--')
    #ax2.plot([27800*0.08,27800*0.08],[0,1],'k--')
    #ax2.plot([29800*0.08,29800*0.08],[0,1],'k--')
    ax2.plot([30900*0.08,30900*0.08],[0,1],'k--')
    handles, labels = plt.gca().get_legend_handles_labels()
    #ax2.legend([handles[idx] for idx in np.arange(len(rename)-1,-1,-1)],[labels[idx] for idx in np.arange(len(rename)-1,-1,-1)],loc='upper left',bbox_to_anchor=(1,-0.5),fontsize=14,ncol = 3)
    if iflegend: 
        rename = rename[:len(order)]
        ax2.legend([handles[idx] for idx in len(rename)-1-np.argsort(rename)],[labels[idx] for idx in len(rename)-1-np.argsort(rename)],loc='upper left',bbox_to_anchor=(1,-0.5),fontsize=14,ncol = 1)
    plt.tight_layout()
    plt.show()
    plt.close()


class MidpointNormalize(Normalize):
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return np.ma.masked_array(np.interp(value, x, y))

def draw_contactmap_from_H(this_H,v_name2,up_H_id,down_H_id,colormap):
    v_map = []
    for i in v_name2:
        res1 = int(re.sub('\D','',i.split('-')[0]))
        res2 = int(re.sub('\D','',i.split('-')[1]))
        v_map.append([res1,res2])
    v_map = np.asarray(v_map)
    minres = np.min(v_map)
    maxres = np.max(v_map)
    contactmap = np.zeros((maxres-minres+1,maxres-minres+1))*np.nan
    for i in range(v_map.shape[0]):
        contactmap[v_map[i][0]-minres,v_map[i][1]-minres] = this_H[i]
        contactmap[v_map[i][1]-minres,v_map[i][0]-minres] = this_H[i]
 

''' make variable name
v_name = np.loadtxt("check/variable_name.txt",dtype="bytes").astype(str)
sequence_s = "ASISTKLQNTLIQYHSIKEDEWRVAKKVKDVTVWRKPSEEFNGYLYKAQGVMDDVVNNVIDHIRPGPWRLDWDRLMTSLDVLEHFEENCCVMRYTTAGQLLNIISPREFVDFSYTVGYEEGLLSCGVSVEWSETRPEFVRGYNHPCGWFCVPLKDSPSQSLLTGYIQTDLRGMIPQSAVDTAMASTLANFYSDLRKGLRKA"
resid_list = np.array([24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224])
v_name2 = np.zeros(v_name.shape).astype("str")
for i in range(len(v_name)):
    resid1 = int(v_name[i].split("-")[0])
    resid2 = int(v_name[i].split("-")[1])
    if resid1==resid2:
        resname1 = sequence_s[resid1-24]
        v_name2[i] = "PIP2-"+resname1+str(resid1)
    else:
        resname1 = sequence_s[resid1-24]
        resname2 = sequence_s[resid2-24]
        v_name2[i] = "PIP2-"+resname1+str(resid1)+"-"+resname2+str(resid2)
np.savetxt("check/variable_name2.txt",v_name2,fmt="%s")
'''
 
''' load protein map '''
v_name_1 = np.loadtxt("../RED/check/variable_name2.txt",dtype="bytes").astype(str)
v_name_2 = np.loadtxt("../REDlipid/check/variable_name2.txt",dtype="bytes").astype(str)
v_name = np.concatenate((v_name_1,v_name_2))
sequence_s = "ASISTKLQNTLIQYHSIKEDEWRVAKKVKDVTVWRKPSEEFNGYLYKAQGVMDDVVNNVIDHIRPGPWRLDWDRLMTSLDVLEHFEENCCVMRYTTAGQLLNIISPREFVDFSYTVGYEEGLLSCGVSVEWSETRPEFVRGYNHPCGWFCVPLKDSPSQSLLTGYIQTDLRGMIPQSAVDTAMASTLANFYSDLRKGLRKA"
resid_list = np.array([24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224])
structure_list = np.array("C H H H H H H H H H H H H H H C C C G G G C E E E E E E T T E E E E E E E T T T T T T E E E E E E E E E C C C H H H H H H H H C C T T T T T T T T T T E E E E E E E E E E E T T E E E E E E E E C C B T T T T B T C E E E E E E E E E E E E T T E E E E E E E E C C C C C C T T T T E E C E E C C E E E E E E E E T T E E E E E E E E E E E C E E C C C C C C H H H H H H H H H H H H H H H H H H H H H H H C C C".split(" "))
motif_region=[[25,38],[46,51],[54,60],[67,75],[79,86],[90,97],[99,109],[112,119],[120,129],[130,141],[144,151],[162,176],[179,192],[199,221]]
motif_name=["H1","β1","β2","β3","H2","Hlp3","β4","β5","Ω1","β6","β7","β8","β9","H4"]
    
''' load trajectories of interest '''
to_load = [[7,0,50000]]
stride = 5
n_traj = len(to_load)
contactdata=[]
for i in range(n_traj):
    traj = to_load[i][0]
    startf = int(to_load[i][1]/stride)
    endf = int(to_load[i][2]/stride)
    print("traj %d"%traj)
    d = np.loadtxt("../RED/data/trimed_CHLp2b2_traj%04d.dat.gz"%to_load[i][0],dtype="bool")
    d1 = d[startf:endf,:]
    d = np.loadtxt("../REDlipid/data/trimed_CHLp2b2_traj%04d.dat.gz"%to_load[i][0],dtype="bool")
    d2 = d[startf:endf,:]
    d = np.concatenate((d1,d2),axis=1)
    contactdata.append(d)
    del d,d1,d2
contactdata_sm = []
for i in range(n_traj):
    contactdata_sm.append(smooth_contacts(contactdata[i], frame_step=0.08*stride, target_window_size=20))

''' load CV for comparison '''
parm1_1 = 'CHL-topSERs-mindist'
parm1_2 = 'Beta1-Cahelix-t7red-scdist'
stride2 = 5
Ser_Dist,Ser_Dist_cov = load_data_in_range(parm1_1,to_load,stride2,251)
B1H4_Dist,B1H4_Dist_cov = load_data_in_range(parm1_2,to_load,stride2,251)

''' define the trajectory fragments of interest '''
colorlist = [u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']
#range_list = [[0000,50000],[0000,50000],[0000,50000],[0000,50000],[0000,50000],[0000,50000]]
range_list = [[0000,50000]]
#range_list = [[12500,35000]]
cat_data = np.zeros((0,contactdata[0].shape[1]))
range_list_add = [0]
for i in range(len(contactdata)):
    startf = int((range_list[i][0]-to_load[i][1])/stride)
    endf = int((range_list[i][1]-to_load[i][1])/stride)
    range_list_add.append(range_list_add[-1]+endf-startf)
    cat_data = np.concatenate((cat_data,contactdata_sm[i][startf:endf,:]))

''' NMF '''
N_component = 6
model = NMF(n_components=N_component, init='nndsvd', random_state=0, l1_ratio=0, max_iter=3000, alpha=1)
W = model.fit_transform(cat_data)
H = model.components_

''' plot concatenated traj '''
#plt.figure(figsize=(6,5))
plt.figure(figsize=(7.8,5))
linewidth = 2 
for component in range(N_component):
    plt.plot(W[:,component], label="comp " + str(component),linewidth=linewidth)
plt.title("Weights of Components over CHL-PI35P2 p2b2 Traj%02d" %(traj))
plt.ylabel("component weight")
plt.xlabel("frame")
x_ticks = (np.array([0,0.25,0.5,0.75,1])*(range_list[0][1]-range_list[0][0])+range_list[0][0]).astype(int)
plt.xticks(find_real_x(x_ticks,range_list[0][0]),find_real_x(x_ticks,range_list[0][0]).astype(int))
#plt.plot([find_real_x(22900,range_list[0][0]),find_real_x(22900,range_list[0][0])],[0,np.max(W)],'k--')
plt.plot([find_real_x(23900,range_list[0][0]),find_real_x(23900,range_list[0][0])],[0,np.max(W)],'k--')
#plt.plot([find_real_x(26000,range_list[0][0]),find_real_x(26000,range_list[0][0])],[0,np.max(W)],'k--')
#plt.plot([find_real_x(29800,range_list[0][0]),find_real_x(29800,range_list[0][0])],[0,np.max(W)],'k--')
plt.plot([find_real_x(30900,range_list[0][0]),find_real_x(30900,range_list[0][0])],[0,np.max(W)],'k--')
#plt.plot([10000,10000],[0,np.max(W)*1.1],'c--')
#plt.plot([20000,20000],[0,np.max(W)*1.1],'c--')
#plt.plot([30000,30000],[0,np.max(W)*1.1],'c--')
#plt.plot([40000,40000],[0,np.max(W)*1.1],'c--')
#plt.plot([50000,50000],[0,np.max(W)*1.1],'c--')
plt.legend(loc='center left', bbox_to_anchor=(1.05,0.7),fontsize=14)
#plt.legend(loc='upper right',fontsize=10,ncol = 3)
plt.tight_layout()
plt.show()
plt.close()

'''  plot every H  '''
thres_fold = 3
v_list = []
for comp in range(7):
    thres = np.mean(H[comp])+np.std(H[comp])*thres_fold
    v_contribute = v_name[H[comp]>thres]
    plt.plot(H[comp],marker='.',linestyle='',color=colorlist[comp],markersize=5)
    plt.ylabel("contribution",fontsize=14)
    plt.plot([-10,len(v_name)+10],[thres,thres],"k--")
    plt.xlim([-10,len(v_name)+10])
    plt.title("Variable contribution in component %d" %comp)
    plt.show()
    plt.close()
    v_list.append(v_contribute)

'''  rewight based on constructive contact '''
constructive_pairs = np.sum(find_constructive_contact(H),axis=0)==N_component
H_new = H*1
W_new = W*1
for i in range(N_component):
    reweight = np.mean(H[i][constructive_pairs])
    H_new[i] = H[i]/reweight
    W_new[:,i] = W[:,i]*reweight
H_old = H*1
W_old = W*1
H = H_new*1
W = W_new*1

''' plot temporal weight '''
traj_id=0
W_test = W_old[range_list_add[traj_id]:range_list_add[traj_id+1],:]
#dominance = np.argmax(W_test,axis=1)
#peak = np.argmax(W_test,axis=0)
#comp_order = np.argsort(np.argmax(W_test,axis=0))
#stack_temporal_comp(W_test,comp_order,traj_id)
comp_order = [0,3,6,1,4,2,5]
rename = ["4","3","2","6","0","1","5"]
colorlist = plt.cm.tab20(list(np.array([0,6,2,1,4,5,3])/20))
plot_temporal_comp(W_test,comp_order,traj_id,rename,0)


'''  stacked temporal weight '''
'''  the rename and reorder is unnecessary, only for a nice sequential naming to publish '''
traj_id=0
W_test = W[range_list_add[traj_id]:range_list_add[traj_id+1],:]
#dominance = np.argmax(W_test,axis=1)
#peak = np.argmax(W_test,axis=0)
#comp_order = np.argsort(np.argmax(W_test,axis=0))
#stack_temporal_comp(W_test,comp_order,traj_id)
#rename = ["4","3","2","6","0","1","5"]
rename = ["3","2","0","5","1","4"] 
comp_order = [5,2,4,1,0,3]
colorlist = plt.cm.tab20(list(np.array([0,6,2,1,4,5,3])/20))
stack_temporal_comp(W_test,comp_order,traj_id,rename,0)

traj_id=0
W_test = W[range_list_add[traj_id]:range_list_add[traj_id+1],:]
#rename = ["B","D","C","F","E","A"]
rename = ["0","1","2","3","4","5"]
comp_order = [0,1,2,3,4,5]
colorlist = plt.cm.tab20(list((np.array([0,6,2,1,4,5,3])+10)/20))
stack_temporal_comp(W_test,comp_order,traj_id,rename,0)



''' H difference '''
plt.figure(figsize=(6,4))
comp_before_name = "0"
comp_after_name = "3"
comp_before = np.array(comp_order)[np.array(rename)==comp_before_name][0]
comp_after = np.array(comp_order)[np.array(rename)==comp_after_name][0]
diff_H = H[comp_after] - H[comp_before]
up_H_id = np.where(diff_H>(+np.std(diff_H)*4))[0]
down_H_id = np.where(diff_H<(-np.std(diff_H)*4))[0]
plt.plot(diff_H,marker='.',linestyle='',color="gray",markersize=3)
plt.plot(up_H_id,diff_H[up_H_id],marker='.',linestyle='',color=colorlist[np.where(np.array(comp_order)==comp_after)[0][0]],markersize=5, label="comp " + comp_after_name)
plt.plot(down_H_id,diff_H[down_H_id],marker='.',linestyle='',color=colorlist[np.where(np.array(comp_order)==comp_before)[0][0]],markersize=5, label="comp " + comp_before_name)
plt.title("transition from component %s to %s" %(comp_before_name,comp_after_name),fontsize=18)
plt.legend(loc='lower left',fontsize=14,ncol = 3)
plt.show()
plt.close()
eventAC_form = up_H_id*1
eventAC_break = down_H_id*1
after_transit = v_name[up_H_id]
before_transit = v_name[down_H_id]
print("Forming contact pairs (%d pairs)"%len(after_transit))
print(", ".join(after_transit))
print("Breaking contact pairs (%d pairs)"%len(before_transit))
print(", ".join(before_transit))




