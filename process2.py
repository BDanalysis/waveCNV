import warnings
import sys
import numpy as np
import mytime
import draw as D
np.set_printoptions(suppress=True)

from readTool import *
from waveCluster import *
from scipy.stats import rv_continuous, norm
warnings.simplefilter("ignore")

def calculate_answer(label):
    label = label.reshape(-1)
    index1 = np.argwhere(label==1)
    if len(index1)==0:
        answer1 = [[0,0]]
        return answer1
    begin = index1[0]
    answer1 = []
    for i in range(len(index1)):
        if i == len(index1)-1:
            if index1[i]-index1[i-1]<=1:
                answer1.append([int(begin), int(index1[i])])
            else:
                answer1.append([int(index1[i]), int(index1[i])])
        elif index1[i+1]-index1[i]<=5:
            continue
        else:
            answer1.append([int(begin),int(index1[i]+1)])
            begin = index1[i+1]
    return np.array(answer1)

def ensure_bad_pery(perp,RD):
    length = len(perp)
    pery_list = np.zeros(RD.shape[0])
    for i in range(length):
        start = perp[i][0]
        end = perp[i][1]
        pery_list[start:end+1]=1
    return pery_list

def ensure_bad_bin(truth_start,truth_end,RD,binSize):
    length = len(truth_start)
    truth_list = np.zeros(np.shape(RD))
    top = len(RD)
    for i in range(length):
        start = np.ceil(truth_start[i]/binSize).astype(int)
        end = np.ceil(truth_end[i]/binSize).astype(int)
        if start <= top and end <= top:
            truth_list[start:end]=1
    return truth_list

def gcCheck(RD,GC):
    RD = RD/1000
    RD = RD.reshape(-1)
    GC = GC.reshape(-1)
    #RD[RD == 0] = modeRD(RD)
    RD = gc_correct(RD, GC)
    return RD*1000

def merge_CNV_segment(RD,index,balance,gap):
    segment = [-1]
    start = -1
    for i in range(len(index)):
        if start == -1:
            start = index[i]
        elif index[i]-index[i-1]<=10 and i != len(index)-1:
            continue
        else:
            segment.append([start,index[i-1]])
            start = index[i]
    segment.pop(0)
    answer = []
    for start,end in segment:
        if end - start > gap:
            tumor_mean = np.mean(RD[start:end])
            if tumor_mean > balance:
                answer.append([start,end,1])
            elif balance == -1:
                answer.append([start, end, 0])
            else:
                answer.append([start, end, -1])
    return answer

def gc_correct(RD, GC):
    # correcting gc bias
    bincount = np.bincount(GC)
    global_rd_ave = np.mean(RD)
    for i in range(len(RD)):
        if bincount[GC[i]] < 2:
            continue
        mean = np.mean(RD[GC == GC[i]])
        if RD[i] != 0:
            RD[i] = global_rd_ave * RD[i] / mean
    #print(np.argwhere(np.isnan(RD)))
    return RD

def alignment(RD, GC, truth_list):
    index = np.argwhere(RD>=0)
    RD = RD[index]
    GC = GC[index]
    truth_list = truth_list[index]
    return RD,GC,truth_list,index

def normalization(data):
    means = np.nanmean(data)
    data = data/means
    return data.reshape(1,-1)

def add_position(data):
    length = len(data)
    x = np.arange(1, length+1).reshape(-1)
    x = x[0:length].reshape([-1,1])
    answer = np.c_[data, x]
    return answer

def calculating_CN(RD,tumor_segment):
    tumor_segment = np.array(tumor_segment)
    mode = np.mean(RD[normal_index])
    mapp_index = np.argwhere(tumor_segment[:,2]==-1)
    loss = tumor_segment[np.where(tumor_segment[:,2]==-1)]
    loss_mean = base_num = 0
    for start,end,type in loss:
        loss_mean += np.sum(RD[start:end+1])
        base_num += (end-start+1)
    loss_mean = loss_mean/base_num
    homoRD = base_num1 = 0
    hemiRD = base_num2 = 0
    for idx in range(len(loss)):
        start, end, type = loss[idx]
        tmp_mean = np.mean(RD[start:end+1])
        if tmp_mean > loss_mean:
            ss = int(mapp_index[idx])
            tumor_segment[ss,2] = -1
            hemiRD += np.sum(RD[start:end+1])
            base_num2 += (end-start+1)
        else:
            ss = int(mapp_index[idx])
            tumor_segment[ss, 2] = -2
            homoRD += np.sum(RD[start:end+1])
            base_num1 += (end - start + 1)
    if base_num1==0:
        homoRD = 0
    else:
        homoRD = homoRD / base_num1
    if base_num2 == 0:
        hemiRD = 1
    else:
        hemiRD = hemiRD / base_num2
    length = len(tumor_segment)
    CN = np.full(length, 0)
    purity = 2 * (homoRD - hemiRD) / (homoRD - 2 * hemiRD)
    print(purity)
    for i in range(len(tumor_segment)):
        begin,end,type = tumor_segment[i]
        error = np.mean(RD[begin:end])
        CN[i] = int(2*error / (mode * purity) - 2 * (1-purity) / purity)
    return CN,tumor_segment

def mapping_index(all_index,unnormal_RD,index_filter):
    max_length = np.max(all_index)
    RD_raw = np.zeros([max_length+1])
    index_raw = all_index[index_filter]
    index_raw = index_raw.reshape([-1])
    all_index = all_index.reshape([-1])
    RD_raw[all_index] = unnormal_RD
    return RD_raw,index_raw

refpath = sys.argv[1]
targe = sys.argv[5]
binSize = int(sys.argv[3])
bam = sys.argv[2]
groundTruth = "./GroundTruthCNV"#sys.argv[5]
outpath = sys.argv[4]
flag_draw = False

try:
    truth_start, truth_end = read_truth_file(groundTruth)
except IOError:
    truth_start = []
    truth_end = []
mid1 = bam.split("/")
mid = mid1[-1].split(".")[0]+mid1[-1].split(".")[1]
outfile = outpath + '/'+mid+".result.txt"
outsorcefile = outpath + '/' + mid +".sorce.txt"

ref = [[] for i in range(25)]
try:
    refList = read_bam_file(bam)
except ValueError:
    print("ValueError:"+bam+"may be the file is error")
    sys.exit(1)

for i in range(len(refList)):
    chr = refList[i]
    chr_num = chr.strip('chr')
    if chr_num != targe:
        continue
    reference = refpath + '/chr' + str(chr_num) + '.fa'
    if chr_num == 'X' or chr_num == 'x':
        chr_num = str(23)
    elif chr_num == 'Y' or chr_num == 'y':
        chr_num = str(24)
    if chr_num.isdigit():
        chr_num = int(chr_num)
        ref = read_ref_file(reference, chr_num, ref)
chrLen = np.full(25, 0)
for i in range(1, 25):
    chrLen[i] = len(ref[i])
RD_raw,GC,out_chr = Binning(ref, binSize, chrLen, bam)

truth_list = ensure_bad_bin(truth_start,truth_end,RD_raw,binSize)
RD, GC, truth_list,all_index = alignment(RD_raw, GC, truth_list)
unnormal_RD = gcCheck(RD,GC)
times = mytime.MyTimer()
times.start()
#plt.plot(unnormal_RD ,lw=0,marker='o',alpha=0.5,markersize=2)
#plt.show()
aim = np.zeros_like(unnormal_RD )
print(unnormal_RD .shape[0])
for i in range(unnormal_RD.shape[0]):
    tmp = np.mean(unnormal_RD [i:i+10])
    aim[i:i+10] = tmp
aim = np.log10(aim+1)


show = add_position(aim)
tags = waveCluster(show, scale=1024,wavelet="db2", threshold=-99, plot=False)#{'db1':0,'db2':1,'bior1.3':2}

# drow
if flag_draw:
    D.drawRD(aim,tags)
show = aim


mean_S = tags.shape[0] / np.unique(tags).shape[0]
means_main = np.zeros_like(np.unique(tags))
for c in range(means_main.shape[0]):
    index = np.argwhere(tags == c).reshape([-1])
    means_main[c] = index.shape[0]
means_main = -1*means_main
sort_clust = np.argsort(means_main)
means_main = np.sort(means_main)*-1
total = 0
i = 0
for i in range(len(means_main)):
    total += means_main[i]
    if total >= (0.9*tags.shape[0]):
        total = i
        break
main_clust = sort_clust[0:total+1]
means_main = []
minor_clust = sort_clust[total+1:]
if flag_draw:
    D.drawA(main_clust,tags,show,minor_clust)
for main in main_clust:
    index = np.argwhere(tags == main).reshape([-1])
    tmp = np.mean(show[index])
    tmp_x = np.mean(index)
    if tmp!=0:
        dis = np.sqrt(tmp**2+tmp_x**2)
        means_main.append([tmp,dis])
means_main= np.array(means_main)

for c in np.unique(tags):
    index = np.argwhere(tags==c).reshape([-1])
    tag = tags[index]
    tmp = np.mean(show[index])
    tmp_x = np.mean(index)
    dis = np.sqrt(tmp**2+tmp_x**2)
    ind = np.argmin(np.abs(means_main[:,1]-dis))

    show[index] =  show[index]/(means_main[ind,0])
if flag_draw:
    D.drawProfile(main_clust,show,tags,minor_clust)


length = (np.max(show)-np.min(show))
show_data = show[np.where(show!=0)]
mu,sigma = rv_continuous.fit(norm,show_data,loc=1)
print(mu,sigma)
results = (1 - norm.cdf(show, mu, sigma*0.9))
check = (results < 0.01) | (results > 0.99)

tumor_index = np.argwhere(check).reshape([-1])
normal_index =np.argwhere(~check).reshape([-1])
mode = np.mean(RD[normal_index])

RD_raw_after,tumor_index_raw = mapping_index(all_index,unnormal_RD,tumor_index)
RD_raw_after,normal_index_raw = mapping_index(all_index,unnormal_RD,normal_index)
tumor_segment = merge_CNV_segment(RD_raw_after,tumor_index_raw,mode,15)
normal_segment = merge_CNV_segment(RD_raw_after,normal_index_raw,-1,1)
tumor_CN,tumor_segment = calculating_CN(RD_raw_after,tumor_segment)

write_CNV_file(outfile,tumor_segment,tumor_CN,binSize)
x = range(len(RD_raw))
times.stop()
if flag_draw:
    D.drawAnswer(x, RD_raw,normal_segment,tumor_segment,tumor_CN)
    D.drawB(aim,length)