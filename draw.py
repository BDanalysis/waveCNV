import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import warnings
config = {
    "font.family": 'serif',
    "font.size": 15,
    "mathtext.fontset": 'stix',
    'axes.unicode_minus': False
}
from matplotlib import rcParams
rcParams.update(config)
warnings.simplefilter("ignore")

def drawAnswer(x, RD_raw,normal_segment,tumor_segment,tumor_CN):
    plt.scatter(x, RD_raw, s=2, c='black')
    for ind in range(len(normal_segment)):
        start, end, type = normal_segment[ind]
        CN = 2
        print(str(start) + '-' + str(end) + ':' + str(CN))
        index = np.array(range(start, end)).astype(int)
        value = RD_raw[index]
        plt.scatter(index, value, s=2, c='green')
    for ind in range(len(tumor_segment)):
        start, end, type = tumor_segment[ind]
        CN = tumor_CN[ind]
        print(str(start) + '-' + str(end) + ':' + str(CN))
        index = np.array(range(start, end)).astype(int)
        value = RD_raw[index]
        plt.scatter(index, value, s=2, c='red')
    plt.show()

def drawProfile(main_clust,show,tags,minor_clust):
    plt.subplot(121)
    plt.title("Sample A Score")
    for main in main_clust:
        index = np.argwhere(tags == main)
        value = show[index]
        p1 = plt.scatter(index, value, marker='o',c='gray',s=5,alpha=0.2,label='normal clusters')#colors[main % 10]
    for minor in minor_clust:
        index = np.argwhere(tags == minor)
        value = show[index]
        p2=plt.scatter(index, value, marker='o',c='black',s=5,alpha=0.8,label='outlier')
    plt.legend(handles=[p1,p2],labels=['normal clusters','outlier'], loc='best')
    plt.subplot(122)
    plt.title("Score profile")
    sns.distplot(show,bins=200,kde_kws={"color": "r", "lw": 1.5,"bw_adjust":0.25, "label": "KDE"},
                      hist_kws={"rwidth":0.5,"alpha": 1, "color": "k"}
                 )
    plt.show()

def drawRD(show,tags):
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    clust = np.unique(tags)
    for i in clust:
        index = np.argwhere(tags == i)
        value = show[index]
        plt.scatter(index, value, c=colors[i % 10], s=1, alpha=0.4)
    plt.show()

def drawA(main_clust,tags,show,minor_clust):
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    plt.title("Sample A")
    for main in main_clust:
        index = np.argwhere(tags == main)
        value = show[index]
        p1 = plt.scatter(index, value, marker='o',c='gray',s=5,alpha=0.2,label='normal clusters')#colors[main % 10]
    for minor in minor_clust:
        index = np.argwhere(tags == minor)
        value = show[index]
        p2=plt.scatter(index, value, marker='o',c='black',s=5,alpha=0.8,label='outlier')
    plt.legend(handles=[p1,p2],labels=['normal clusters','outlier'], loc='best')

    for main in main_clust:
        index = np.argwhere(tags == main)
        nums_ = index.shape[0]
        x = np.mean(index)
        value = np.mean(show[index])
        plt.scatter(x, value, c=colors[main % 10], marker='*',s=5*(nums_/100), alpha=1)
    for minor in minor_clust:
        index = np.argwhere(tags == minor)
        nums_ = index.shape[0]
        x = np.mean(index)
        value = np.mean(show[index])
        plt.scatter(x, value, c=colors[minor % 10],s=5*(nums_/20), alpha=0.9)
    plt.show()

def drawB(aim,length):
    plt.subplot(212)
    s = []
    aim= aim.reshape([-1])-np.min(aim)
    aim = aim/length*1000
    fq = np.zeros([1001])
    for i in range(1001):
        step = 1
        aim -= step
        index = np.argwhere(aim>0).reshape([-1])
        num = aim.shape[0]-index.shape[0]
        s.append(num)
        aim =aim[index]
    plt.plot(s,lw = 0.1,marker='o',markersize=2)
    plt.fill_betweenx(s, 0, range(1001),alpha=0.5)
    plt.show()

