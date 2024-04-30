from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-t","--tag"    , help="tag. output files will be saved as <tag>.links , <tag>.png" ,default="hf_raw")
parser.add_argument("-s","--seed"   , help="Seed for the random number generator",type=int , default=799 )
args = parser.parse_args()

np.random.seed(args.seed)


def fillCluster(HF,Energy=30.0,width=0.8,verbose=10):

    vals=stats.norm.pdf(np.array([-2,-1,0,1,2])/width)/stats.norm.pdf(0)
    en_dist=np.outer(vals,vals)

    selct_=np.random.random(25).reshape(5,5)
    mask = selct_> en_dist
    en_fraction=np.array(en_dist)
    en_fraction[mask]=0.0
    en_fraction=en_fraction/np.sum(en_fraction)
    en_fraction=(en_fraction + en_fraction*np.random.normal(scale=0.05) )*Energy
    x0,y0=np.random.randint(11),np.random.randint(18*2)
    if verbose >0:
        print("Adding cluster at ",x0,y0," of energy ",Energy)
    if verbose >2:
        print("en_fraction\n",en_fraction)
    for i in range(-2,2+1):
        for j in range(-2,2+1):
    #         print(i,j,x[2+i][2+j],en_fraction[i+2][j+2])
            xval=x0+i
            yval=y0+j
            if xval >= 11:
                xval-=11
            if yval >= 18*2:
                yval-=(18*2)
            if verbose >1:
                print("  > Adding cell at ",xval,yval," with ",en_fraction[i+2][j+2])
            HF[xval][yval]+=en_fraction[i+2][j+2]
            HF[xval][yval]=min(HF[xval][yval],256)

def getLinkStrings(HF):
    link_strings=["" for i in range(18)]
    for wIdx in range(18):
        for i in range(11):
            et=min(int(HF[i][2*wIdx]),256)
            link_strings[wIdx]+=f"{et:10b}".replace(" ","0")
        for i in range(11):
            et=min(int(HF[i][2*wIdx+1]),256)
            link_strings[wIdx]+=f"{et:10b}".replace(" ","0")
    return link_strings        


def writePatternFile(link_strings,ofname="output.links"):
    with open(ofname,"w") as f:
        print("writing out ",ofname)
        header=[f"{i:^256}" for i in range(18) ]
        f.write(" ".join(header))
        for link_string in link_strings:
            f.write("\n"+" ".join(link_string))

        f.close()



def writeDescriptionFile(HF,ofname="output.csv"):
    with open(ofname,"w") as f:
        print("writing out ",ofname)
        for wIdx in range(18):
            for i in range(11):

                et=min(int(HF[i][2*wIdx]),256)
                if et < 1:
                    continue
                ostr=f"{et},{i},{2*wIdx}"
                ostr+=f",{et:10b}".replace(" ","0")
                f.write(ostr+"\n")
                
                et=min(int(HF[i][2*wIdx+1]),256)
                if et < 1:
                    continue
                ostr=f"{et},{i},{2*wIdx+1}"
                ostr+=f",{et:10b}".replace(" ","0")
                f.write(ostr+"\n")


HF=np.zeros(shape=(11,18*2))

#########     CLUSTER ENERGY DEFINITION ############
fillCluster(HF,Energy=380.0,width=0.8,verbose=1)
fillCluster(HF,Energy=130.0,width=0.8,verbose=1)
fillCluster(HF,Energy=20.0,width=3,verbose=1)
fillCluster(HF,Energy=3.0,width=10,verbose=1)
for et in [10.0,12,30,23,50,30,60]:
    fillCluster(HF,Energy=et,width=2,verbose=1)
    
####################################################



f,ax=plt.subplots(figsize=(14,6))
HF_toPlot=np.array(HF)
HF_toPlot[HF_toPlot < 0.5 ] = np.inf
cbar=ax.matshow(HF_toPlot,cmap='Reds',vmax=251,vmin=0.0)
_=ax.set_xticks(np.arange(0,36,2)-0.5,minor=False,labels=[])
_=ax.set_xticks(np.arange(0,36,2)+0.5,minor=True,labels=np.arange(0,18))

_=ax.set_yticks(np.arange(0,11,1),minor=True ,labels=np.arange(0,11,1)+30)
_=ax.set_yticks(np.arange(-0.5,11,1),minor=False,labels=[])
ax.grid(which='major',alpha=0.8)
ax.grid(which='minor',axis='x',alpha=0.2)
plt.colorbar(cbar,fraction=0.015, pad=0.01)
ax.set_xlabel('link',fontsize=12,loc='right')
ax.set_ylabel('$I\eta$',fontsize=12,loc='top')
ax.xaxis.set_ticks_position('bottom')
ofname=f"../data/png/{args.tag}_seed{args.seed}.png"
print("exporting image ",ofname)
f.savefig(ofname,bbox_inches='tight')


link_strings=getLinkStrings(HF)
writePatternFile([link_strings],f"../data/inputPatterFiles/{args.tag}_seed{args.seed}.links")

writeDescriptionFile(HF,f"../data/csv/{args.tag}_seed{args.seed}.csv")
