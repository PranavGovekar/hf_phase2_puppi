
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-p","--prefix", help="Output prefix", default="hf_raw")
parser.add_argument("-m","--modes", help="Comma-separated list of modes", default="l")
parser.add_argument("-s","--seeds", nargs=2, type=int, 
                    help="Seed offset and number of vectors", default=[69, 1])
args = parser.parse_args()

modes = args.modes.split(',')
seed_offset, num_vectors = args.seeds
vectors_per_mode = num_vectors // len(modes)
remainder = num_vectors % len(modes)
seed_counter = seed_offset

def fillCluster(HF,Energy=30.0,width=0.8,verbose=10, kind="None"):

    vals=stats.norm.pdf(np.array([-2,-1,0,1,2])/width)/stats.norm.pdf(0)
    en_dist=np.outer(vals,vals)

    selct_=np.random.random(25).reshape(5,5)
    mask = selct_> en_dist
    en_fraction=np.array(en_dist)
    en_fraction[mask]=0.0
    en_fraction=en_fraction/np.sum(en_fraction)
    en_fraction=(en_fraction + en_fraction*np.random.normal(scale=0.05) )*Energy
    x0,y0=np.random.randint(11),np.random.randint(18*2)

    with open(f"./data/inputPatternFiles/{args.tag}.des","a") as f:
        f.write(f"@@TrueClusters | {kind}, {Energy:.2f}, {x0}, {y0*2}\n")

    if verbose >0:
        print("Adding cluster at ",x0,y0,f" of energy {Energy:.2f}",f" and width {width:.3f}")
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
            link_strings[wIdx]=f"{et:10b}".replace(" ","0") + link_strings[wIdx]
        for i in range(11):
            et=min(int(HF[i][2*wIdx+1]),256)
            link_strings[wIdx]=f"{et:10b}".replace(" ","0") + link_strings[wIdx]
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


# --------------------------
# Main Loop
# --------------------------


for mode_idx, mode in enumerate(modes):
    mode_count = vectors_per_mode + (1 if mode_idx < remainder else 0)

    if mode == 's':
        nPV_particle = 1
        nPU_particle = 2
        highEdepProb = 0.85
    elif mode == 'm':
        nPV_particle = 2
        nPU_particle = 4
        highEdepProb = 0.9
    elif mode == 'l':
        nPV_particle = 4
        nPU_particle = 8
        highEdepProb = 0.95
    elif mode == 'e':
        nPV_particle = 10
        nPU_particle = 12
        highEdepProb = 0.2

    for _ in range(mode_count):
        current_seed = seed_counter
        args.tag = f"{args.prefix}_{mode}_{current_seed}"
        np.random.seed(current_seed)

        HF = np.zeros(shape=(11, 18*2))
        nV = max(1, np.random.poisson(nPV_particle))
        for _ in range(nV):
            energy = np.random.uniform(10, 150)
            w = np.random.normal(0.50, 0.25)
            while w < 0.3:
                w = np.random.normal(0.5, 0.25)
            fillCluster(HF, Energy=energy, width=w, verbose=1, kind="PV")

        nV = max(2, np.random.poisson(nPU_particle))
        for _ in range(nV):
            energy = np.random.uniform(10, 80)
            w = np.random.normal(2.0, 0.7)  
            while w < 2.0:
                w = np.random.normal(2.0, 0.7)
            fillCluster(HF, Energy=energy, width=w, verbose=1, kind="PU")

        if np.random.uniform() > highEdepProb:
            nV = max(1, np.random.poisson(1))
            for _ in range(nV):
                energy = np.random.uniform(100, 400)
                w = max(2.0, np.random.normal(0.8, 0.3))
                fillCluster(HF, Energy=energy, width=w, verbose=1)

        f, ax = plt.subplots(figsize=(14, 6))
        HF_toPlot = np.array(HF)
        HF_toPlot[HF_toPlot < 0.5] = np.inf
        cbar = ax.matshow(HF_toPlot, cmap='Reds', vmax=251, vmin=0.0, origin='lower')
        ax.set_xticks(np.arange(0, 36, 2)-0.5, minor=False, labels=[])
        ax.set_xticks(np.arange(0, 36, 2)+0.5, minor=True, labels=np.arange(0, 18))
        ax.set_yticks(np.arange(0, 11, 1), minor=True, labels=np.arange(0, 11, 1))
        ax.set_yticks(np.arange(-0.5, 11, 1), minor=False, labels=[])
        ax.grid(which='major', alpha=0.8)
        ax.grid(which='minor', axis='x', alpha=0.2)
        plt.colorbar(cbar, fraction=0.015, pad=0.01)
        ax.set_xlabel('link', fontsize=12, loc='right')
        ax.set_ylabel('$I\eta$', fontsize=12, loc='top')
        ax.xaxis.set_ticks_position('bottom')

        ofname = f"./data/inputPatternFiles/{args.tag}.png"
        f.savefig(ofname, bbox_inches='tight')
        plt.close(f)
        
        link_strings = getLinkStrings(HF)
        writePatternFile([link_strings], f"./data/inputPatternFiles/{args.tag}.links")
        
        seed_counter += 1