import numpy as np
import argparse
import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('-i',"--inputFile" , help="Input file" , default=None)
parser.add_argument('-t',"--tag" , help="Tag" , default='opt')
args = parser.parse_args()

def unpackJet(data):
    jdata=np.uint64(data)
    jet={}
    jet['eta'] = (jdata >> np.uint64(12) ) & np.uint64(0x3)
    jet['phi'] = (jdata >> np.uint64(15) ) & np.uint64(0x1f)
    jet['et'] = jdata & np.uint64(0xfff)
    jet['seedEt'] = (jdata >> np.uint64(27) ) & np.uint64(0xfff)
    return jet

def unpackCalocluster(data):
    jdata=np.uint64(data)
    calo={}
    calo['eta'] = (jdata >> np.uint64(12) ) & np.uint64(0x5)
    calo['phi'] = (jdata >> np.uint64(17) ) & np.uint64(0xff)
    calo['et']  = jdata & np.uint64(0xfff)
    return calo

prefix='./'
tag=args.tag
with open(args.inputFile) as f:
    txt=f.readlines()
cmap=mpl.colormaps["jet"]

hf_tower_grid={}
emax=0
for l in txt:
    l=l[:-1]
    if not l.startswith('@@HFTowers'):
        continue
#     print(l)
    items=l.split("|")
    link,offset=items[1].split(",")
    iphi=int(link)*4+int(offset)
    hf_tower_grid[iphi]={}
    hf_tower_grid[iphi][-1]=0
    for i in range(12):
        hf_tower_grid[iphi][i]=int(items[2+i])
        emax=max(emax,int(items[2+i]))
    hf_tower_grid[iphi][12]=0

hf_supertower_grid={}
emaxST=0
for l in txt:
    l=l[:-1]
    if not l.startswith('@@HFSuperTowers'):
        continue
#     print(l)
    items=l.split("|")
    ieta=int(items[1])
    hf_supertower_grid[ieta]={}
    for i in range(24):
        hf_supertower_grid[ieta][i]=int(items[2+i])
        emaxST=max(emaxST,int(items[2+i]))

eta=np.arange(12+2)+10
phi=np.linspace(0.0,360.0,72,endpoint=False)

artists=[]
v0=np.log2(emax)
for p in hf_tower_grid:
    for e in hf_tower_grid[p]:
        r =eta[e+1]
        ph=phi[p]
        v=hf_tower_grid[p][e]
        if v > 0.0:
            col=cmap(0.15+np.log2(v)*(0.9-0.15)/v0)
#             col=cmap(np.log2(v)*(0.8)/v0)
        else:
            col='gainsboro'
#         print(col)
        w=mpatches.Wedge((0, 0), r+0.07, ph+0.5, ph+5.0-0.7,ec="k",width=0.8)
        w.set(color=col)
#         else:
#             w.set(color='grey')
        artists.append(w)
#     break
f,ax = plt.subplots(figsize=(5, 5), layout="constrained")
for i,  artist in enumerate( artists):
    ax.add_artist(artist)

Rmax=max(eta)
Rmin=max(eta)

# ax.set(title=type(artist).__name__,
#        aspect=1, xlim=(-1.2, 1.2), ylim=(-1.2, 1.2))
ax.set_axis_off()
ax.text(-4.5,-1,"HF Towers",fontsize=12)
ax.set(aspect=1)
ax.set_xlim(-1.05*Rmax,1.05*Rmax)
ax.set_ylim(-1.05*Rmax,1.05*Rmax)
foutname=f'{prefix}/hfTowers_{tag}.png'
print("Exporting : ",foutname)
plt.savefig(foutname,bbox_inches='tight')


eta=np.arange(6)+5
phi=np.linspace(0.0,360.0,24,endpoint=False)

artists=[]
v0=np.log2(emaxST)
for e in hf_supertower_grid:
    for p in hf_supertower_grid[e]:
        r =eta[e]
        ph=phi[p]
        v=hf_supertower_grid[e][p]
        if v > 0.0:
            col=cmap(np.log2(v)*0.8/v0)
            col=cmap(0.15+np.log2(v)*(0.9-0.15)/v0)
        else:
            col='gainsboro'
#         print(col)
        w=mpatches.Wedge((0, 0), r+0.07, ph+0.5, ph+15.0-0.7, ec="k",width=0.8)
        w.set(color=col)
        artists.append(w)
#     break
f,ax = plt.subplots(figsize=(5, 5), layout="constrained")
for i,  artist in enumerate( artists):
    ax.add_artist(artist)
ax.set_axis_off()
ax.text(-2.25,-0.5,"Super Towers",fontsize=12)
Rmax=max(eta)
Rmin=max(eta)
ax.set(aspect=1)
ax.set_xlim(-1.05*Rmax,1.05*Rmax)
ax.set_ylim(-1.05*Rmax,1.05*Rmax)
foutname=f'{prefix}/hfSuperTowers_{tag}.png'
print("Exporting : ",foutname)
plt.savefig(foutname,bbox_inches='tight')


jets=[]
emax=0
for l in txt:
    l=l[:-1]
    if not l.startswith('@@JETS'):
        continue
    items=l.split("|")
    for i in range(9):
        jt=unpackJet(items[2+i])
        jets.append(jt)
jets.reverse()

eta=np.arange(6)+5
phi=np.linspace(0.0,360.0,24,endpoint=False)

artists={}
v0=np.log2(emaxST)
for e in hf_supertower_grid:
    artists[e]={}
    for p in hf_supertower_grid[e]:
        r =eta[e]
        ph=phi[p]
        
        w=mpatches.Wedge((0, 0), r+0.07, ph+0.5, ph+15.0-0.7, ec="k",width=0.8)
        w.set(color='gainsboro')
        artists[e][p]=w

v0=max([ jt['et'] for jt in jets] )
emaxJ=v0
if v0 > 0:
    emaxJ=np.log2(v0)
    
for jt in jets:
    e=jt['eta'] + 1
    p=jt['phi']
    v=jt['et']
    print(f"Got jet with et = {v} , eta = {e} , phi = {p} ")
    if v > 0.0:
        x=0.15+np.log2(v)*(0.9-0.15)/emaxJ
        col=cmap(x)
    else:
        continue
    artists[e+1][p+1].set(color=col)    
    artists[e  ][p+1].set(color=col)        
    artists[e-1][p+1].set(color=col)        
    artists[e+1][p  ].set(color=col)    
    artists[e  ][p  ].set(color=col)        
    artists[e-1][p  ].set(color=col)
    artists[e+1][p-1].set(color=col)    
    artists[e  ][p-1].set(color=col)        
    artists[e-1][p-1].set(color=col)         
        
#     break
f,ax = plt.subplots(figsize=(5, 5), layout="constrained")
for e in artists:
    for p in artists[e]:
        ax.add_artist(artists[e][p])


        
Rmax=max(eta)
Rmin=max(eta)

# ax.set(title=type(artist).__name__,
#        aspect=1, xlim=(-1.2, 1.2), ylim=(-1.2, 1.2))
# ax.set_axis_off()
ax.set_axis_off()
ax.text(-0.75,-0.25,"Jets",fontsize=12)
ax.set(aspect=1)
ax.set_xlim(-1.05*Rmax,1.05*Rmax)
ax.set_ylim(-1.05*Rmax,1.05*Rmax)
foutname=f'{prefix}/hfJets_{tag}.png'
print("Exporting : ",foutname)
plt.savefig(foutname,bbox_inches='tight')

caloClusters=[]
emax=0
for l in txt:
    l=l[:-1]
    if not l.startswith('@@CALOCLUSTER'):
        continue
    items=l.split("|")
    for i in range(12):
        cl=unpackCalocluster(items[2+i])
        caloClusters.append(cl)
        
caloClusters.reverse()

eta=np.arange(12+2)+10
phi=np.linspace(0.0,360.0,72,endpoint=False)

artists={}
v0=np.log2(emax)
for p in hf_tower_grid:
    artists[p]={}
    for e in hf_tower_grid[p]:
        r =eta[e+1]
        ph=phi[p]
        v=hf_tower_grid[p][e]
        if v > 0.0:
            col=cmap(0.15+np.log2(v)*(0.9-0.15)/v0)
#             col=cmap(np.log2(v)*(0.8)/v0)
        else:
            col='gainsboro'
#         print(col)
        w=mpatches.Wedge((0, 0), r+0.07, ph+0.5, ph+5.0-0.7,fc='gainsboro',width=0.8)
        artists[p][e]=w
#         else:
#             w.set(color='grey')
#         artists.append(w)
#     break


v0=max([ ct['et'] for ct in caloClusters] )
emaxJ=v0
if v0 > 0:
    emaxJ=np.log2(v0)
    
for cc in caloClusters:
    e=cc['eta']
    p=cc['phi']
    v=cc['et']
    if v >0:
        print(f"Got calo cluster with et = {v} , eta = {e} , phi = {p} ")
    if v > 0.0:
        x=0.15+np.log2(v)*(0.9-0.15)/emaxJ
        col=cmap(x)
    else:
        continue
    artists[p+1][e+1].set(color=col)       
    artists[p  ][e+1].set(color=col)    
    artists[p-1][e+1].set(color=col)    
    artists[p+1][e  ].set(color=col)    
    artists[p  ][e  ].set(color=col) 
    artists[p-1][e  ].set(color=col) 
    artists[p  ][e-1].set(color=col)      
    artists[p-1][e-1].set(color=col)        
    artists[p+1][e-1].set(color=col)         
    


f,ax = plt.subplots(figsize=(5, 5), layout="constrained")
for p in artists:
    for e in artists[p]:
        ax.add_artist(artists[p][e])

Rmax=max(eta)
Rmin=max(eta)

# ax.set(title=type(artist).__name__,
#        aspect=1, xlim=(-1.2, 1.2), ylim=(-1.2, 1.2))
# ax.set_axis_off()
ax.set(aspect=1)
ax.set_axis_off()
ax.text(-5.25,-0.5,"Calo Clusters",fontsize=12)
ax.set_xlim(-1.05*Rmax,1.05*Rmax)
ax.set_ylim(-1.05*Rmax,1.05*Rmax)


foutname=f'{prefix}/hfCaloClusters_{tag}.png'
print("Exporting : ",foutname)
plt.savefig(foutname,bbox_inches='tight')




