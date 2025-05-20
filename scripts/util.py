import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches

TOWER_EDGE_COLOR=(0.05, 0.2, 0.5, 0.3)
TOWER_EDGE_COLOR=(0.05, 0.2, 0.5, 0.1)
PADDING=(213./256., 22/256., 238/256., 0.05)

cmap=mpl.colormaps["jet"]
sub_cmap=mpl.colormaps["Greys"]


eta_tow_map=np.arange(12+2)+10
phi_tow_map=np.linspace(0.0,360.0,72,endpoint=False) 

eta_st_map=np.arange(6)+4
phi_st_map=np.linspace(0.0,360.0,24,endpoint=False)
def assignIssuClass( hf_tower_grid  ):
    return

def unpackCalocluster(data):
    jdata=np.uint64(data)
    calo={}
    calo['eta'] = int((jdata >> np.uint64(12) ) & np.uint64(0x1f))
    calo['phi'] = int((jdata >> np.uint64(17) ) & np.uint64(0xff))
    calo['et']  = int(jdata & np.uint64(0xfff))
    return calo

def unpackJet(data):
    jdata=np.uint64(data)
    jet={}
    jet['et']  = int(jdata & np.uint64(0xfff))
    jet['eta'] = int((jdata >> np.uint64(12) ) & np.uint64(0x7))
    jet['phi'] = int((jdata >> np.uint64(15) ) & np.uint64(0x1f))
    jet['seedEt'] = int((jdata >> np.uint64(27) ) & np.uint64(0x3fff))
    return jet

def unpackPuppiCluster(data):
    items = data.split(',')
    if len(items) < 4:
        return None
    cluster = {
        'et': float(items[0]),
        'eta': int(items[1]),
        'phi': int(items[2]),
    }
    return cluster if cluster['et'] > 0 else None

def unpackTrueCluster(data):
    items = data.split(',')
    if len(items) < 4:
        return None
    cluster = {
        'type': items[0].strip(),
        'et': float(items[1]),
        'eta': int(items[2]),
        'phi': int(items[3])
    }
    return cluster if cluster['et'] > 0 else None

def isClusterMatch(c1,c2):
    if int(c1['eta'])!=int(c2['eta']) : return False
    if int(c1['phi'])!=int(c2['phi']) : return False
    if int(c1['et'])!=int(c2['et']) : return False
    return True

def getTrueCaloJets(hf_supertower_grid):
    print("HERE")
    mask={}
    for e in hf_supertower_grid:
        mask[e]={}
        for p in hf_supertower_grid[e]:
            mask[e][p]=False
    caloJet_truths=[]
    caloJet_truthsE=[]
    etaL=[]
    phiL=[]
    for e in hf_supertower_grid:
        etaL.append(e)
        if len(phiL) <1:
            for p in hf_supertower_grid[e]:        
                phiL.append(p)
    for idx in range(9):        
        emax=-1
        pC=-1
        eC=-1

        for p in phiL:
            for e in etaL:      
                if e < 0 :
                    continue
                if mask[e][p]:
                    continue
                if emax <= hf_supertower_grid[e][p]:
                    emax=hf_supertower_grid[e][p]
                    pC=p
                    eC=e
        if emax<0:
            emax=0
            pC=0
            eC=0
            break
        esum=0
        for i in [-1,0,1]:
            ii=pC+i
            if ii <0:
                ii=23
            if ii >23:
                ii=0
            for j in [-1,0,1]:
                jj=eC+j
                esum+=hf_supertower_grid[jj][ii]*(not mask[jj][ii])
                mask[jj][ii]=True
        if esum<=0:
            continue
        caloC={'et':esum,'eta':eC,'phi':pC,'seedEt':emax}
        print(f"Adding true calo jet {idx} : ",caloC)
        caloJet_truths.append(caloC)
        caloJet_truthsE.append(emax)
    srtidx=np.argsort(np.array(caloJet_truthsE)*-1)
    caloJet_truths=[caloJet_truths[i] for i in srtidx]
    
    return caloJet_truths

def getTrueCaloClusters(hf_tower_grid,nClusMax=-1):
    mask={}
    for p in hf_tower_grid:
        mask[p]={}
        for e in hf_tower_grid[p]:
            mask[p][e]=False
    caloCluster_truths=[]
    caloCluster_truthsE=[]
    caloCluster_truthsEseed=[]
    
    #print(" Search order ! (p,e) : ")
    #for p in hf_tower_grid:
    #    for e in hf_tower_grid[p]:        
    #        print(f"({p},{e}) | ",end="")
    #    print()
    print(" Cluster Making !  " )

    pList=[]
    etaList=[]
    for p in hf_tower_grid:
        pList.append(p)
    for e in hf_tower_grid[pList[0]]:       
        etaList.append(e)

    for idx in range(64):        
        emax=0
        pC=-1
        eC=-1
        for e in etaList:        
            for p in pList:
                if e < 0 :
                    continue
                if mask[p][e]:
                    continue
                if emax <= hf_tower_grid[p][e]:
                    #print(f"    --> ( {emax=} , {pC=} , {eC=} ) to ",end="")
                    emax=hf_tower_grid[p][e]
                    pC=p
                    eC=e
                    #print("  ( {emax=} , {pC=} , {eC=} )   ")
        if not(emax>4):
            emax=0
            pC=0
            eC=0
            break
        esum=0
        # print(f"CID {idx:>2} [p{pC:>2},eta{eC:>2},eMax {emax:>3}] --> ",end="")
        eseed=emax
        for i in [-1,0,1]:
            ii=pC+i
            if ii <0:
                ii=71
            if ii >71:
                ii=0
            for j in [-1,0,1]:
                jj=eC+j
                esum+=hf_tower_grid[ii][jj]*(not mask[ii][jj])
                # print(f"{hf_tower_grid[ii][jj]:>3}/{int(mask[ii][jj])} | ",end="")
                mask[ii][jj] = True
        # print()
        # print(f"{esum=} {esum/2=}")
        caloC={'eta':eC,'phi':pC,'et':esum,'sector' : int(pC/4),'cid':idx,'seed':emax}
    #     print(f"Adding cluster {idx} : ",caloC)
        caloCluster_truths.append(caloC)
        caloCluster_truthsE.append(esum)
        caloCluster_truthsEseed.append(eseed)
    # srtidx=np.argsort(np.array(caloCluster_truthsE)*-1)
    srtidx=np.argsort(np.array(caloCluster_truthsEseed)*-1)
    caloCluster_truths=[caloCluster_truths[i] for i in srtidx]
    if nClusMax>0:
        if nClusMax < len(caloCluster_truths):
            caloCluster_truths=caloCluster_truths[:nClusMax]
    return caloCluster_truths


def add_flushed_st_grid(hf_supertower_grid,ax):
    st_artists=[]
    emaxST=0.0
    for e in hf_supertower_grid:
        for p in hf_supertower_grid[e]:
            emaxST=max(emaxST,hf_supertower_grid[e][p])     
    v0=np.log2(emaxST)
    for e in hf_supertower_grid:
        for p in hf_supertower_grid[e]:
            r =eta_st_map[e+1]
            ph=phi_st_map[p]
            v=hf_supertower_grid[e][p]
            if v > 0.0:
                col=sub_cmap(0.2+np.log2(v)*(0.7-0.15)/v0)
            else:
                col='none'
    #         print(col)
            w=mpatches.Wedge((0, 0), r+0.07, ph+0.5, ph+15.0-0.7, ec="k",width=0.8)
            w.set(color=col)
            w.set_edgecolor(TOWER_EDGE_COLOR)
            if (e==-1) or (e==4):
                w.set_color(PADDING)
                w.set_edgecolor(TOWER_EDGE_COLOR)
            st_artists.append(w)   
        
    for st in st_artists:
        ax.add_artist( st )


def add_st_grid(hf_supertower_grid,ax):
    st_artists=[]
    emaxST=0.0
    for e in hf_supertower_grid:
        for p in hf_supertower_grid[e]:
            emaxST=max(emaxST,hf_supertower_grid[e][p])
    v0=np.log2(emaxST)
    for e in hf_supertower_grid:
        for p in hf_supertower_grid[e]:
            r =eta_st_map[e+1]
            ph=phi_st_map[p]
            v=hf_supertower_grid[e][p]
            w=mpatches.Wedge((0, 0), r+0.07, ph+0.5, ph+15.0-0.7, ec="k",width=0.8)
            ecol='none'
            if v > 0.0:
                col=cmap(np.log2(v)*0.8/v0)
                col=cmap(0.15+np.log2(v)*(0.9-0.15)/v0)
            else:
                col='gainsboro'
            if (e==-1) or (e==4):
                col =PADDING
                ecol=TOWER_EDGE_COLOR
    #         print(col)
            w.set(color=col)
            w.set_edgecolor(ecol)
            odc={'artist' : w ,'phi' :ph}
            if v >0:
                odc['label'] = f'{v}'
            st_artists.append( odc )
    #     break
    
    for i,  tow in enumerate( st_artists):
        ax.add_artist(tow['artist'])
        if 'label' in tow:
            if ph<270:
                ph+=270+20
            ax.annotate( tow['label'],(0.2,0.3), 
                        xycoords=tow['artist'],
                        rotation=ph)
def add_hf_towers(hf_tower_grid,ax,addLabel=False):
    artists=[]
    emax=0
    for p in hf_tower_grid:
        for e in hf_tower_grid[p]:
            emax=max(emax,hf_tower_grid[p][e])
    v0=np.log2(emax)
    for p in hf_tower_grid:
        for e in hf_tower_grid[p]:
            r =eta_tow_map[e+1]
            ph=phi_tow_map[p]
            v=hf_tower_grid[p][e]
            if v > 0.0:
                col=cmap(0.15+np.log2(v)*(0.9-0.15)/v0)
    #             col=cmap(np.log2(v)*(0.8)/v0)
            elif (e==-1) or (e==12):
                col=PADDING
            else:
                col='gainsboro'
    #         print(col)
            w=mpatches.Wedge((0, 0), r+0.07, ph+0.5, ph+5.0-0.7,ec="k",width=0.8)
            w.set(color=col)
    #         else:
    #             w.set(color='grey')
            ddict={'artist' : w ,'phi' :ph}
            if addLabel:
                if v > 0:
                    if True : #p%2:
                        ddict['label']=f"{v}"
            if (e==12) and (ph%2):
                ddict['label']=f"{int(ph/5)}"
            artists.append(ddict)
    #     break
    for i,  artistD in enumerate( artists):
        ax.add_artist(artistD['artist'])
        if 'label' in artistD:
            ax.annotate( artistD['label'],(0.2,0.2), 
                        xycoords=artistD['artist'],fontsize=8,
                        rotation=ph)

def add_flushed_hf_towers(hf_tower_grid,ax):
    artists=[]
    emax=0
    for p in hf_tower_grid:
        for e in hf_tower_grid[p]:
            emax=max(emax,hf_tower_grid[p][e])
    v0=np.log2(emax)
    for p in hf_tower_grid:
        for e in hf_tower_grid[p]:
            r =eta_tow_map[e+1]
            ph=phi_tow_map[p]
            v=hf_tower_grid[p][e]
            if v > 0.0:
                col=sub_cmap(0.15+np.log2(v)*(0.6-0.15)/v0)
    #             col=cmap(np.log2(v)*(0.8)/v0)
            elif (e==-1) or (e==12):
                col=PADDING
            else:
                col='none'
    #         print(col)
            w=mpatches.Wedge((0, 0), r+0.07, ph+0.5, ph+5.0-0.7,ec="k",width=0.8)
            w.set(color=col)
            w.set_edgecolor(TOWER_EDGE_COLOR)
            artists.append(w)
    #     break
    for i,  artist in enumerate( artists):
        ax.add_artist(artist)        

def add_hf_tower_sector(hf_tower_grid,ax,sector=1):
    addLabel=True
    artists=[]
    emax=0
    pmin=1e9
    plist_sector=[]
    fiducial_sector=[]
    for p in hf_tower_grid:
        skip=False
        if p <( sector*12      -2) : skip=True
        if p >( sector*12 +12 + 1) : skip=True
        if (sector==0) and (p > 69 ) :
            skip=False
        if (sector==5) and (p < 2 ):
            skip=False
        if skip : continue    
        plist_sector.append(p)
        fiducial_sector.append( (p >= sector*12) and ( p <sector*12 + 12 ) )
        pmin=min(pmin,phi_tow_map[p])
    emax=0
    for p in plist_sector:  
        for e in hf_tower_grid[p]:
            emax=max(emax,hf_tower_grid[p][e])           
    v0=np.log2(emax)
    phMap={}
    poff=1e9
    for p in plist_sector:  
        ph=phi_tow_map[p]
        ph-=pmin
        if ph > 180:
            ph-=360
        while (ph > 80 ):
            ph -=80
        ph+=5    
        phMap[p]=ph
        poff=min(poff,ph)
    phMap={p:phMap[p] - poff for p in phMap}
#     print(phMap)
#     return
    for p,isFiducial in zip(plist_sector,fiducial_sector):  
        ph=phMap[p]+5
        for e in hf_tower_grid[p]:
            r =eta_tow_map[e+1]
            v=hf_tower_grid[p][e]
            
            col ='none'
            ecol='none'
            
            if v > 0.0:
                col=cmap(0.15+np.log2(v)*(0.9-0.15)/v0)
            elif (e==-1) or (e==12):
                col=PADDING
            else:
                col='gainsboro'
            if not isFiducial:
                ecol=TOWER_EDGE_COLOR
                col='none'
                if v > 0.0:
                    col=sub_cmap(0.15+np.log2(v)*(0.7-0.15)/v0)
                elif (e==-1) or (e==12):
                    col=PADDING
            w=mpatches.Wedge((0, 0), r+0.07, ph+0.5, ph+5.0-0.7,ec="k",width=0.8)
            w.set(color=col)
            if ecol !='none':
                w.set_edgecolor(color=ecol)
    #         else:
    #             w.set(color='grey')
            ddict={'artist' : w ,'phi' :ph}
            if addLabel:
                if v > 0:
                    ddict['label']=f"{v}"
            if (e==12):
                ddict['label']=f"{int(p)}"
            artists.append(ddict)
    #     break
    for i,  artistD in enumerate( artists):
        ax.add_artist(artistD['artist'])
        if 'label' in artistD:
            ax.annotate( artistD['label'],(0.2,0.2), 
                        xycoords=artistD['artist'],fontsize=8,
                        rotation=ph)
    ax.set_ylim([0,25])
    ax.set_xlim([0,25])
    
def add_flushed_hf_tower_sector(hf_tower_grid,ax,sector=1,addLabel=False):
    artists=[]
    emax=0
    pmin=1e9
    plist_sector=[]
    fiducial_sector=[]
    for p in hf_tower_grid:
        skip=False
        if p <( sector*12      -2) : skip=True
        if p >( sector*12 +12 + 1) : skip=True
        if (sector==0) and (p > 69 ) :
            skip=False
        if (sector==5) and (p < 2 ):
            skip=False
        if skip : continue    
        plist_sector.append(p)
        fiducial_sector.append( (p >= sector*12) and ( p <sector*12 + 12 ) )
        pmin=min(pmin,phi_tow_map[p])
    emax=0
    for p in plist_sector:  
        for e in hf_tower_grid[p]:
            emax=max(emax,hf_tower_grid[p][e])           
    v0=np.log2(emax)
    phMap={}
    poff=1e9
    for p in plist_sector:  
        ph=phi_tow_map[p]
        ph-=pmin
        if ph > 180:
            ph-=360
        while (ph > 80 ):
            ph -=80
        ph+=5    
        phMap[p]=ph
        poff=min(poff,ph)
    phMap={p:phMap[p]+5 - poff for p in phMap}
#     print(phMap)
#     return
    for p,isFiducial in zip(plist_sector,fiducial_sector):  
        ph=phMap[p]
        for e in hf_tower_grid[p]:
            r =eta_tow_map[e+1]
            v=hf_tower_grid[p][e]
            
            col ='none'
            ecol=TOWER_EDGE_COLOR
            
            if v > 0.0:
                col=sub_cmap(0.15+np.log2(v)*(0.7-0.15)/v0)
            elif (e==-1) or (e==12):
                col=PADDING
            w=mpatches.Wedge((0, 0), r+0.07, ph+0.5, ph+5.0-0.7,ec="k",width=0.8)
            w.set(color=col)
            if ecol !='none':
                w.set_edgecolor(color=ecol)
    #         else:
    #             w.set(color='grey')
            ddict={'artist' : w ,'phi' :ph}
            if addLabel:
                if v > 0:
                    ddict['label']=f"{v}"
            if (e==12):
                ddict['label']=f"{int(p)}"
            artists.append(ddict)
    #     break
    for i,  artistD in enumerate( artists):
        ax.add_artist(artistD['artist'])
        if 'label' in artistD:
            ax.annotate( artistD['label'],(0.2,0.2), 
                        xycoords=artistD['artist'],fontsize=8,
                        rotation=ph)
    ax.set_ylim([0,25])
    ax.set_xlim([0,25])    
    return phMap

def add_sector_overlay(ax,n_sectors=6,rmin=5,rmax=10,addLabel=True,color=TOWER_EDGE_COLOR,phi_offset=0):
    sector_artists=[]
    phi_edges=np.linspace(0.0,360,num=n_sectors,endpoint=False)+phi_offset
    dPh=phi_edges[1]-phi_edges[0]
    width=(rmax-rmin)
    r    =rmax

    for idx,ph in enumerate(phi_edges):
        w=mpatches.Wedge((0, 0), r+0.07, ph+1, ph+dPh-1, ec="k",width=width)
        w.set(color='none')
        w.set_edgecolor(color)
        odc={'artist' : w ,'phi' :ph}
        if addLabel:
            odc['label'] = f'S{idx}'
        sector_artists.append( odc )

    for i,  sec in enumerate( sector_artists):
        ax.add_artist(sec['artist'])
        if 'label' in sec:
            ph=sec['phi']+dPh/8
            x=r*1.05*np.cos(np.pi*ph/180.0)
            y=r*1.05*np.sin(np.pi*ph/180.0)
            ax.annotate( sec['label'],(x,y), 
                        horizontalalignment='center',
                        verticalalignment='center',
                        rotation=0,color=color)
#     ax.set_xlim([-r*1.08,r*1.08])        
#     ax.set_ylim([-r*1.08,r*1.08])        
