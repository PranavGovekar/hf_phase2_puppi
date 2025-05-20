import numpy as np
import argparse,json
import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

import util as utl


eta_st_map=utl.eta_st_map
eta_tow_map=utl.eta_tow_map
phi_st_map=utl.phi_st_map
phi_tow_map=utl.phi_tow_map
cmap=utl.cmap


parser = argparse.ArgumentParser()
parser.add_argument('-i',"--inputFile" , help="Input file" , default=None)
parser.add_argument('-t',"--tag" , help="Tag" , default='opt')
parser.add_argument('-c',"--doPFCluster" , help="make pfcluster validation" , action='store_true')
parser.add_argument('-j',"--doJets" , help="make jet validation" , action='store_true')
parser.add_argument(    "--processST", help="Process Super TOwer", action='store_true')
parser.add_argument( "--exportToJson", help="Export validation data to json", action='store_true')
parser.add_argument(     "--pHFT" , help="Print HF tower grid" , action='store_true')
parser.add_argument(     "--pST" , help="Print ST grid" , action='store_true')
parser.add_argument('-p',"--doPlots" , help="make all plots" , action='store_true')
parser.add_argument(     "--allObjects" , help="Print all objects" ,action='store_true')
parser.add_argument( "-e", "--expandSector" , help="Expand the sector" , default=None)
args = parser.parse_args()

exportData={"METADATA":{"tag" : args.tag}}

sector=None    
if args.expandSector is not None:
    sector=int(args.expandSector)

# prefix='./'
tag=args.tag
with open(args.inputFile) as f:
    txt=f.readlines()

## READING THE HF TOWERS
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
print("HF TOWERS -> IPHI : ",len(hf_tower_grid)," | IETA ",len(hf_tower_grid[0]))

exportData["hfTowers"]=hf_tower_grid

## READING THE CALO CLUSTERS    
caloClusters=[]
unpackedCaloClusters=[]
unpackedCaloClusterE=[]
emax=0
for l in txt:
    l=l[:-1]
    if not l.startswith('@@CALOCLUSTER'):
        continue
    items=l.split("|")
    sector_=int(items[1])
    sector_=int(items[1])
    caloClusters.append([])
    for i in range(12):
        cl=utl.unpackCalocluster(items[2+i])
        if cl['et'] < 1:
            continue
        cl['sector']=sector_
        cl['sector']=int( cl['phi'] /4 )
        caloClusters[-1].append(cl)
        unpackedCaloClusters.append(cl)
        unpackedCaloClusterE.append(int(cl['et']))

eta=np.arange(12+2)+10
phi=np.linspace(0.0,360.0,72,endpoint=False)
srtIdx=np.argsort(np.array(unpackedCaloClusterE)*-1)
unpackedCaloClusters=[unpackedCaloClusters[i] for i in srtIdx] 

## READING THE Number of calo clusters in overlap found per sectors
up_overlap_per_sector={}
for l in txt:
    l=l[:-1]
    if not l.startswith('@@CLUSTERSINOVFOUND'):
        continue
    items=l.split("|")
    count=int(items[1])
    sec  =int(items[2])
    up_overlap_per_sector[sec]=count


## READ PUPPI CLUSTERS
puppiClusters = []
for l in txt:
    l = l[:-1]
    if not l.startswith('@@PUPPI'):
        continue
    items = l.split("|")
    sector = int(items[1].strip())
    for i in range(2, 10):
        cluster_data = items[i].strip()
        cluster = utl.unpackPuppiCluster(cluster_data)
        if cluster is not None:
            cluster['sector'] = sector
            puppiClusters.append(cluster)

## READ TRUE CLUSTERS
trueClusters = []
for l in txt:
    l = l[:-1]
    if not l.startswith('@@TrueClusters'):
        continue
    items = l.split("|")
    cluster_data = items[1].strip()
    cluster = utl.unpackTrueCluster(cluster_data)
    if cluster is not None:
        trueClusters.append(cluster)

exportData["trueClusters"] = trueClusters
print(f"Found {len(trueClusters)} True clusters")

exportData["puppiClusters"] = puppiClusters
print(f"Found {len(puppiClusters)} PUPPI clusters")

## CALO CLUSTERS CLOSURE CHECKS
if args.pHFT:
    print()
    print("HF Tower Grid")
    pList=[]
    etaList=[]
    for p in hf_tower_grid:
        pList.append(p)
    for p in pList:
        for e in hf_tower_grid[p]:
            etaList.append(e)
        break
    
    _plist=[]
    offset=-1
    for i,p in enumerate(pList):
        if i> offset:
            offset+=24
            _plist.append([])
        _plist[-1].append(p)
    
    for i in range(len(_plist)):
        _txt="PHI -> | "
        _txt=f"{_txt:>13}"
        print(_txt,end="")
        for p in _plist[i]:
            _txt=f"{p:>3}|"
            print(_txt,end="")
        print()
        for e in etaList:
            _txt=f"eta = {e:>3} | "
            _txt=f"{_txt:>13}"
            print(_txt,end="")
            for p in _plist[i]:
                en=hf_tower_grid[p][e]
                _txt=f"{en:>3}|"
                print(_txt,end="")
            print()
        print()
    print()

    pass
caloClusterValidationSummary={}
if args.doPFCluster:
    caloCluster_truths=utl.getTrueCaloClusters(hf_tower_grid,nClusMax=48)
    exportData["trueCaloClusters"]=caloCluster_truths
    print()
    print("Calo Cluster validation ")
    print("  Number of true calo clusters = ",len(caloCluster_truths))
    print("  Number of unpacked calo clusters from FW = ",len(unpackedCaloClusters))
    true_sector_count={}
    for cc in caloCluster_truths:
        if cc['sector'] not in true_sector_count:
            true_sector_count[ cc['sector']  ] =0
        true_sector_count[ cc['sector']  ]+=1
    up_sector_count={}
    for cc in unpackedCaloClusters:
        if cc['sector'] not in up_sector_count:
            up_sector_count[ cc['sector']  ] =0
        up_sector_count[ cc['sector']  ]+=1
    
    secs=list(up_sector_count.keys()) ; secs.sort()
    for sec in secs:
        break
        print(f" -- sector {sec:>2} | true clusters : {true_sector_count[sec]} , up clusters : {up_sector_count[sec]:>2} , up clus from overlap : {up_overlap_per_sector[sec]} | [ DIFF : {up_sector_count[sec] - true_sector_count[sec]} ]")
    for sec in true_sector_count:
        if sec not in up_sector_count:
            print("sector {sec} in true cluster set not found in unpacked")
    
    if args.allObjects:
        for c1 in caloCluster_truths:
            print(f"    True :  et {c1['et']} , eta {c1['eta']} , phi {c1['phi']}  , sector  : {c1['sector']} , seed : {c1['seed']}   [cid:{c1['cid']}]")
        for c1 in unpackedCaloClusters:
            print(f"      UP :  et {c1['et']} , eta {c1['eta']} , phi {c1['phi']}  , sector  : {c1['sector']}  ")

    if args.expandSector is not None:
        smin= sector*12-2 
        smax= sector*12+12+1 
        print("Expanding sector ",args.expandSector, f" phi in [{smin},{smax}] ")
        has=False
        phiList=[]
        for p in hf_tower_grid:
            skip=False
            if p < smin : skip=True
            if p > smax : skip=True
            if skip: continue
            phiList.append(p)
        if sector==0: phiList=[70,71] + phiList
        if sector==5: phiList=phiList+[0,1]
        for i,p in enumerate(phiList):
            if i in [0,2,14]:
                print("---"*32)
            if i in [6,10]:
                print("- -"*32)
            if not has:
                print(f"phi x eta | ",end="")
                for e in hf_tower_grid[p]:
                    print(f"{e:>3} | ",end="")
                print("\n"+"---"*32)
                has=True
            print(f"phi = {p:>3} | ",end="")
            for e in hf_tower_grid[p]:
                x=hf_tower_grid[p][e]
                if x==0: x=' '
                print(f"{x:>3} | ",end="")
            print()
        print("---"*32)
        print("")
        print("  Unpacked calo clusters : ")
        for c1 in unpackedCaloClusters:
            if c1['sector']==sector:
                print(f"    UP :  et {c1['et']} , eta {c1['eta']} , phi {c1['phi']}  , sector  : {c1['sector']}")

        print("")
        print("  True calo clusters : ")
        for c1 in caloCluster_truths:
            if c1['sector']==sector:
                print(f"    True :  et {c1['et']} , eta {c1['eta']} , phi {c1['phi']}  , sector  : {c1['sector']} , seed : {c1['seed']}  [cid:{c1['cid']}]")

    k=0
    print("\n  Checking closure ! ")
    closure=True
    print("   True clusters not found in Unpacked set : ")
    listOfTrueClustersMissingInUP=[]
    for c1 in caloCluster_truths:
        isFound=False
        if sector is not None:
            if c1['sector']!=sector: continue
        for c2 in unpackedCaloClusters:
            if utl.isClusterMatch(c1,c2):
                isFound=True
                break
        if not isFound:
            closure=False
            listOfTrueClustersMissingInUP.append(c1)
            k+=1
            print(f"    {k} | et {c1['et']} , eta {c1['eta']} , phi {c1['phi']}  , sector  : {c1['sector']} seed : {c1['seed']}  [cid:{c1['cid']}]" )
    caloClusterValidationSummary["nTrueCaloCluterMissingInUP"]=k
    exportData["METADATA"]["nTrueMissingCaloClustersFromUP"] = len(listOfTrueClustersMissingInUP)
    exportData["misingTrueCaloCalusters"]=listOfTrueClustersMissingInUP
    print()        
    k=0
    print("   Unpacked clusters not found in Truth set")
    for c1 in unpackedCaloClusters:
        isFound=False
        if sector is not None:
            if c1['sector']!=sector: continue
        for c2 in caloCluster_truths:
            if utl.isClusterMatch(c1,c2):
                isFound=True
                break
        if not isFound:
            closure=False
            k+=1
            print(f"    {k} | et {c1['et']} , eta {c1['eta']} , phi {c1['phi']}  , sector  : {c1['sector']}" )
    
    caloClusterValidationSummary["nUPCaloCluterMissingInTrue"]=k
    if closure:
        print("  PFCluster Closure sucessfull")
    else:
        print("  PFCluster Closure Failed")

if args.processST:
    ## READING THE HF SUPER TOWERS
    hf_supertower_grid = {}
    emaxST = 0
    iSTPhi=-1
    for l in txt:
        l = l[:-1]
        if not l.startswith('@@HFSuperTowers'):
            continue
    #     print(l)
        iSTPhi+=1
        items = l.split("|")
        # ieta = int(items[1])
        # if (ieta > 5) or (ieta < 0):
            # print(
                # f"Warning !  eta value in input file for the Super Tower. eta={ieta}  , not inside [0,5]")
        p=iSTPhi
        for i in range(len(items[2:])):
            word = items[2+i]
            ieta = i-1
            if len(word.strip()) < 1:
                continue
            if ieta not in hf_supertower_grid :
                hf_supertower_grid[ieta] = {}
            if (ieta > 4) or (ieta < -1):
                print(           f"UNEXPETCTED eta value in input file for the Super Tower. eta={ieta}  , not inside [-1,4]")
                exit(1)
            hf_supertower_grid[ieta][p] = int(word)
            # print(ieta,p, i, word)
            emaxST = max(emaxST, int(items[2+i]))
    print("HF Super Towers -> IETA : ", len(hf_supertower_grid),
        " | IPHI ", len(hf_supertower_grid[0]))
    if args.pST:
        print()
        print("Super Tower Grid")
        _txt = "PHI -> | "
        _txt = f"{_txt:>13}"
        print(_txt, end="")
        for e in hf_supertower_grid:
            for p in hf_supertower_grid[e]:
                _txt = f"{p:>3}|"
                print(_txt, end="")
            break
        print()
        for e in hf_supertower_grid:
            _txt = f"eta = {e:>3} | "
            _txt = f"{_txt:>13}"
            print(_txt, end="")
            for p in hf_supertower_grid[e]:
                en = hf_supertower_grid[e][p]
                _txt = f"{en:>3}|"
                print(_txt, end="")
            print()

        print()
        print()

    jets = []
    emax = 0
    for l in txt:
        l = l[:-1]
        if not l.startswith('@@JETS'):
            continue
        items = l.split("|")
        for i in range(9):
            jt = utl.unpackJet(items[2+i])
            if jt['et'] >0:
                jets.append(jt)
                print(f"UP Jets {i} | {jt}")
    unpackedJets = jets


    caloJetValidationSummary={}
    if args.doJets:
        caloJet_truths=utl.getTrueCaloJets(hf_supertower_grid)
        print()
        print("Calo Jet validation ")
        print("  Number of true jets = ",len(caloJet_truths))
        print("  Number of unpacked jets from FW = ",len(unpackedJets))
        print("  Checking closure ! ")
        k=0
        closure=True
        #for jt in caloJet_truths:
        #    print(f"True Jets {i} | {jt}")
        for c1 in caloJet_truths:
            isFound=False
            for c2 in unpackedJets:
                if utl.isClusterMatch(c1,c2):
                    isFound=True
                    break
            if not isFound:
                closure=False
                k+=1
                print(f"    {k} True Jet not found in UP! ",c1)
        caloJetValidationSummary["nTureMissingInUP"]=k
        print()        
        k=0
        for c1 in unpackedJets:
            isFound=False
            for c2 in caloJet_truths:
                if utl.isClusterMatch(c1,c2):
                    isFound=True
                    break
            if not isFound:
                closure=False
                k+=1
                print(f"    {k} UP Jet did not find a match ! ",c1)
        caloJetValidationSummary["nUPMissingInTrue"]=k
        if closure:
            print("  JET Closure sucessfull")
        else:
            print("  JET Closure Failed")
    print()




if args.doPlots:
    ## SUPER TOWER SUMMARY PLOT
    if True:
        f,axes_to_plot = plt.subplots(1,3,figsize=(21, 6), layout="tight",dpi=120)
        Rmax=max(eta_st_map)
        Rmin=min(eta_st_map)-1
        for ax in axes_to_plot:
            ax.set_axis_off()
            ax.set(aspect=1)
            ax.set_xlim(-1.05*Rmax,1.05*Rmax)
            ax.set_ylim(-1.05*Rmax,1.05*Rmax)
        
        ax0,ax,axTrue=axes_to_plot
        utl.add_st_grid(hf_supertower_grid,ax0)
        utl.add_sector_overlay(ax0,2,Rmin,Rmax,color=(0,0,0,0.5))
        utl.add_sector_overlay(ax,2,Rmin,Rmax,color=(0,0,0,0.5))
        utl.add_flushed_st_grid(hf_supertower_grid,ax)        
        
        utl.add_sector_overlay(axTrue,2,Rmin,Rmax,color=(0,0,0,0.5))
        utl.add_flushed_st_grid(hf_supertower_grid,axTrue)        

        caloJet_truths = utl.getTrueCaloJets(hf_supertower_grid)

        for jetC,ax_ in zip([jets,caloJet_truths],[ax,axTrue]):
            v0=max([ jt['et'] for jt in jetC] )
            emaxJ=v0
            if v0 > 0:
                emaxJ=np.log2(v0)
                
            jet_coll=[]
            for jid,jt in enumerate(jetC):
                if jid>7: break
                e=jt['eta']
                p=jt['phi']
                v=jt['et']
                if v > 0.0:
                    x=0.15+np.log2(v)*(0.9-0.15)/emaxJ
                    col=cmap(x)
                else:
                    continue
                r =eta_st_map[e+1]
                ph=phi_st_map[p]
            #     print(e,p,v,r,ph)
                
                w=mpatches.Wedge((0, 0), r+1+0.07, ph-15*1 + 0.1, ph+15.0*2-0.7, 
                                edgecolor="k",width=3-0.2)
                w.set(color=col,alpha=0.5)
                w.set_edgecolor("k")
                jet_coll.append(
                        {'artist' : w ,'label' : f'[{jid+1}] {v}' ,'phi' :ph}
                )
            jet_coll.reverse()     
                    
            
            for ja in jet_coll:
                ax_.add_artist( ja['artist'] )
                ph=ja['phi']-90
                ph=0
                ax_.annotate( ja['label'],(0.4,0.4), 
                            xycoords=ja['artist'],
                            rotation=ph)
            
        
        ax0.text(-2.2,-0.25," Super Towers",fontsize=12)
        
        txt=f"Unpacked [{len(jets)}]"
        if "nUPMissingInTrue" in caloJetValidationSummary:
            txt+=f"\n nMissing : {caloJetValidationSummary['nUPMissingInTrue']}"
        ax.text(-2.0,-0.25,txt,fontsize=12)
        
        txt=f"True [{len(caloJet_truths)}] "
        if "nTureMissingInUP" in caloJetValidationSummary:
            txt+=f"\n nMissing : {caloJetValidationSummary['nTureMissingInUP']}"
        axTrue.text(-2.0,-0.25,txt,fontsize=12)
        axTrue.text(8.0,-9.0, f"{tag} " , bbox=dict(boxstyle="square",ec=None, fc='yellow'),fontsize=18)
        f.subplots_adjust(wspace=-0.5,hspace=-1.0)
        
        foutname=f'./outputs/stJets/stJets_{tag}.png'
        print("Exporting : ",foutname)
        plt.savefig(foutname,bbox_inches='tight')
    ## CALO CLUSTER SUMMARY PLOT
    
    if True:
        f,axes_to_plot = plt.subplots(1,3,figsize=(27*1.5, 7.5*1.5), layout="tight",dpi=150)
                
        Rmax=max(eta_tow_map)
        Rmin=min(eta_tow_map)
        
        for ax in axes_to_plot:
            ax.set_axis_off()
            ax.set(aspect=1)
            ax.set_xlim(-1.1*Rmax,1.1*Rmax)
            ax.set_ylim(-1.1*Rmax,1.1*Rmax)
        
        ax0,ax,axTrueClusters=axes_to_plot
        
        FONTSIZE_C = 12
        
        utl.add_hf_towers(hf_tower_grid,ax0,addLabel=True)       

        utl.add_flushed_hf_towers(hf_tower_grid,ax)
        utl.add_flushed_hf_towers(hf_tower_grid,axTrueClusters)
        
        utl.add_sector_overlay(ax0,18,Rmin-0.8,Rmax+0.8,color=(0,0,0,0.5))        
        utl.add_sector_overlay(ax ,18,Rmin-0.8,Rmax+0.8,color=(0,0,0,0.5))        
        utl.add_sector_overlay(axTrueClusters,18,Rmin-0.8,Rmax+0.8,color=(0,0,0,0.5))        
        
        v0=0
        for caloClusters_in_region in caloClusters:
            for cid,cc in enumerate(caloClusters_in_region):
                v0=max(v0,cc['et'])
        emaxJ=v0
        if v0 > 0:
            emaxJ=np.log2(v0)
        
        caloCluster_truths=[utl.getTrueCaloClusters(hf_tower_grid,nClusMax=48)]
        for axX, caloColl in zip([ax, axTrueClusters], [caloClusters, caloCluster_truths]):

            calo_collection=[]
            for rid, caloClusters_in_region in enumerate(caloColl):
                for cid,cc in enumerate(caloClusters_in_region):
                    e=cc['eta']
                    p=cc['phi']
                    v=cc['et']
                    if v > 0.0:
                        x=0.15+np.log2(v)*(0.9-0.15)/emaxJ
                        col=cmap(x)
                    else:
                        continue
                    r =eta_tow_map[e+1]
                    ph=phi_tow_map[p]
                    w=mpatches.Wedge((0, 0), r+1, ph-5 + 0.1, ph+10.0-0.7, 
                                    edgecolor="k",width=3-0.2)
                    w.set(color=col,alpha=0.5)
                    w.set_edgecolor("k")
                    ddict={'artist' : w ,'phi' :ph}
                    # if cid < 8:
                    ddict['label']=f'{v}'
                    calo_collection.append(
                            ddict
                    )
            #     break
            print(f"Adding {len(calo_collection)} clusters to image")
            for cc in calo_collection:
                axX.add_artist(cc['artist'])
                if 'label' in cc:
                    ph=0
                    axX.annotate( cc['label'],(0.4,0.4), 
                                xycoords=cc['artist'],
                                rotation=ph, fontsize=FONTSIZE_C)
        
        Rmax=max(eta)
        Rmin=max(eta)
        
        ax0.text(-5,-0.25,"HF Towers",fontsize=20)

        k=0
        for a in caloClusters:
            k+=len(a)
        txt=f"   Unpacked\nCalo Clusters [{k}]"
        if "nTrueCaloCluterMissingInUP" in caloClusterValidationSummary:
            txt+=f"\n nMissing {caloClusterValidationSummary['nTrueCaloCluterMissingInUP']}"
        ax.text(-5.0,-0.25,txt,fontsize=20)
        k=0
        for a in caloCluster_truths:
            k+=len(a)
        txt=f"   True\nCalo Clusters {k}"
        if "nUPCaloCluterMissingInTrue" in caloClusterValidationSummary:
            txt+=f"\n nMissing {caloClusterValidationSummary['nUPCaloCluterMissingInTrue']}"
        axTrueClusters.text(-5.0,-0.25,txt,fontsize=20)
        
        f.subplots_adjust(wspace=-0.55,hspace=-1.0)

        foutname=f'./outputs/caloClusterSummary/caloClusterSummary_{tag}.png'
        print("Exporting : ",foutname)
        plt.savefig(foutname,bbox_inches='tight')
        # exit()
    ## CALO CLUSTER DESCRIPTIVE
    if False:
        for sector in range(6):
            f,axs_to_plot = plt.subplots(1,3,figsize=(15, 5), layout="tight")
            axTW  =axs_to_plot[0]
            axCalo=axs_to_plot[1]
            axCaloTrue=axs_to_plot[1]
            utl.add_hf_tower_sector(hf_tower_grid,axTW,sector=sector)
            phMap=utl.add_flushed_hf_tower_sector(hf_tower_grid,axCalo,sector=sector)
            # phMap=utl.add_flushed_hf_tower_sector(hf_tower_grid,ax,sector=sector)
            calo_collection=[]
            ax=axCalo
            for cid,cc in enumerate(caloClusters[sector]):
                e=cc['eta']
                p=cc['phi']
                v=cc['et']
                if v > 0.0:
                    x=0.15+np.log2(v)*(0.9-0.15)/emaxJ
                    col=cmap(x)
                else:
                    continue
                r =eta_tow_map[e+1]
                ph=phMap[p]
                w=mpatches.Wedge((0, 0), r+1, ph-5 + 0.1, ph+10.0-0.7, 
                                 edgecolor="k",width=3-0.2)
                w.set(color=col,alpha=0.5)
                w.set_edgecolor("k")
                ddict={'artist' : w ,'phi' :ph}
                if cid < 8:
                    ddict['label']=f'{v}'
                calo_collection.append(
                        ddict
                )
            #     break
            for cc in calo_collection:
                ax.add_artist(cc['artist'])
                if 'label' in cc:
                    ph=0
                    ax.annotate( cc['label'],(0.4,0.4), 
                                xycoords=cc['artist'],
                                rotation=ph)

            print(Rmax, Rmin)
            utl.add_sector_overlay(axTW,18,Rmax-14,Rmax+0.8,color=(0,0,0,0.5),phi_offset=-5)
            utl.add_sector_overlay(axCalo,18,Rmax-14,Rmax+0.8,color=(0,0,0,0.5),phi_offset=-5)
            axTW.text(17,20,f"Sector {sector}\nCaloTowers",fontsize=15)
            axCalo.text(17,20,f"Sector {sector}\nCaloClusters",fontsize=15)
            axTW.set_xticks([])
            axCalo.set_xticks([])
            foutname=f'{prefix}/caloClusters_S{sector}_{tag}.png'
            print("Exporting : ",foutname)
            plt.savefig(foutname,bbox_inches='tight')

## PUPPI POLT
if True and (len(caloClusters) > 0 or len(trueClusters) > 0 or len(puppiClusters) > 0):
    f = plt.figure(figsize=(36, 11.25), dpi=150)
    gs = mpl.gridspec.GridSpec(1, 3, width_ratios=[1,1,1],
                             wspace=-0.22,
                             left=0.05, right=0.95,
                             bottom=0.05, top=0.95)
    
    Rmax = max(eta_tow_map)
    Rmin = min(eta_tow_map)
    
    axCalo = plt.subplot(gs[0])
    axTrue = plt.subplot(gs[1])
    axPuppi = plt.subplot(gs[2])
    
    for ax in [axCalo, axTrue, axPuppi]:
        ax.set_axis_off()
        ax.set(aspect=1)
        ax.set_xlim(-1.1*Rmax, 1.1*Rmax)
        ax.set_ylim(-1.1*Rmax, 1.1*Rmax)
        utl.add_flushed_hf_towers(hf_tower_grid, ax)
        utl.add_sector_overlay(ax, 18, Rmin-0.8, Rmax+0.8, color=(0,0,0,0.5))
 
    v0 = max([cc['et'] for region in caloClusters for cc in region], default=1)
    emaxJ = np.log2(v0) if v0 > 0 else 1
    
    for region in caloClusters:
        for cc in region:
            if cc['et'] > 0:
                e = cc['eta']
                p = cc['phi']
                r = eta_tow_map[e+1]
                ph = phi_tow_map[p]
                x = 0.15 + np.log2(cc['et'])*(0.9-0.15)/emaxJ
                w = mpatches.Wedge((0,0), r+1, ph-5+0.1, ph+10-0.7, 
                                 width=3-0.2, ec="k", fc=cmap(x), alpha=0.5)
                axCalo.add_artist(w)
                axCalo.annotate(f"{cc['et']:.1f}", (0.4,0.4), 
                               xycoords=w, fontsize=12)
    
    axCalo.text(-5.0, -0.25, f"Calo Clusters\n{sum(len(r) for r in caloClusters)}", fontsize=20)

    emax_true = max([c['et'] for c in trueClusters], default=1)
    emaxJ = np.log2(emax_true) if emax_true > 0 else 1
    
    pv_count = pu_count = 0
    for cc in trueClusters:
        if cc['et'] > 0:
            e = cc['eta']
            p = cc['phi']
            r = eta_tow_map[e+1]
            ph = phi_tow_map[p]
            
            if cc['type'] == 'PV':
                x = 0.15 + np.log2(cc['et'])*(0.9-0.15)/emaxJ
                color = cmap(x)
                pv_count += 1
            else:  # PU
                color = (0.5, 0.7, 0.8, 0.5)
                pu_count += 1
            
            w = mpatches.Wedge((0,0), r+1, ph-5+0.1, ph+10-0.7,
                             width=3-0.2, ec="k", fc=color)
            axTrue.add_artist(w)
            axTrue.annotate(f"{cc['et']:.1f}", (0.4,0.4), 
                           xycoords=w, fontsize=12)
    
    axTrue.text(-5.0, -0.25, f"True Clusters\nPV: {pv_count}, PU: {pu_count}", fontsize=20)
    
    emax_puppi = max([c['et'] for c in puppiClusters], default=1)
    emaxJ = np.log2(emax_puppi) if emax_puppi > 0 else 1
    
    for cc in puppiClusters:
        if cc['et'] > 0:
            e = cc['eta']
            p = cc['phi']
            r = eta_tow_map[e+1]
            ph = phi_tow_map[p]
            x = 0.15 + np.log2(cc['et'])*(0.9-0.15)/emaxJ
            w = mpatches.Wedge((0,0), r+1, ph-5+0.1, ph+10-0.7,
                             width=3-0.2, ec="k", fc=cmap(x), alpha=0.5)
            axPuppi.add_artist(w)
            axPuppi.annotate(f"{cc['et']:.1f}", (0.4,0.4), 
                           xycoords=w, fontsize=12)
    
    axPuppi.text(-5.0, -0.25, f"PUPPI Clusters\n{len(puppiClusters)}", fontsize=20)
    
    f.text(0.5, 0.02, f"{tag}", 
         bbox=dict(boxstyle="square", ec=None, fc='yellow'),
         fontsize=18, ha='center', va='bottom')
    
    print(f'Exporting : ./outputs/clusterSummary/clusterSummary_{tag}.png')
    plt.savefig(f'./outputs/clusterSummary/clusterSummary_{tag}.png', bbox_inches='tight')


if args.exportToJson:
    fname=f'validation_data_{args.tag}.json'
    with open(fname,'w') as f:
        print("Exporting ",fname)
        json.dump(exportData,f,indent=4)
print()
