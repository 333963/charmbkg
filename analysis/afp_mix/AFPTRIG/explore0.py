import pandas as pd, numpy as np, sys
D='../runs/lxplus/'
for run in (428770,427929):
    ev=pd.read_parquet(D+f'data_{run}_events.parquet')
    print('=====',run,len(ev))
    grl=ev.passGRL_MinBiasAFP==1
    tc=[c for c in ev.columns if c.startswith('trigPassed_')]
    print(ev[tc+[c for c in ev.columns if c.startswith('trigTBP')]].apply(lambda s:s.value_counts().to_dict()).to_string())
    afpc=[c for c in tc if 'L1AFP' in c]; non=[c for c in tc if 'L1AFP' not in c]
    afp=(ev[afpc]==1).any(axis=1); nonafp=(ev[non]==1).any(axis=1)
    ev['cls']=np.where(afp&~nonafp,'afp_only',np.where(afp&nonafp,'both',np.where(nonafp,'nonafp','none')))
    tag=(ev.afp_tag_A==1)|(ev.afp_tag_C==1)
    anyp=(ev.afp_nProtons_A>0)|(ev.afp_nProtons_C>0)
    print(pd.DataFrame({'n':ev[grl].groupby('cls').size(),'tag':tag[grl].groupby(ev.cls[grl]).mean(),'anyp':anyp[grl].groupby(ev.cls[grl]).mean(),
      'npu0':(ev.nPileupVertices[grl]==0).groupby(ev.cls[grl]).mean(),'hasD':(ev.n_cand_sel[grl]>0).groupby(ev.cls[grl]).mean(),'ncand>0':(ev.n_cand[grl]>0).groupby(ev.cls[grl]).mean(),
      'nsel':ev.nSelBL[grl].groupby(ev.cls[grl]).mean(),'TBP_AFP':(ev.trigTBP_L1_AFP_A_OR_C[grl]==1).groupby(ev.cls[grl]).mean()}))
    # combos of afp chains
    print(ev[grl].groupby([ev[c]==1 for c in afpc]+[nonafp[grl]]).size())
    print(ev[['nPileupVertices','nSelBL','n_cand','n_cand_sel','afp_nProtons_A','afp_nProtons_C','xi_cal_A','afp_xi1_A']].describe().T)
    a=ev[grl&(ev.cls=='afp_only')]
    print('afp_only nprot A x C'); print(pd.crosstab(a.afp_nProtons_A.clip(upper=3),a.afp_nProtons_C.clip(upper=3)))
    print('afp_only nSelBL quantiles',a.nSelBL.quantile([0,.05,.25,.5,.75,.95,1]).tolist())
    print('LB afp_only',a.lumiBlock.min(),a.lumiBlock.max(),a.lumiBlock.nunique())
    print('xi1_A quantiles (np>=1)',a.afp_xi1_A[a.afp_nProtons_A>=1].quantile([0,.01,.05,.25,.5,.75,.95,.99,1]).round(4).tolist())
    print('xi1_C quantiles (np>=1)',a.afp_xi1_C[a.afp_nProtons_C>=1].quantile([0,.01,.05,.25,.5,.75,.95,.99,1]).round(4).tolist())
    print('xical A q',a.xi_cal_A.quantile([0,.01,.05,.25,.5,.75,.95,.99,1]).round(4).tolist())
