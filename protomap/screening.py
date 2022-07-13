#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@file: sreening.py
@comment: 
@created: 2022/04/14 18:12:01
@auther: Zhu, Yinyue
@version: 1.0
'''

import numpy as np
import scipy as sp
from .proteins import ProteinGroups, Protein


def validate(protein:Protein,invert,method):
    markers=protein.g.markers
    ms=sorted((markers["m"].to_list()+[protein.m]),reverse=True)
    index=ms.index(protein.m)
    
    if invert==True:
        ctrl=np.flip(protein.ctrl())#上面是0，下面是20
        expr=np.flip(protein.expr())
    else:
        ctrl=protein.ctrl()
        expr=protein.expr()
    cm=ctrl.argmax()
    #ideal=np.zeros(len(ctrl))
    #ideal[cm]=1
    #ctrl@ideal/
    if cm>=len(ctrl)-2:
        return False
# mass vs migration
    if index==0:
        c1=bool(cm+1<=markers["f"][index]+1)
    elif index==len(ms)-1:
        c1=bool(cm+1>=markers["f"][index-1]-1)
    else:
        c1=bool(cm+1>=markers["f"][index-1]-1 and cm+1<=markers["f"][index]+1)
    if c1:
# noise peaks 
        c2=bool(len(protein.ints()[lambda df: df>0.01*protein.maxint])<20)
        if c2:
# peak shift 
            diff=ctrl.to_numpy()-expr.to_numpy()
            c3=bool(diff[cm-1:cm+2].sum()>0.1*ctrl[cm-1:cm+2].sum()\
                        and (diff[cm+2:].sum()<-0.1*expr[cm+2:].sum()))
            #try:
            #    c3=bool(diff[cm-1:cm+2].sum()>0.1*ctrl[cm-1:cm+2].sum()\
            #            and (diff[cm+2:].sum()<-0.1*expr[cm+2:].sum()))# 改成存在一个条带，上调水平超过什么？
            #except:
            #    c3=bool(diff[cm-1:cm+2].sum()>0.1*ctrl[cm-1:cm+2].sum())
            if c3:
# fraction shift
                if method=='mean':
                    v=np.arange(protein.g.nfracs)+1
                    cp=ctrl[cm:]@v[cm:]/ctrl[cm:].sum()
                    ep=expr[cm:]@v[cm:]/expr[cm:].sum()
                elif method=='max':
                    cp=ctrl.argmax()+1
                    ep=expr.argmax()+1
                c4=bool(ep>cp+1)
                if c4:
                    return True
    #print(locals())
    return False

def screen(proteins:ProteinGroups,method='mean',invert=False):
    cleavable=[]
    for protein in proteins:
        if validate(protein,method=method,invert=invert):
            cleavable.append(protein)
    return cleavable        
