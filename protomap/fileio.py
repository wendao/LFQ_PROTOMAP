#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@file: fileio.py
@comment: 
@created: 2022/04/14 18:05:17
@auther: Zhu, Yinyue
@version: 1.0
'''

import os
import re
import pandas as pd
from Bio import SeqIO
from .proteins import ProteinGroups
from .screening import screen
from .uniprot_API import id_mapping

ONLINE=False

# processing input files
def filehandle(files:dict,**kargs):
    #if os.path.splitext(files['proteins'])
    proteins = pd.read_csv(files['proteins'],
                           sep="\t",
                           index_col=0,
                           low_memory=False)
    peptides = pd.read_csv(files['peptides'], sep="\t", header=0, index_col=0)
    if ONLINE:
        fasta = SeqIO.to_dict((SeqIO.parse(files['fasta'].temporary_file_path(), "fasta")))
    else:
        fasta = SeqIO.to_dict((SeqIO.parse(files['fasta'], "fasta")))
    try:
        markers = pd.read_csv(files['markers'], sep="\t", names=("m", "f"))
    except:
        markers = files['markers']
    return ProteinGroups(proteins, peptides, fasta, markers=markers,**kargs)

 
def file_to_set(path):
    file=open(path)
    return set(file.read().splitlines())

def dir_to_list(path):
    imgs=os.listdir(path)
    ids=[]
    for i in imgs:
        if os.path.splitext(i)[1]=='.png':
            name=i.split("_")[1]
            ids.append(name)
    return ids

def id_to_link(name:str,proteins:ProteinGroups):
    try:
        if name==None:
            return ''
        else:
            link=re.sub('\|', '_',str(proteins[name]))
            return f'<a href="./{link}.png">{name}</a>'
    except:
        return name
    
def screen(path:str,proteins:ProteinGroups,fmt:str,method='mean',invert=False):
    file=open(path,'w')
    
    if fmt=='html':
        s=0
        for protein in proteins:
            if validate(protein,method=method,invert=invert):
                s+=1
                name=re.sub('\|', '_',protein.name.split(';')[0])
                file.write(f'<p><a href="./{name}.png">{protein.name}</a></p>\n')
        file.write(f'<p>valid proteins:{s}/{len(proteins)}</p>')
    if fmt=='csv':
        for protein in proteins:
            if validate(protein,method=method,invert=invert):
                file.write(f'{protein.id()}\n')
    file.close()
    return 
