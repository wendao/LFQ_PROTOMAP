#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@file: command.py
@comment: 
@created: 2022/04/14 18:00:19
@auther: Zhu, Yinyue
@version: 1.0
'''
import os
import pandas as pd
import argparse
from .mapping import ProtoMap
from .fileio import filehandle

parser = argparse.ArgumentParser(
    description="Select MaxQuant report files and database FASTA file.")
parser.add_argument('proteins',
                    metavar='pro',
                    help="MaxQuant output proteinGroups.txt")
parser.add_argument('peptides',
                    metavar='pep',
                    help="MaxQuant output peptides.txt")
parser.add_argument('fasta',
                    metavar='f',
                    help="FASTA database file which MaxQuant searches against")
parser.add_argument('-m',
                    '--markers',
                    default=None,
                    help="Tab-separated .txt file recording migrations of markers")
parser.add_argument('-o',
                    '--output',
                    default='./protomap/',
                    help="Protomap img output directory")
parser.add_argument('--invert',
                    action='store_true',
                    help="Invert fraction numbers")
parser.add_argument('--heatmap',
                    choices=["intensity", "spectral_count"],
                    default=None,
                    help="Heatmap type of peptides mapping")

def main():
    args =parser.parse_args()
    
    if os.path.isdir(args.output) == False:
        os.mkdir(args.output)
    proto=filehandle(vars(args))
    renderer = ProtoMap(invert_fraction=args.invert,
                    side_axis=True if isinstance(proto.markers,pd.DataFrame) else False,
                    heatmap=args.heatmap,
                    save_path=args.output)

    for i in proto:
        renderer.plot(i)
    return
    
if __name__ == '__main__':
    main()
