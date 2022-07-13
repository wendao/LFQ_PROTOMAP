#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@file: proteins.py
@comment: 
@created: 2022/04/14 17:11:21
@auther: Zhu, Yinyue
@version: 1.0
'''

import re
import numpy as np
import pandas as pd
from typing import Dict, List, Union, NewType

#import statsmodels.formula.api as smf
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# filtering proteins by intensity
MININT = 1E7
#plt.rc('font', family='Arial', size=12)
ID = NewType('ID', str)


class ProteinGroups():
    '''
    Class for MaxQuant outputs from a combined search.

    pros: proteinGroups file as DataFrame
    peps: peptides file as DataFrame
    fasta: dictionary-like, keys: str(protein ID), values: SeqRecord object
    markers(optional): marker migration file as DataFrame, line: <molecular weight>\t <fraction>
    filt: filter contaminants and reverse sequences
    filtints: filter off proteins with low intensity
    invertgroup: if True, control is the first.
    '''

    def __init__(self,
                 pros: pd.DataFrame,
                 peps: pd.DataFrame,
                 fasta: Dict[str, SeqIO.SeqRecord],
                 *,
                 markers: pd.DataFrame = None,
                 #engine: str = "MaxQuant",
                 filt: bool = True,
                 filtints: bool = False,
                 invertgroup: bool = False):
        # engine

        try:
            #self.keyword = re.complie('(?<={}\s)[^\s_-]*(?=[ _-]\d+)'.format("LFQ intensity"))
            self.keyword = "LFQ intensity"
            self.cols = self._columns(pros, self.keyword)
        except:
            #self.keyword = re.compile('[^\s_-]*(?=\s{})'.format("MaxLFQ Intensity"))
            self.keyword = "MaxLFQ Intensity"
            self.cols = self._columns(pros, self.keyword)

        # keywords
        self.name = tuple(self.cols.keys())
        # column names
        self.expr_ = self.cols[self.name[0]]
        try:
            self.ctrl_ = self.cols[self.name[1]]
        except:  # single group
            self.ctrl_ = []

        if invertgroup:
            self.expr_, self.ctrl_ = self.ctrl_, self.expr_

        self.nfracs = len(self.expr_)

        if filt:
            self.pros = pros.loc[self._filt(pros, filtints)]
        else:
            self.pros = pros

        self.peps = peps
        self.db = fasta
        self.markers = markers

    def __len__(self):
        return len(self.pros.index)

    def __iter__(self):
        for i in self.pros.index:
            # yield self.val(i)
            yield Protein(self, i)

    def __getitem__(self, index):
        return Protein(self, index)

# filter proteins

    def _filt(self, pros, filtints) -> List[bool]:
        # filtering contaminants and reverse sequences.
        pattern = re.compile('sp\|')
        # pattern=re.compile('(?<!REV__)sp\|')
        crit1 = list(pattern.match(i) != None for i in pros.index)

        # filtering proteins of which all intensities in any group are lower than MININT.
        if filtints:
            crit2 = (pros[self.expr_] > MININT).any(axis=1).to_list()
            crit3 = (pros[self.ctrl_] > MININT).any(axis=1).to_list()
            filtered = list(x and y and z
                            for x, y, z in zip(crit1, crit2, crit3))
        # filtering proteins of which all intensities are 0.
        else:
            crit2 = (pros[self.expr_ + self.ctrl_] > 0).any(axis=1).to_list()
            filtered = list(x and y for x, y in zip(crit1, crit2))

        return filtered

# detect columns with a certain keyword

    def _columns(self, df, keyword) -> Dict[str, List[str]]:
        # begin with keyword and end up with digits
        #pattern = re.compile('{}\s.*[ _-]\d+$'.format(keyword))
        # before which is keyword and after which are digits
        #prefix = re.compile('(?<={}\s)[^\s]*(?=[ _-]\d+)'.format(keyword))
        # prefix=re.compile('\w*(?=[\s_-]\d+$)')
        prefix = re.compile('\w*(?=[\s_-]\d+[\s_-]\d+$)')
        cols = {}
        for i in df.columns:
            # if pattern.match(i):
            if keyword in i:
                target = re.sub(' ?{} ?'.format(keyword), '', i)
                try:
                    #name = prefix.search(i).group()
                    name = prefix.match(target).group()
                except:
                    name = 0
                try:
                    cols[name].append(i)
                except:
                    cols[name] = []
                    cols[name].append(i)
        print(cols)

        if len(cols) != 2:
            raise ValueError("invalid group number!")
        s = 0
        for i in cols.values():
            if s != len(i) and s != 0:
                raise ValueError(
                    "The experiment and control have different fraction number!"
                )
            s = len(i)
            # if pros[i].isna().any(axis=None):
            #raise ValueError("Invalid values in dataframe!")
            # sorting columns by experiment number(for both cases: _1 and _01).
            i.sort(key=lambda x: int(re.split('[_-]', x)[-2]))
        return cols

    def ints(self):
        return self.pros[self.ctrl_ + self.expr_]

    def norm_ints(self):
        return self.ints().div(self.ints().max(axis=1), axis='index')

    def expr(self):
        return self.pros[self.expr_]

    def ctrl(self):
        return self.pros[self.ctrl_]

# validate protein name and convert it into Protein IDs

    def val(self, name: Union[str, int]) -> ID:
        if isinstance(name, int):
            name = self.pros.index[name]
            return ID(name)
        if name in self.pros.index:
            return ID(name)
        else:
            for i in self.pros.index:
                # UniProt id, e.g. "P42166"
                if '|%s|' % name in i:
                    return ID(i)
                # Fasta header, e.g. "sp|P42166|LAP2A_HUMAN"
                elif name in i.split(";"):
                    return ID(i)
            raise ValueError("Protein ID not found")

    def seq(self, row: ID):
        return self.db[row.split(';')[0]].seq

# retrieve all peptides of a protein

    def pepts(self, row: ID, cols: List[str] = []) -> pd.DataFrame:
        # row=self.val(row)
        if cols == []:
            cols = self.expr_+self.ctrl_
        try:
            df = self.peps.iloc[self.pros["Peptide IDs"][row].split(
                ';')][["Start position", "End position"] + cols]
        except:
            df = self.peps.loc[lambda df: df["Protein"]
                               == row][["Start", "End"]+cols]
        # start position is 1-based
        #df["Length"] = df["End position"] - df["Start position"] + 1
        return df

# return sequence representation matrices of a protein

    def seq_repr(self, row, repr='intensity', rasterize=True, normalized=False) -> np.ndarray:

        if repr == 'spectral_count':
            expr, ctrl = self._columns(self.peps, "Experiment").values()
        elif repr == 'intensity':
            expr, ctrl = self.expr_, self.ctrl_
        else:
            raise NameError(
                "Representation of heatmap must be 'intensity' or 'spectral_count'.")

        df = self.pepts(row, expr+ctrl)
        l = len(self.seq(row))
        # rasterize sequence length to 100, in order to control time complexity and compatible to CV processing
        if rasterize:
            coef = 100 / l
            spix = (coef * df["Start position"]).round().astype(np.int8)
            epix = (coef * df["End position"]).round().astype(np.int8) + 1
        # caution: if False, rendering might be wrong
        else:
            coef = 1
            spix = df["Start position"].astype(np.int8) - 1
            epix = df["End position"].astype(np.int8)

        mx = np.zeros(
            (2, int(coef * l), self.nfracs))  # default: (2,100,nfracs)
        for i in df.index:
            mx[0, spix[i]:epix[i], :] += df.loc[i, expr].fillna(0).to_numpy()
            mx[1, spix[i]:epix[i], :] += df.loc[i, ctrl].fillna(0).to_numpy()
        # normalize each pixel comlumn with maximium over all fractions
        if normalized:
            # denominator: (1,100,1)
            mx = mx/(mx.max(axis=(0, 2), keepdims=True)+1)
        return mx

    def get_prot(self, row):
        return Protein(self, row)

    def to_list(self):
        return list(x.split(";")[0] for x in self.pros.index)

    def to_idlist(self):
        return list(x.split(";")[0].split("|")[1] for x in self.pros.index)


class Protein():
    '''
    Class for a protein(group).
    '''

    def __init__(self, proteingroup: ProteinGroups, name: str):
        self.g = proteingroup
        #self.name = name
        self.name = proteingroup.val(name)
        #self.m=mass(self.name, self.g.db)
        #self.seq=seq(self.name, self.g.db)
        self.seq = self.g.seq(self.name)
        self.l = len(self.seq)
        self.m = ProteinAnalysis(str(self.seq)).molecular_weight()/1000
        self.maxint = self.ints().max()

    def __repr__(self):
        return self.name.split(";")[0]

    def id(self):
        return str(self).split("|")[1]

    def minor_prot(self):
        return self.name.split(";")[1:]

    def ints(self):
        return self.g.ints().loc[self.name]

    def ctrl(self):
        return self.g.ctrl().loc[self.name]

    def expr(self):
        return self.g.expr().loc[self.name]

    def norm_ints(self):
        return self.ints() / self.maxint

    def pepts(self):
        return self.g.pepts(self.name)

    def seq_repr(self, **kwargs):
        return self.g.seq_repr(self.name, **kwargs)
