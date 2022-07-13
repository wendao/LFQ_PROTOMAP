#!/usr/bin/env python
# -*- encoding: utf-8 -*-
'''
@file: mapping.py
@comment: 
@created: 2022/04/14 17:28:42
@auther: Zhu, Yinyue
@version: 1.0
'''

from mpld3 import plugins
import mpld3
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import ScalarFormatter, PercentFormatter
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from .proteins import ProteinGroups, Protein
from .fileio import file_to_set

import os
import re
import base64
import numpy as np
from typing import Union
from io import BytesIO

import matplotlib
matplotlib.use('Agg')
#from matplotlib.backends.backend_agg import FigureCanvasAgg

font = dict(size=12, family='Arial')
# filtering fractions of where intensity is higher than CUTOFF*max intensity among all fractions and groups
CUTOFF = 0.01


class ProtoMap():
    '''
    Class for plotting configuration.
    save_path: directory where saves images.
    normalize_intensity: if True, intensities will be normalized where maximum = 1 and represented as percent abundances.
    invert_fraction: if True, fraction numbers start at light proteins(long migrations on the gel).
    side_axis: plot a side axis identifying migrations of markers.
    heatmap: type of peptides mapping, 3 options: None, "intensity" and "spectral_count".
    interactive: if True, render interactive 3d.js images in the browser.
    legend: manually set group legends, List[str].
    dpi: resolution of output images.
    '''

    def __init__(
        self,
        save_path: str = None,
        *,
        normalize_intensity=False,
        invert_fraction=False,
        side_axis=False,
        heatmap: Union[str, None] = None,
        interactive=False,
        legend=None,
        dpi=300,
    ):
        self.save = save_path
        self.norm = normalize_intensity
        self.invert = invert_fraction
        self.sideax = side_axis
        self.heatmap = heatmap
        self.inter = interactive
        self.legend = legend
        self.dpi = 100
        
        if self.save != None and os.path.isdir(self.save) == False:
            os.mkdir(self.save)

# generate figure and save image

    def plot(self, protein: Protein):
        fig = Figure(figsize=(10, 5))
        ax1, ax2 = fig.subplots(
            1,
            2,
            sharey=True,
            #figsize=(10, 5),
            gridspec_kw=dict(width_ratios=[2, 1], wspace=0.05))

        if self.heatmap == "intensity":
            self.hmaplot(protein, ax1, spec_count=False)
        elif self.heatmap == "spectral_count":
            self.hmaplot(protein, ax1, spec_count=True)
        else:
            self.peplot(protein, ax1)
        self.proplot(protein, ax2)
        # static image
        if self.save != None:
            try:
                path = os.path.join(self.save,
                                    re.sub('\|', '_', str(protein)) + '.png')
                fig.savefig(path, dpi=self.dpi)
            except:
                print("Save failed.")
            return os.path.abspath(path)
        # interactive d3.js encoding
        elif self.inter:
            return mpld3.fig_to_html(fig)
        # base64 encoding
        else:
            buffer = BytesIO()
            fig.savefig(buffer, dpi=300)
            img = base64.b64encode(buffer.getvalue()).decode()
            buffer.close()
            return img
            # return f"<img src='data:image/png;base64,{img}' width='800'/>"

        # plt.show()
        # plt.close("all")

# plot protein map

    def proplot(self, prot: Protein, ax: Axes):
        if self.norm:
            height_e = prot.expr() / prot.maxint
            height_c = prot.ctrl() / prot.maxint
            ax.xaxis.set_major_formatter(PercentFormatter(xmax=1.0))
            ax.set_xlabel(prot.g.keyword + " normalized", fontdict=font)
        else:
            height_e = prot.expr()
            height_c = prot.ctrl()
            # ax.xaxis.set_major_formatter(ScalarFormatter())
            ax.set_xlabel(prot.g.keyword, fontdict=font)

        yticks = np.arange(prot.g.nfracs) + 1
        ax.barh(yticks - 0.1, height_e, height=0.2, color='red')  # experiment
        ax.barh(yticks + 0.1, height_c, height=0.2, color='blue')  # control
        # ax.grid(axis="y",color='gray')

        ax.set_title("Mass: {:.1f} kDa    Length:{:d}".format(prot.m, prot.l),
                     fontdict=font)
        if isinstance(self.legend, list):
            legend = self.legend
        else:
            legend = prot.g.name
        ax.legend(legend, loc="upper right")
        # marker ticks
        if self.sideax:
            try:
                ax_ = ax.twinx()
                # ax_.invert_yaxis()
                # ax_.set_ylim(0,self.nfracs+1)
                ax_.set_ylim(prot.g.nfracs + 1, 0)
                ax_.set_yticks(np.array(prot.g.markers["f"]))
                ax_.set_yticklabels(prot.g.markers["m"])
                ax_.set_ylabel("Molecular weight (kDa)", fontdict=font)
            except:
                pass
        return

# plot peptides map

    def peplot(self, prot: Protein, ax: Axes):
        peptides = prot.pepts()
        expr = prot.expr()
        ctrl = prot.ctrl()
        # only plotting peptides map for fractions where intensity is visible on protein-level map.
        for i in range(prot.g.nfracs):
            if expr[i] > CUTOFF * prot.maxint:
                # e = peptides.loc[lambda df: df[prot.g.expr_[i]] > 0][
                #    ["Start position", "End position"]].to_numpy()
                # 0:2 start and end columns
                e = peptides.loc[lambda df: df[prot.g.expr_[i]]
                                 > 0].iloc[:, 0:2].to_numpy()
                # broken_barh() accepts arguments as (start,length)
                e[:, 1] = e[:, 1] - e[:, 0] + 1
                ax.broken_barh(e, (i + 1, 0.2), color="blue")
            if ctrl[i] > CUTOFF * prot.maxint:
                c = peptides.loc[lambda df: df[prot.g.ctrl_[i]]
                                 > 0].iloc[:, 0:2].to_numpy()
                c[:, 1] = c[:, 1] - c[:, 0] + 1
                ax.broken_barh(c, (i + 1 - 0.2, 0.2), color="red")

        yticks = np.arange(prot.g.nfracs) + 1
        ax.set_xlim(1, prot.l + 1)
        ax.set_ylim(0, prot.g.nfracs + 1)
        ax.set_yticks(yticks)
        # ax.set_xticks(np.arange(0,l,100))
        if self.invert:
            ax.set_yticklabels(prot.g.nfracs + 1 - yticks)
        else:
            ax.invert_yaxis()
        ax.set_title(prot.name.split(";")[0], fontdict=font)
        ax.set_ylabel("Fraction", fontdict=font)
        ax.set_xlabel("Sequence position", fontdict=font)
        return


# plot peptides heatmap

    def hmaplot(self, prot: Protein, ax: Axes, **kwargs):
        arr = prot.seq_repr(**kwargs)  # (2,l,nfracs)
        _, l, f = arr.shape
        v1 = np.insert(arr[0].T, np.arange(f), np.zeros(l),
                       axis=0)  # (2*nfracs,100)
        v2 = np.insert(arr[1].T, np.arange(f) + 1, np.zeros(l), axis=0)
        # grid location of matrix elements
        X, Y = np.meshgrid(
            np.arange(l + 1) + 1,
            np.arange(0, f + 0.5, 0.5) + 0.5)  # (43,101)
        p1 = ax.pcolormesh(X, Y, v1, vmax=arr.max(), cmap="Blues")
        p2 = ax.pcolormesh(X, Y, v2, vmax=arr.max(), cmap="Reds", alpha=0.5)

        ax.set_xlim((1, prot.l))
        xticks = ax.get_xticks()  # get l-scale ticks
        yticks = np.arange(prot.g.nfracs) + 1
        ax.set_xticks(100 / prot.l * xticks)
        ax.set_xticklabels(xticks.astype(int))
        ax.set_xlim((1, 100))
        ax.set_yticks(yticks)
        if self.invert:
            ax.set_yticklabels(prot.g.nfracs + 1 - yticks)
        else:
            ax.invert_yaxis()
        ax.set_title(prot.name.split(";")[0], fontdict=font)
        ax.set_ylabel("Fraction", fontdict=font)
        ax.set_xlabel("Sequence position", fontdict=font)

        # colorbars
        divider = make_axes_locatable(ax)
        c1 = divider.append_axes("left", size="5%", pad=0.5)
        c2 = divider.append_axes("left", size="5%", pad=0.35)

        ax.figure.colorbar(p1, cax=c1)
        ax.figure.colorbar(p2, cax=c2)
        c1.yaxis.set_ticks_position('left')
        c2.yaxis.set_ticks_position('left')
        # ax.figure.colorbar(p1, ax=ax,#label=prot.g.name[1],
        #                  fraction=0.08,location='left',pad=0.1)
        # ax.figure.colorbar(p2, ax=ax,#label=prot.g.name[0],
        #                  fraction=0.08,location='left',pad=0.15)
        return


def scatterplot(x, y, proteins: ProteinGroups, screenlist, save):
    screened = file_to_set(screenlist)
    labels = proteins.to_idlist()
    colors = []
    targets = []

    for i in labels:
        if i in screened:
            colors.append('red')
        else:
            colors.append('blue')
        path = './{}.png'.format(re.sub('\|', '_', i))
        targets.append(path)

    fig = Figure()
    ax = fig.add_subplot(1, 1, 1)
    sct = ax.scatter(x, y, s=1, c=colors)
    tooltip = plugins.PointHTMLTooltip(sct, labels=labels, targets=targets)
    plugins.connect(fig, tooltip)
    if save:
        mpld3.save_html(fig, save)
    mpld3.show()
    # return mpld3.fig_to_html()