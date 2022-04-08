import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt
from adjustText import adjust_text
import seaborn as sns

from matplotlib.colors import LinearSegmentedColormap
#myColors = ((0.8, 0.0, 0.0, 1.0), (0.0, 0.8, 0.0, 1.0), (0.0, 0.0, 0.8, 1.0))
#cmap = matplotlib.colors.LinearSegmentedColormap.from_list("", colors)
colors=['#F7828A',"#F9C7C6","#FDFAF3","#D4E3D0","#9CCCA4",]
c = LinearSegmentedColormap.from_list('Custom', colors, len(colors))
colors.reverse()
c_r=LinearSegmentedColormap.from_list('Custom', colors, len(colors))
import os
import sys

class deg_desep2(object):

    def __init__(self,count_data,deseq2_data,meta_data,ensemble):

        data1=deseq2_data.copy()
        data1['sig'] = 'normal'
        data1.loc[(data1.log2FoldChange> 1 )&(data1.padj < 0.05),'sig'] = 'up'
        data1.loc[(data1.log2FoldChange< -1 )&(data1.padj < 0.05),'sig'] = 'down'
        data1['log(padj)'] = -np.log10(data1['padj']+0.000000001)
        data1['Abs_FC']=abs(data1['log2FoldChange'])
        data1=data1.sort_values('Abs_FC',ascending=False)

        self.count_data=count_data
        self.deseq2_data=data1
        self.ensemble=ensemble
        self.meta_data=meta_data
        self.vol_c=['#3B8989','#a1d8d3','#c5e8e2']
        #print(os.getcwd())
        #print(sys.path)
        curr_path=os.path.dirname(__file__)
        if ensemble=='danRer7':
            pair=pd.read_csv(curr_path+'/pair/pair_danRer7.tsv',sep='\t')
        elif ensemble=='danRer11':
            pair=pd.read_csv(curr_path+'/pair/pair_danRer11.tsv',sep='\t')
        elif ensemble=='GRCh37':
            pair=pd.read_csv(curr_path+'/pair/pair_GRCh37.tsv',sep='\t')
        elif ensemble=='GRCh38':
            pair=pd.read_csv(curr_path+'/pair/pair_GRCh38.tsv',sep='\t')
        elif ensemble=='GRCm39':
            pair=pd.read_csv(curr_path+'/pair/pair_GRCm39.tsv',sep='\t')
        pair.set_index(pair.columns[0],inplace=True)
        self.pair=pair

    def change_log2Fold_threshold(self,threshold):
        data1=self.deseq2_data.copy()
        data1['sig'] = 'normal'
        data1.loc[(data1.log2FoldChange> threshold )&(data1.padj < 0.05),'sig'] = 'up'
        data1.loc[(data1.log2FoldChange< (0-threshold) )&(data1.padj < 0.05),'sig'] = 'down'
        self.deseq2_data=data1


    def plot_hist(self):
        plt.hist(self.deseq2_data['log2FoldChange'])

    def plot_vol(self,title):
        font1={
            'size':15,
        }
        pp=plt.figure(figsize=(4,4))
        ax=pp.add_subplot(1,1,1)
        #data1=data1.sort_values
        data1=self.deseq2_data
        vol_c=self.vol_c

        plt.scatter(x=data1[data1['sig']=='up']['log2FoldChange'],y=data1[data1['sig']=='up']['log(padj)'],color=vol_c[0],label='up')
        plt.scatter(x=data1[data1['sig']=='down']['log2FoldChange'],y=data1[data1['sig']=='down']['log(padj)'],color=vol_c[1],label='down')
        plt.scatter(x=data1[data1['sig']=='normal']['log2FoldChange'],y=data1[data1['sig']=='normal']['log(padj)'],color=vol_c[2],label='normal')

        test=data1.index
        texts=[plt.text(data1.loc[i,'log2FoldChange'], data1.loc[i,'log(padj)'],
                self.pair.loc[i,'symbol'],fontdict={'size':12,'weight':'bold'}) for i in test[:10]]

        adjust_text(texts,only_move={'text': 'xy'},lim=5000,precision=0.0001,arrowprops=dict(arrowstyle='->',color='red',
                                    lw=1))

        plt.yticks(fontsize=15)
        plt.xticks(fontsize=15)
        plt.legend(['up:{0}'.format(len(data1[data1['sig']=='up'])),'down:{0}'.format(len(data1[data1['sig']=='down']))], \
            loc='center', bbox_to_anchor=(0.5, -0.3), ncol=2,fontsize=15)
        ax.set_ylabel('-log(padj)',font1)                                    
        ax.set_xlabel('log2FC',font1)
        ax.plot([data1['log2FoldChange'].min(),data1['log2FoldChange'].max()],[-np.log10(0.05),-np.log10(0.05)],linewidth=2, linestyle="--",color='black')
        ax.plot([2,2],[data1['log(padj)'].min(),12],linewidth=2, linestyle="--",color='black')
        ax.plot([-2,-2],[data1['log(padj)'].min(),12],linewidth=2, linestyle="--",color='black')
        plt.title(title,fontsize=15)
        plt.savefig("{0}_volcano.png".format(title),dpi=300,bbox_inches = 'tight')
        return ax

    def plot_heatmap(self,col_c,title):
        diff_gene=self.deseq2_data[self.deseq2_data['sig']!='normal'].index.values
        diff_count=self.count_data.loc[diff_gene]
        diff_count.index=self.pair.loc[diff_gene,'symbol']
        diff_meta=self.meta_data
        

        a=sns.clustermap(diff_count, 
                    cmap=c_r, 
                    standard_scale = 0,figsize=(4,12),
                    col_colors=diff_meta['Type'].map(col_c),
        )
        a.ax_heatmap.yaxis.set_tick_params(labelsize=12)
        a.ax_heatmap.xaxis.set_tick_params(labelsize=12)
        #a.ax_heatmap.xaxis.set_ticklabels(['','HUVEC','','','','HUVEC+ADSC'])
        labels=a.ax_heatmap.xaxis.get_ticklabels()
        #a.ax_heatmap.xaxis.set_text(['','RAW','','','RAW+EV',''])
        plt.setp(labels, rotation=45, horizontalalignment='right',)
        for label in diff_meta['Type'].unique():
            a.ax_col_dendrogram.bar(0, 0, color=col_c[label],
                                    label=label, linewidth=0)
        a.ax_col_dendrogram.legend(loc="best", ncol=1,bbox_to_anchor=(-0.5, 0., 0.5, 0.5),fontsize=12)
        a.cax.set_position([-.15, .2, .03, .45])
        plt.setp(a.cax.yaxis.get_majorticklabels(), fontsize=12)
        plt.savefig("{0}_heatmap.png".format(title),dpi=300,bbox_inches = 'tight')
        return a

        