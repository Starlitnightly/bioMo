import gseapy as gp
import pandas as pd 
import numpy as np 
import matplotlib.pyplot as plt

class enrichr(object):

    def __init__(self,gene_list,organism,title):
        self.gene_list=gene_list
        self.organism=organism
        self.title=title

    def lazy(self):
        self.enrichr_kegg()
        self.enrichr_go_bio()
        self.enrichr_go_cell()
        self.enrichr_go_mol()
        self.enrichr_wiki()

        self.enrichr_log(self.kegg_enr)
        self.enrichr_log(self.go_bio_enr)
        self.enrichr_log(self.go_cell_enr)
        self.enrichr_log(self.go_mol_enr)
        self.enrichr_log(self.wiki_enr)

        self.enrichr_plot(self.kegg_enr,'{0}_KEGG'.format(self.title),path='{0}/kegg_plot.png'.format(self.title))
        self.enrichr_plot(self.go_bio_enr,'{0}_GO_Bio'.format(self.title),path='{0}/go_bio_plot.png'.format(self.title))
        self.enrichr_plot(self.go_cell_enr,'{0}_GO_Cell'.format(self.title),path='{0}/go_cell_plot.png'.format(self.title))
        self.enrichr_plot(self.go_mol_enr,'{0}_Go_Mol'.format(self.title),path='{0}/go_mol_plot.png'.format(self.title))
        self.enrichr_plot(self.wiki_enr,'{0}_Wiki'.format(self.title),path='{0}/wiki_plot.png'.format(self.title))

    def enrichr_kegg(self):
        print('......Calculate kegg enrichr')
        if self.organism=='human':
            gene_set='KEGG_2021_Human'
        elif self.organism=='mouse':
            gene_set='KEGG_2019_Mouse'
        elif self.organism=='zebrafish':
            gene_set='KEGG_2019'
        enr = gp.enrichr(gene_list=self.gene_list,
                 gene_sets=[gene_set],
                 organism=self.organism,
                 description=self.title,
                 outdir='{0}/{1}'.format(self.title,gene_set),
                 no_plot=True,
                 cutoff=0.5 # test dataset, use lower value from range(0,1)
                )
        self.kegg_enr=enr
        return enr.res2d
    
    def enrichr_wiki(self):
        print('......Calculate wiki enrichr')
        if self.organism=='human':
            gene_set='WikiPathways_2019_Human'
        elif self.organism=='mouse':
            gene_set='WikiPathways_2019_Mouse'
        elif self.organism=='zebrafish':
            gene_set='WikiPathways_2018'
        enr = gp.enrichr(gene_list=self.gene_list,
                 gene_sets=[gene_set],
                 organism=self.organism,
                 description=self.title,
                 outdir='{0}/{1}'.format(self.title,gene_set),
                 no_plot=True,
                 cutoff=0.5 # test dataset, use lower value from range(0,1)
                )
        self.wiki_enr=enr
        return enr.res2d
    
    def enrichr_go_bio(self):
        print('......Calculate go bio enrichr')
        if self.organism=='human':
            gene_set='GO_Biological_Process_2021'
        elif self.organism=='mouse':
            gene_set='GO_Biological_Process_2021'
        elif self.organism=='zebrafish':
            gene_set='GO_Biological_Process_2018'
        enr = gp.enrichr(gene_list=self.gene_list,
                 gene_sets=[gene_set],
                 organism=self.organism,
                 description=self.title,
                 outdir='{0}/{1}'.format(self.title,gene_set),
                 no_plot=True,
                 cutoff=0.5 # test dataset, use lower value from range(0,1)
                )
        self.go_bio_enr=enr
        return enr.res2d

    def enrichr_go_mol(self):
        print('......Calculate go mol enrichr')
        if self.organism=='human':
            gene_set='GO_Molecular_Function_2021'
        elif self.organism=='mouse':
            gene_set='GO_Molecular_Function_2021'
        elif self.organism=='zebrafish':
            gene_set='GO_Molecular_Function_2018'
        enr = gp.enrichr(gene_list=self.gene_list,
                 gene_sets=[gene_set],
                 organism=self.organism,
                 description=self.title,
                 outdir='{0}/{1}'.format(self.title,gene_set),
                 no_plot=True,
                 cutoff=0.5 # test dataset, use lower value from range(0,1)
                )
        self.go_mol_enr=enr
        return enr.res2d

    def enrichr_go_cell(self):
        print('......Calculate go cell enrichr')
        if self.organism=='human':
            gene_set='GO_Cellular_Component_2021'
        elif self.organism=='mouse':
            gene_set='GO_Cellular_Component_2021'
        elif self.organism=='zebrafish':
            gene_set='GO_Cellular_Component_2018'
        enr = gp.enrichr(gene_list=self.gene_list,
                 gene_sets=[gene_set],
                 organism=self.organism,
                 description=self.title,
                 outdir='{0}/{1}'.format(self.title,gene_set),
                 no_plot=True,
                 cutoff=0.5 # test dataset, use lower value from range(0,1)
                )
        self.go_cell_enr=enr
        return enr.res2d

    def enrichr_log(self,enr):
        print('......Calculate log')
        enr.res2d['logp']=-np.log(enr.res2d['P-value'])
        enr.res2d['logadjp']=-np.log(enr.res2d['Adjusted P-value'])
        enr.res2d['logc']=np.log(enr.res2d['Combined Score'])
        return enr

    def enrichr_plot(self,enr,title,num=10,adjp=False,figsize=(4,6),fontsize=12,path=None):
        print('......Plot {0}'.format(title))
        pp=plt.figure(figsize=figsize)
        vol_c=['#3B8989','#a1d8d3','#c5e8e2']
        if adjp==True:
            kegg_result=enr.res2d.sort_values('logadjp',ascending=False)[:num]
        else:
            kegg_result=enr.res2d.sort_values('logp',ascending=False)[:num]
        #用ax控制图片
        ax=pp.add_subplot(1,1,1)
        if adjp==True:
            plt.barh(kegg_result['Term'],
                kegg_result['logadjp'],color=vol_c[0])
            plt.xlabel('-log(adj_pvalue)',fontsize=fontsize) 
        else:
            plt.barh(kegg_result['Term'],
                    kegg_result['logp'],color=vol_c[0])
            plt.xlabel('-log(pvalue)',fontsize=fontsize)                   
        
        plt.title(title,fontsize=fontsize)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(True)
        ax.spines['left'].set_visible(True)
        if path==None:
            plt.savefig("{0}.png".format(title),dpi=300,bbox_inches = 'tight')
        else:
            plt.savefig(path,dpi=300,bbox_inches = 'tight')
        return ax