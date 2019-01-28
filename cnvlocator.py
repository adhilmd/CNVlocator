#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  8 12:20:30 2018

@author: Mohamood Adhil
"""
import sys
import pandas as pd
import numpy as np
from sklearn import preprocessing
import argparse
import pysal
import scipy as sp
from scipy import signal
from scipy.stats import nbinom

def argparser():
    parser = argparse.ArgumentParser(prog = '''\ncnvlocator.py''', description='''\n----------Calls CNV by comparing test and control sample--------- \n
    \n[Date: 7th March 2018], \n[help: python cnvlocator.py -h]\n''', usage = 'cnvlocator.py *args')
    parser.add_argument('-t','--test', type=str, dest='t', help="Path for Test interval file containig read accumlation (Mandatory)", action = 'store', required = True)
    parser.add_argument('-r','--ref', type=str, dest='r', help="Path for Refernece interval file containig read accumlation (Mandatory)", action = 'store', required = True)
    parser.add_argument('-lfc','--logfc', type=float, dest='lfc', help="Log fold change cut off Default = 0.6", action = 'store', required = True, default = 0.6)
    parser.add_argument('-mn','--minp', type=float, dest='mn', help="minimum percentage of reads in both sample in fraction Deafult = 0.05", action = 'store', default = 0.05)
    parser.add_argument('-ws','--windowskip', type=int, dest='ws', help="Maximum windows to skip if there are possible continuous CNV calls default=2", action = 'store', default = 2)
    parser.add_argument('-len','--minlength', type=int, dest='lm', help="minimum length of CNV default=length of the window size", action = 'store', default = 0)
    parser.add_argument('-p','--prefix', type=str, dest='pf', help="Output prefix (Mandatory)", action = 'store', required = True)
    parser.add_argument('-o','--outdir', type=str, dest='ot', help="outdir (Mandatory)", action = 'store', required = True)
    args = parser.parse_args()
    return(args)

def cnvcalls(sc,gapvalue):
    gp = gapvalue+1
    allte = sc
    cnt = 1
    gte = []
    mg = []
    for j in range(0, len(allte)):
        if (j < len(allte)-1 and cnt == 1):
            gte.append(allte[j])
            cnt = 0
        if (j < len(allte)-1 and cnt == 0):
            diff = allte[j+1] - allte[j]
            if (diff <= gp):
                gte.append(allte[j+1])
            if (diff > gp):
                cnt = 1
                mi = min(gte)
                ma = max(gte)
                te = [i for i in range(mi, ma+1)]
                mg.append(te)
                gte = list()
    return(mg)

def mergecalls(dfall,mc,cnvtype):
    chrc = []
    start = []
    end = []
    lfc = []
    pvaltest = []
    pvalref = []
    for i in mc:
        te = dfall.loc[min(i):max(i)]
        cr = te.chr.unique().tolist()
        for k in cr:
            re = te[te.chr == k]
            chrc.append(re.chr.unique().tolist()[0])
            start.append(re.start.min())
            end.append(re.end.max())
            lfc.append(round(re.lfc.median(),2))
            pvaltest.append(round(re.testpval.min(),5))
            pvalref.append(round(re.refpval.min(),5))	
    dfcalls = pd.DataFrame({'Chromosome':chrc, 'Start':start, 'End':end, 'LogFoldChange':lfc, 'Type': [cnvtype]*len(chrc), 'Pvaltest':pvaltest, 'Pvalref':pvalref})
    dfcalls = dfcalls[['Chromosome','Start','End','LogFoldChange','Type','Pvaltest','Pvalref']]
    return(dfcalls)

def main():
    #Reading Intervals file
    print (args)
    test = pd.read_csv(args.t, sep="\t", header=None)
    test.columns = ["chr", "start", "end", "treads", "tRPFHT", "tRPFHT-med"]
    ref = pd.read_csv(args.r, sep="\t", header=None)
    ref.columns = ["chr", "start", "end", "rreads", "rRPFHT", "rRPFHT-med"]
    #Merging both files with intervals
    df = test.merge(ref,on=["chr","start","end"])
    df['te'] = df['tRPFHT']/df['rRPFHT']
    df['lfc'] = np.log2(sp.signal.medfilt(df['te'],3))
    nreads = df[['tRPFHT-med','rRPFHT-med']].as_matrix()
    min_max_scaler = preprocessing.MinMaxScaler()
    normnp = min_max_scaler.fit_transform(nreads)
    normdata = pd.DataFrame(data=normnp,columns=['tnorm','rnorm'])
    dfall = df.merge(normdata,right_index=True,left_index=True)
    #Filtering    
    ampg = dfall[(dfall.lfc >= args.lfc) & ((dfall.tnorm >= args.mn) | (dfall.rnorm >= args.mn))]
    delg = dfall[(dfall.lfc <= -args.lfc) & ((dfall.tnorm >= args.mn) | (dfall.rnorm >= args.mn))]
    scg = ampg.index.tolist()
    scl = delg.index.tolist()
    ig = cnvcalls(scg,args.ws)
    il = cnvcalls(scl,args.ws)
    #negative binomial stats
    rpval = 1/(dfall['rRPFHT'].mean()+1)
    tpval = 1/(dfall['tRPFHT'].mean()+1)
    dfall['refpval'] = nbinom.cdf(dfall['rRPFHT'].tolist(), 1, rpval)
    dfall['testpval'] = nbinom.cdf(dfall['tRPFHT'].tolist(), 1, tpval)
    # Amplification calls
    cnvtype = "Amplification"
    dfamp = mergecalls(dfall,ig,cnvtype)
    # Deletion calls
    cnvtype = "Deletion"
    dfdel = mergecalls(dfall,il,cnvtype)
    frames = [dfamp, dfdel]
    result = pd.concat(frames)
    result = result[['Chromosome','Start','End','LogFoldChange','Type','Pvalref','Pvaltest']]
    result['Length'] = result['End'] - result['Start']
    if (args.lm != 0):
        result = result[(result.Length >= args.lm)]
    dfall = dfall[["chr","start","end","treads","tRPFHT","tRPFHT-med","rreads","rRPFHT","rRPFHT-med","lfc","tnorm","rnorm","refpval","testpval"]]
    outfile = args.ot+"/"+args.pf+".alldat"
    dfall.to_csv(outfile,sep="\t",header=True,index=False)
    outfile = args.ot+"/"+args.pf+".cnvcalls"
    result.to_csv(outfile,sep="\t",header=True,index=False)
    return 0

if (__name__ == "__main__"):
    args = argparser()
    status = main()
