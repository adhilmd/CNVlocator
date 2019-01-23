#Author: Mohamood Adhil, Date:15/10/2017
import argparse
from pybedtools import BedTool  
import pandas as pd
import glob, os

def argparser():
    parser = argparse.ArgumentParser(prog = '''\nintervals_sw.py''', description='''\n----------Sliding windows with read accumlation and fold change--------- \n
    \n[Date: 7th March 2018], \n[help: python intervals_sw.py -h]\n''', usage = 'intervals_sw.py *args')
    parser.add_argument('-w','--intervalfile', type=str, dest='w', help="File containing chromosome with size, two coloumn chromosome and size (Optional)", action = 'store', default="NA")
    parser.add_argument('-m','--mainpath', type=str, dest='m', help="mainpath where the bedgraph file and count files are present (Mandatory)", action = 'store', required = True)
    parser.add_argument('-f','--fname', type=str, dest='fn', help="Comma seperated prefix input file names (Mandatory)", action = 'store', required = True)
    parser.add_argument('-wn','--window', type=int, dest='wn', help="window size", action = 'store', default = 200)
    parser.add_argument('-s','--slide', type=int, dest='s', help="sliding window default=50", action = 'store', default = 50)
    parser.add_argument('-p','--prefix', type=str, dest='pf', help="Comma seperated output prefix (Mandatory). Same number of items as -f argument", action = 'store', required = True)
    parser.add_argument('-o','--outdir', type=str, dest='ot', help="outdir (Mandatory)", action = 'store', required = True)
    args = parser.parse_args()
    return(args)

def main():
    print (args)
    mainpath=args.m
    fnam = args.fn.split(',')
    bdg = [args.m + '/' + s + '.bdg' for s in fnam]
    cnt = [args.m + '/' + s + '.cnt' for s in fnam]
    pref = args.pf.split(',')
    #make windows with sliding windows
    if args.w != "NA":
        a = BedTool()
        windows = a.window_maker(b=args.w, w=args.wn, s=args.s)
    for z in range(0, len(pref)):
        if args.w == "NA":
            cntf = pd.read_csv(cnt[z], sep="\t", skiprows=1, header=None)
            cntf = cntf[:-1]
            cntf = cntf[[0,1]]
            outfile = args.ot+"/"+pref[z]+".temp"
            cntf.to_csv(outfile,sep="\t",header=False,index=False)
            a = BedTool()
            windows = a.window_maker(g=outfile, w=args.wn, s=args.s)
            os.remove(outfile)
        b=BedTool(bdg[z])
        #Map the read accumlation in the interval
        mp = windows.map(b,c=4,o="median")
        df = pd.read_table(mp.fn, names=['chrom', 'start', 'stop', 'score'])
        df['score'] = df['score'].replace(0, 0.25)
        cntf = pd.read_csv(cnt[z], sep="\t", skiprows=1, header=None)
        cntf.columns = ["chrom", "t", "actreads", "n"]
        cntf = cntf[["chrom","actreads"]]
        df = df.merge(cntf, how="inner",on="chrom")
        #Reads Per Fragment normalized to Hundered Thousand (RPFHT)  
        df['RPFHT'] = df['score']/(df['actreads']/100000)
        te = df.groupby(df.chrom)[['RPFHT']].median()
        te1 = pd.DataFrame({'chrom':te.index, 'Median':te.RPFHT})
        df = df.merge(te1,how="inner",on="chrom")
        df['RPFHT-Median'] = df['RPFHT']/df['Median']
        df['difference'] = df["stop"] - df['start']
        df = df[df.difference >= args.wn]
        df = df[["chrom","start","stop","score","RPFHT","RPFHT-Median"]]
        outfile = args.ot+"/"+pref[z]+".int"
        df.to_csv(outfile,sep="\t",header=False,index=False)
    return 0

if (__name__ == "__main__"):
    args = argparser()
    status = main()
    for f in glob.glob("/tmp/pybedtools*"):
        os.remove(f)
