# import required packages
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import os
import glob
import subprocess
from StringIO import StringIO
from scipy.optimize import curve_fit

# get paths to reference and bams and sample names
ref_index = glob.glob('/lustre/data/ANC/NGS/ref/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fai')[0]
bams = {}
exps = ['20190516_GA12']
for exp in exps:
    bams[exp] = glob.glob('/lustre/data/ANC/NGS/'+exp+'/bam_pipeline/*.bam')

print bams

def get_coverage(bams, ref_index, blocksize=10000, checkbams=False, outfile=None):
    """Generate coverage from bam files."""
    bams_tested = bams[:]
    # removing bams that cannot be uncompressed -- not ideal, could replace it at some point
    if checkbams:
        for i,bam in zip(range(len(bams)), bams):
            try:
                print 'checking file '+bam+'...'
                subprocess.check_output('gunzip -t '+bam, shell=True)
            except:
                print 'error in file '+bam+'\nwill be removed from list of .bam files'
                bams_tested.pop(i)
                continue 
            print 'valid'
    bam_arg = ' '.join(bams_tested)
    file_names = [os.path.basename(b) for b in bams_tested]
    sample_names = [f.split('.')[0] for f in file_names]
    # create empty data frame to collect results
    df = pd.DataFrame(columns=['chr','pos']+sample_names)
    # read info about chromosomes from ref index file
    df_index = pd.read_csv(ref_index, sep='\t', names=['chr', 'length'], usecols=[0,1])
    # sort by chromosome length
    df_index = df_index.sort_values(by='length', ascending=False)
    # get depth for all samples by looping over chromosomes
    for i,r in df_index.iterrows():
        chr_i = r.chr
        length_i = r.length
        cmd = 'samtools depth '+ bam_arg +' -r '+chr_i
        try:
            depth_i_str = subprocess.check_output(cmd, shell=True)
        except subprocess.CalledProcessError:
            print 'Samtools produced an error processing chromosome '+chr_i+'. Consider re-running with check_bams=True'
            raise
        depth_i_file = StringIO(depth_i_str) # enable samtools output to be read as file by pandas
        df_i = pd.read_csv(depth_i_file, sep='\t', names=['chr', 'pos']+sample_names) # data frame for chr i     
        blocks_i = pd.Series(np.repeat(range(int(length_i/blocksize)), blocksize)) # divide chr into blocks 
        blocks_i.append(pd.Series(np.repeat(blocks_i.iloc[-1],length_i-len(blocks_i))), ignore_index=True) # assign leftover to last block
        grouped = df_i.groupby(by=blocks_i) # group by blocks
        df_i_med = grouped.agg('median') # calculate blockwise median (also for position)
        df_i_med.insert(loc=0, column='chr', value=chr_i) # insert chr column
        df = df.append(df_i_med, ignore_index=True) # append to results data frame
    if outfile is not None:
        df.to_csv(outfile, sep='\t') # save results as csv
    return df


for exp in exps:
    print 'getting coverage for experiment '+exp
    df = get_coverage(bams[exp], ref_index, checkbams=True, outfile='/lustre/home/fgross/sequencing_nongit/coverage_analysis/results/coverage_'+exp+'.csv')

def correct_length_bias(df, drop_mito=True):
    df_corr = df.copy()
    if drop_mito:
        df_corr = df_corr[df_corr.chr!='Mito']
    sample_names = df_corr.columns[2::]
    for s in sample_names:
        df_corr[s] = df_corr[s].div(df_dropped_mito[s].median())
        x = df_corr[['chr','pos',s]].groupby(by='chr', sort=False).count().pos
        y = df_corr[['chr','pos',s]].groupby(by='chr', sort=False).agg('median')[s]
        popt, pcov = curve_fit(exp_fit, x, y, bounds=(0, [100., 1., 50.]))
        y_fit = exp_fit(x, *popt) - popt[2]
        df_corr[s] = df_corr[s].sub(df_corr.chr.map(y_fit))
    return df_corr

def correct_smile(df, cutoff=0.5):
    df_norm = df.copy()
    df_norm['norm_pos'] = df_norm.groupby('chr').transform(lambda x: np.linspace(0,1,x.count())).pos
    for s in sample_names:
        df_s = df_norm[np.abs(df_norm[s]-1)<cutoff] # remove outliers
        x = df_s.sort_values(by='norm_pos').norm_pos
        y = df_s.sort_values(by='norm_pos')[s]
        popt, pcov = curve_fit(quad_fit, x, y, bounds=(0, [1., 1., 1.]))
        y_fit = quad_fit(df_norm.norm_pos, *popt)-popt[-1]
        df_norm[s] = df_norm[s].sub(y_fit)
    df_norm.drop('norm_pos',axis=1,inplace=True)
    return df_norm


def quad_fit(x, x0, a, b):
    return a * (x-x0)**2 + b


def exp_fit(x, a, b, c):
    return a * np.exp(-b * x) + c
