#!/usr/bin/env python

import sys, os

script_path = os.path.dirname(__file__)
localdotplot_path = os.path.join(script_path, '../../mmfold/local_dotplot/')
print(localdotplot_path)
sys.path.append(localdotplot_path)
import local_dotplot_lib as ldp

from matplotlib import pyplot as plt
import seaborn as sns
import numpy as np
import argparse

import matplotlib.pyplot as plt
import matplotlib.colors as colors
import itertools

from Bio.SeqRecord import SeqRecord
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import Alphabet
from subprocess import Popen, PIPE
from os.path import isfile

def is_valid_file(file_name):
    if os.path.isfile(file_name):
        return os.path.abspath(file_name)
    else:
        raise FileNotFoundError(os.path.abspath(file_name))
def is_valid_directory(dir_name):
    if os.path.isdir(dir_name):
        return os.path.abspath(dir_name)
    else:
        raise FileNotFoundError(os.path.abspath(dir_name))

def is_valid_sequence(s):
    rec = SeqRecord(Seq(s, IUPAC.unambiguous_rna), id="RNA")
    if not Alphabet._verify_alphabet(rec.seq):
        raise RuntimeError("Invalid RNA sequence, unknown characters in input string {}".format(s))
    return rec
  
def makeSafeFilename(inputFilename):   
    try:
        safechars = string.letters + string.digits + " -_."
        return filter(lambda c: c in safechars, inputFilename)
    except:
        return ""  


def get_CD_record(full_rec, utr5_len, utr3_len):
    if utr3_len == 0:
        CD_seq = full_rec.seq[utr5_len:]
    else:
        CD_seq = full_rec.seq[utr5_len: -utr3_len]
    CD_rec = SeqRecord(CD_seq, id=full_rec.id+'-CDS', name=full_rec.name, description='utr5-0 utr3-0')
    return CD_rec


def write_dp_from_matrix(matrix, out_dp, template_dp,p_range=[-1,1]):
    from math import sqrt
    import re
    ureg = re.compile(r'^(\d+)\s+(\d+)\s+(\d+\.\d+)\s+[ul]box\s*')
    with open(template_dp) as template, open(out_dp, 'w') as out_handle:
        # Write first part before probs
        for line in template:
            if "ubox" in line or "lbox" in line:
                um = ureg.match(line)
                if um:
                    break
            out_handle.write(line)
        cx = matrix.tocoo()
        #Write alternative given probs
        for i,j,p in zip(cx.row, cx.col, cx.data):# used to have itertoolz.izip in Python2
            if abs(p) >= 1e-5 and p_range[0]<=p<=p_range[1]:
                out_handle.write("{} {} {:1.9f} ubox\n".format(i+1,j+1,sqrt(abs(p))))
        #Write remianing footer after probs
        for line in template:
            if "ubox" in line or "lbox" in line:
                um = ureg.match(line)
                if um:
                    continue
            out_handle.write(line)

def write_diff_dp(dp_wild, dp_mut, out_dp):
    '''Reads two dotplot matrices and write the absoloute difference of basepair probs into out_dp'''
    #     print "write_diff_dp: ", dp_wild, dp_mut, out_dp
    mut_matrix, mfe_dic = ldp.parse_dp_ps_sparse(dp_wild, sparse=True)
    wild_matrix, mfe_dic = ldp.parse_dp_ps_sparse(dp_mut, sparse=True)
    assert mut_matrix.shape == wild_matrix.shape
    diff_mat = (mut_matrix - wild_matrix)
    write_dp_from_matrix(diff_mat, out_dp=out_dp, template_dp=dp_wild)
    
    # dp_removed = '.'.join(
        # out_dp.split('.')[:-1]+['.removed.dp']+out_dp.split('.')[-1:])
    dp_removed = os.path.join(os.path.dirname(out_dp), os.path.basename(out_dp).replace('.ps','_removed.dp'))
    write_dp_from_matrix(diff_mat, out_dp=dp_removed, 
                         template_dp=dp_wild,p_range=[0,1])
    # dp_introduced = '.'.join(
    #     out_dp.split('.')[:-1]+['.introduced.dp']+out_dp.split('.')[-1:])
    dp_introduced = os.path.join(os.path.dirname(out_dp), os.path.basename(out_dp).replace('.ps','_introduced.dp'))

    write_dp_from_matrix(diff_mat, out_dp=dp_introduced, 
                         template_dp=dp_wild,p_range=[-1.0,0])
    return out_dp, dp_removed, dp_introduced 
    



def call_vienna_plfold(sequence, seq_name, do_localfold=False, local_W=200, local_L=150, out_dir='./'):
    '''Runs Vienna RNAfold with partition function for all sequences inside input fasta file 
    # call_RNAfold_pf("ACCGGCUUAAAGG", "seq1")'''

    from subprocess import Popen, PIPE
    dp_file_name = os.path.join(out_dir, "{}_dp.ps".format(seq_name))
    unp_file_name = os.path.join(out_dir, "{}_lunp".format(seq_name))
    if isfile(dp_file_name): # Caution Race condition may occur 
        os.remove(dp_file_name)
    if isfile(unp_file_name): # Caution Race condition may occur 
        os.remove(unp_file_name)
    if not do_localfold:
        local_W, local_L = len(sequence), len(sequence)
    # RNAFOLD = 'RNAfold -p2 '
    RNAPLFOLD = 'RNAplfold -W {} -L {} -u 1 '.format(local_W, local_L)  # -u 1 for unpaired probablitiy 
    assert len(sequence.split()) == 1
    cmd = "cd {}; ".format(out_dir)
    cmd += ('echo ">%s\\n%s\\n" | '%(seq_name, sequence))  
    cmd += RNAPLFOLD 
    print(cmd)
    p = Popen(cmd, stdin=PIPE, shell=True, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    print (out)
    print (err)
    has_error = False
    if err:
        if b"warn" in err.lower():
            print ("Warning in calling call_RNAfold:\n {} {}\n".format(out, err))
        else:
            raise RuntimeError("Error in calling call_RNAfold: {} {}\n".format(out, err))
    if not isfile(dp_file_name):
        raise RuntimeError("Error: Expected dp file: {} is not created!".format(dp_file_name))
    if not isfile(unp_file_name):
        raise RuntimeError("Error: Expected lunp file: {} is not created!".format(unp_file_name))
    
    return dp_file_name, unp_file_name

    
def run_dot2circ(dp_file, prefix, out_dir=""):
    file_name_string = "".join(x for x in dp_file if x not in ['|','<', '>'])
    #print ("dp_file:", dp_file, file_name_string)
    
    cmd = 'cd "{}/dot2circ/"; python dot2circ.py --outputdir \'{}\' --prefix \'{}\' --dp-file \'{}\' --title \'{}\' '.format(
        script_path, out_dir, prefix+'-circos', file_name_string, prefix)
    print (cmd)
    p = Popen(cmd, stdin=PIPE, shell=True, stdout=PIPE, stderr=PIPE)
    out, err = p.communicate()
    print (out.decode('ascii'))
#     print (err.decode('ascii'))
    has_error = False
    if err:
        raise RuntimeError("Error in calling dot2circ.py: {} {}\n".format(out, err))

    
def create_circos_annotation(CDS_len, utr5_len, utr3_len, snp_locs, snp_names):
    '''
    Genes formatted Example:
    seq 193 759 CDS fill_color=green,r0=1.01r,r1=1.01r+20p
    seq 0 192  5UTR fill_color=yellow,r0=1.01r,r1=1.01r+20p
    seq 760 5765  3UTR fill_color=yellow,r0=1.01r,r1=1.01r+20p
    seq 227 228 SN-86 fill_color=vdred,r0=1.03r,r1=1.03r+20p
    '''
    print ("snp_locs:", snp_locs)
    start = 0
    formatted_str = ""
    if utr5_len is not None:
        formatted_str += 'seq {} {} 5UTR fill_color=yellow,r0=1.01r,r1=1.01r+20p\n'.format(start, start+utr5_len)
        start += utr5_len + 1 # Tocheck: maybe plus one not needed?
    if CDS_len is not None:
        formatted_str += 'seq {} {} CDS fill_color=green,r0=1.01r,r1=1.01r+20p\n'.format(start, start+CDS_len)
        start += CDS_len + 1 # Tocheck: maybe plus one not needed?
    if utr3_len is not None:
        formatted_str += 'seq {} {} 3UTR fill_color=blue,r0=1.01r,r1=1.01r+20p\n'.format(start, start+utr3_len)
        start += utr3_len + 1 # Tocheck: maybe plus one not needed?
    if snp_locs is not None:
        for i in range(len(snp_locs)):
            formatted_str += 'seq {} {} c.-{} fill_color=vdred,r0=1.03r,r1=1.03r+20p\n'.format(snp_locs[i], snp_locs[i]+1, snp_names[i])
    # print (formatted_str)
    genes_file =  '{}/dot2circ/data/genes.formatted.txt'.format(script_path)
    with open (genes_file, 'w') as genes_out:
        genes_out.write(formatted_str)

def plot_circos_seq_annotate(rec, annotate_indices, annotate_names, local_fold=False, 
                             plotted_seq_lenght=None,utr5_l=0, utr3_l = 0,color='r',dp_full = None,suffix='' ):
    from Bio import SeqIO
 
    seq = rec.seq

    if len(annotate_indices) != len(annotate_names):
        raise RuntimeError("Mismatch indic/names: {}, {} }".format(len(annotate_indices), len(annotate_names))) 
    
    ID = rec.id
    if dp_full is None:
        dp_full, unp_full = call_vienna_plfold(rec.seq, ID, local_fold)
    create_circos_annotation(len(rec)-utr5_l-utr3_l, utr5_l, utr3_l, annotate_indices, annotate_names)

    run_dot2circ(dp_full, ID+suffix)
        
    #dp_diff = dp_mut.replace('.ps', '_diff.ps')
    #write_diff_dp(dp_full, dp_mut, dp_diff)
    #run_dot2circ(dp_diff, rec_mut.id+suffix+'-diff')
    

from cycler import cycler

plt.rc('axes', prop_cycle=(cycler('color', ['k',(0.8901960784313725, 0.10196078431372549, 0.10980392156862745)]
#                                   sns.color_palette()
#                                   [sns.color_palette("Paired")[1], sns.color_palette("Paired")[-1]]
                                  #                                   sns.color_palette("RdBu_r", 2)
#                                   ['k', 'b', 'r', 'g']
                                 ) 
#                            +                            cycler('linestyle', ['-', '--', ':', '-.']) 
                          ))

   

def get_unpaired_probs(unp_file):
    # Read Vienna RNAplfold unpaired prob file (-u 1) into dict
    with open(unp_file) as unp_in:
        line1 = unp_in.readline()
        if "#unpaired probabilities" not in line1:
            raise IOError('Unexpected header for lunp file: {}'.format(line1))
        line2 = unp_in.readline()
        if "#i$\tl=1" not in line2:
            raise IOError('Unexpected second header for lunp file: {}'.format(line2))
        up_dic = dict()
        for line in unp_in:
            splits = line.split()
            assert len(splits) >= 2
            pos = int(splits[0])
            up_prob = float(splits[1])
            assert pos >=1
            assert up_prob >= 0 and up_prob <= 1
            assert pos not in up_dic
            up_dic[pos] = up_prob
            
    
    return up_dic


def plot_up_dict(up_dic, plot_lims=None, title='XX', fig=None, diff=False,tidy=False):
    if plot_lims is None:
        x, y = up_dic.keys(), up_dic.values()
    else:
        x, y = up_dic.keys()[plot_lims[0]:plot_lims[1]+1], up_dic.values()[plot_lims[0]:plot_lims[1]+1]
    if fig is None:
        fig = plt.figure(figsize=(9, 3))
    x, y = list(x), list(y)
    ax = fig.add_subplot(111) 
    ax.plot(x, y, label=title)
    if not tidy:
        ax.legend(loc='upper left')#, bbox_to_anchor=(0.0, 1.1))


    import numpy as np
    minor_ticks = np.arange(min(x), max(x), 1)                                               
    major_ticks = np.arange(min(x)-min(x)%5, max(x), 5)                                               

    ax.set_xticks(minor_ticks, minor=True)                                           
    ax.set_xticks(major_ticks) 
    ax.set_xlabel('Position')
    if diff:
        ax.set_yticks([-1, 1], minor=True)                                           
        ax.set_ylim([-1.05,1.05])
        ax.set_ylabel('Accessibility(WT) - Accessibility(MUT)')
    else:
        ax.set_yticks([0,1], minor=True)                                           
        ax.set_ylim([-0.05,1.05])
        ax.set_ylabel('Accessibility')
#     ax.grid(which='both')                                                            

    # or if you want differnet settings for the grids:                               
    ax.grid(which='minor', alpha=0.5)
    ax.axhline(0)
#     ax.axhline(0, linestyle='--', color='k', alpha=0.5) # horizontal lines
#     ax.axhline(1, linestyle='--', color='k', alpha=0.5) # horizontal lines
    
    ax.set_xlim([min(x)-1, max(x)+1])
    if not tidy:
        ax.set_title(title)


def heatmap_up_dict(up_dic, plot_lims=None, title='XX', ax=None,fig=None, diff=False,
                    ticklabel=True,legend=True):
    if plot_lims is None:
        x, y = list(up_dic.keys()), list(up_dic.values())
    else:
        x, y = list(up_dic.keys())[plot_lims[0]:plot_lims[1]+1], list(up_dic.values())[plot_lims[0]:plot_lims[1]+1]
    
    
    if ax is None:
        ax = fig.add_subplot(121)
    ax.yaxis.set_label_position("right")

    # Now adding the colorbar
#     pos1 = ax.get_position() # get the original position 
#     pos2 = [pos1.x0 + 0.3, pos1.y0 + 0.3,  pos1.width / 2.0, pos1.height / 2.0]
    if legend:
        cax = fig.add_axes([1.0, 0.5, 0.1, .1])
    else:
        cax = None

#     ax.plot(x, y, label=title)
#     sns.pointplot(x,y,ax=ax)
    sns.heatmap(np.reshape(np.array(y),(len(y),1)), cmap=sns.cubehelix_palette(as_cmap=True), ax=ax, 
                cbar=legend,cbar_ax=cax,vmin=0.0, vmax=1.0,
                yticklabels=x, 
               )
#     cbar = ax.collections[0].colorbar
#     cbar.set_ticks([0., .2, .4, .6, .8, 1.0])
    
#     cbar.set_ticklabels(['low', '20%', '75%', '100%'])
    if legend:
        ax.legend(loc='upper left')#, bbox_to_anchor=(0.0, 1.1))

    minor_ticks = np.arange(min(x), max(x), 1)                                               
    major_ticks = np.arange(min(x)-min(x)%10, max(x), 10)                                               
    yticks = x
    keptticks = yticks[::10]
    yticks = ['' for y in yticks]
    yticks[::10] = keptticks
#     ax.set_yticks(minor_ticks,                   minor=True)
    from matplotlib.ticker import AutoMinorLocator,MultipleLocator,FormatStrFormatter,FuncFormatter,IndexLocator,StrMethodFormatter
    minorLocator = IndexLocator(1,offset=0.5)
    majorLocator = IndexLocator(10,offset=0.5)
    def incer(x, pos):
        'The two args are the value and tick position'
        return '%d' % (x+1)

    majorFormatter = FuncFormatter(incer)
#     majorFormatter = StrMethodFormatter('{x}',use_offset=False                                       )

#     ax.set_yticklabels(yticks, rotation=0,)

    ax.yaxis.set_minor_locator(minorLocator)
    ax.yaxis.set_major_locator(majorLocator)
    ax.tick_params(which='major', length=6)
    ax.tick_params(which='minor', length=3)
    if ticklabel:
        ax.yaxis.set_major_formatter(majorFormatter)
#         ax.set_yticklabels(yticks, rotation=0,)

#         ax.tick_params(which='minor', length=6)

        ax.yaxis.set_label_position("right")
    else:
        ax.set_yticklabels([])
#     ax.set_yticks(major_ticks) 
#     ax.set_ylabel('Position')
#     ax.set_yticks([-1, 1], minor=False, )                                           

    if diff:
        ax.set_xticks([-1, 1], minor=True)                                           
        ax.set_xlim([-1.05,1.05])
        ax.set_xlabel('P_unpaired(wild) - P_unpaired(mut)')
    else:
        ax.set_xticks([], minor=False,)                                           
#         ax.set_xlim([-0.05,1.05])
#         ax.set_xlabel('Accessibility')
#     ax.grid(which='both')                                                            

    # or if you want differnet settings for the grids:                               
#     ax.grid(which='minor', alpha=0.5) 
#     ax.axhline(0, linestyle='--', color='k', alpha=0.5) # horizontal lines
#     ax.axhline(1, linestyle='--', color='k', alpha=0.5) # horizontal lines
    
#     ax.set_ylim([min(x)-1, max(x)+1])
    ax.set_title(title,rotation=90,va='bottom')


def plot_unpaired_probs(up_file_pairs, plot_heatmap=False,rang=None, out_dir='./'):
    for up_file_wild, up_file_mut in up_file_pairs:
        # print (up_file_wild, up_file_mut)
        d_wild = get_unpaired_probs(up_file_wild)
        d_mut = get_unpaired_probs(up_file_mut)
        dict_diff = {f:(d_wild[f] - d_mut[f])  for f in d_wild}
        if plot_heatmap:
            title_key = 'heatband'
            fig, ax = plt.subplots(1,figsize=(0.5, 11))
            #heatmap_up_dict(dict_diff, rang, title=os.path.basename(up_file_wild).replace('-WT-','-').replace('_lunp','-DIFF'), ax=ax, fig=fig,diff=True)
        else:
            title_key = 'ECG'

            fig = plt.figure(figsize=(8, 2))
            plot_up_dict(dict_diff, rang, title=os.path.basename(up_file_wild).replace('-WT-','-').replace('_lunp','-DIFF')
                         , fig=fig, diff=True,tidy=True)
        fig.savefig(os.path.join(out_dir, os.path.basename(up_file_wild)+'-diff-{}.png'.format(title_key)), bbox_inches='tight', pad_inches=0.2
                   )
        fig.savefig(os.path.join(out_dir, os.path.basename(up_file_wild)+'-diff-{}.svg'.format(title_key)), bbox_inches='tight', pad_inches=0.2
                   )
        
        if plot_heatmap:
            fig = plt.figure(figsize=(2.5, 11),tight_layout={'w_pad':2})

        else:
            fig = plt.figure(figsize=(8, 1))
        labeldic = {0:True, 1:False}
        for iup, up_file in enumerate([up_file_wild, up_file_mut]):
            if plot_heatmap:
                ax = fig.add_subplot(131+iup*2) # Skip one ax in dirty way for cbar
                heatmap_up_dict(get_unpaired_probs(up_file), rang, title=os.path.basename(up_file).replace('-MUT-','-').replace('-WT-','-').replace('_lunp',''), 
                                fig=fig,ax=ax,ticklabel=labeldic[iup])            
            else:
                plot_up_dict(get_unpaired_probs(up_file), rang, title=os.path.basename(up_file).replace('-MUT-','-').replace('-WT-','-').replace('_lunp',''), 
                             fig=fig,tidy=True,
#                              ax=ax, ticklabel=labeldic[iup]
                            )            
            
#         fig.tight_layout(pad=0.1)

        fig.savefig(os.path.join(out_dir, os.path.basename(up_file_wild).replace('WT','WTMUT')+'-{}.png'.format(title_key)), bbox_inches='tight', #pad_inches=0.5,
                    dpi=600)
        fig.savefig(os.path.join(out_dir, os.path.basename(up_file_wild).replace('WT','WTMUT')+'-{}.svg'.format(title_key)), bbox_inches='tight', #pad_inches=0.5,
                    )

def plot_circos_seq_SNP(rec_wild, SNP_tag, local_fold=False, plotted_seq_lenght=None,dotplot=True,ECGplot=True,suffix='',annot_locs=[], annot_names=[],out_dir='./'):
    from Bio import SeqIO

    
    ID = '_'.join(rec_wild.id.split('|')[:2])+'_'+SNP_tag#"".join(x for x in rec_wild.id if x not in ['|','<', '>'])
    utr5_l, utr3_l = 0, 0 
    #     print "rec_wild: " , rec_wild
    wild_seq = rec_wild.seq
    # print wild_seq

    matches =  re.match('(\D)(\d+)(\D)', SNP_tag)
    if not matches:
        raise RuntimeError("No matches founs for tag:".format(SNP_tag)) 
    wild_char, loc, mut_char = matches.group(1), int(matches.group(2)), matches.group(3)
    assert(wild_seq[loc-1].upper() == wild_char.upper())
    mut_seq = wild_seq[:loc-1] + mut_char + wild_seq[loc:]
    
    #print mut_seq
    rec_mut = SeqRecord(mut_seq, id=ID+'-MUTANT')

    dp_wild, unp_wild = call_vienna_plfold(rec_wild.seq, ID, local_fold, out_dir=out_dir)
    create_circos_annotation(len(rec_wild), utr5_l, utr3_l, annot_locs, annot_names)
    run_dot2circ(dp_wild, ID+'-WILD'+suffix, out_dir=out_dir)
    

    dp_mut, unp_mut = call_vienna_plfold(rec_mut.seq, rec_mut.id, local_fold, out_dir=out_dir)
    create_circos_annotation(len(rec_mut), utr5_l, utr3_l, [int(SNP_tag[1:-1])]+annot_locs, [SNP_tag]+annot_names)
    run_dot2circ(dp_mut, rec_mut.id+suffix, out_dir=out_dir)
    dp_diff = dp_mut.replace('.ps', '_diff.ps')
    dpabs, dpremove, dpintroduce = write_diff_dp(dp_wild, dp_mut, dp_diff)

    run_dot2circ(dp_diff, rec_mut.id+suffix+'-diff', out_dir=out_dir)
    run_dot2circ(dpremove, rec_mut.id+suffix+'-removed', out_dir=out_dir)
    run_dot2circ(dpintroduce, rec_mut.id+suffix+'-introduced', out_dir=out_dir)
    
    if dotplot is True:
        ldp.plot_heat_maps(None, ldp.parse_dp_ps(dp_wild), filename=ID+'-WILD', what='basepairs',inverse=True, out_dir=out_dir)#, gene_loc=[2,10])
        ldp.plot_heat_maps(None, ldp.parse_dp_ps(dp_mut), filename=ID+'-MUTANT', what='basepairs',inverse=True, out_dir=out_dir)

        ldp.plot_heat_maps(None, ldp.parse_dp_ps(dp_diff), filename=ID+'-DIFF', what='basepairs',inverse=True, out_dir=out_dir)
        ldp.plot_heat_maps(None, ldp.parse_dp_ps(dpremove), filename=ID+'-REMOVED', what='basepairs',inverse=True, out_dir=out_dir)
        ldp.plot_heat_maps(None, ldp.parse_dp_ps(dpintroduce), filename=ID+'-INTRODUCED', what='basepairs',inverse=True, out_dir=out_dir)
        
        fig = ldp.plot_heat_maps(None, ldp.parse_dp_ps(dp_wild)+ldp.parse_dp_ps(dp_mut).transpose(), filename=ID+'-WT-MUT', what='basepairs',
                   inverse=True, interactive=False, gene_loc=None,title_suffix=ID+'-'+SNP_tag, out_dir=out_dir)
        fig.text(x=0.72,y=0.3,s='WT')
        fig.text(x=0.65,y=0.2,s=SNP_tag)
    
    if ECGplot is True:
#         plot_up_dict(u, None, title=ID, fig=myfig,tidy=True)
        plot_unpaired_probs([(unp_wild, unp_mut)], plot_heatmap=False, out_dir=out_dir)
        plot_unpaired_probs([(unp_wild, unp_mut)], plot_heatmap=True, out_dir=out_dir)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='MutaRNA-plot Predict and plot local and global base-pair probabilities if wildtype and mutant RNAs'\
        '\nSample call: \"python MutaRNA-plot.py --sequence ACGGGCACU --SNP-tag G3C'
        )
    parser.add_argument('--sequence-wild', required=True, type=is_valid_sequence, help='Input sequence string wildtype')
    parser.add_argument('--SNP-tag', required=True, type=str, help='SNP tag e.g. "C3G" for mutation at position 3 from C to G')
    parser.add_argument('--out-dir', default="./", type=is_valid_directory, help='output directory')


# Save to file in the current working directory

    args = parser.parse_args()
        
    plot_circos_seq_SNP(args.sequence_wild, args.SNP_tag, local_fold=True, out_dir=args.out_dir)

