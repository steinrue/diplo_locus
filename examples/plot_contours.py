import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def load_surface(tablename):
    """return the matrix and its rolnames & colnames"""
    onGrid_LLs = pd.read_csv(tablename, delimiter="\t", comment="#") #, index_col=0, header=[0]
    # assume the DF has both row names and column names
    if "ID" in list(onGrid_LLs.columns):
        loci_names = list(onGrid_LLs.ID)
        # s_pairs = list(onGrid_LLs.columns[1:])
        s_pairs = [colname for colname in onGrid_LLs.columns if colname != 'ID']
    else:
        loci_names = list(onGrid_LLs.index)
        s_pairs = list(onGrid_LLs.columns)
    # make sure they're numbers
    # print(s_pairs)
    try:
        s_pairs = [list(map(float, pair.strip("()").split(','))) for pair in s_pairs]
    except:
        print(s_pairs[:5])
        print( s_pairs.index("ID"), len(s_pairs))
        print(s_pairs[-5:])
        sys.exit()
    s1_list, s2_list = map(lambda x: np.array(sorted(list(set(x)))),
                           zip(*s_pairs))
    LLcontainer = np.array(onGrid_LLs.loc[:, onGrid_LLs.columns != "ID"])
    # make sure s_pairs have tupples
    s_pairs = [tuple(pair) for pair in s_pairs]
    return np.array(LLcontainer), np.array(loci_names), s1_list, s2_list, s_pairs, len(s_pairs)


def main():
    # read commandline args
    inputfile, outprefix = sys.argv[1:]

    # SNPs = pd.read_csv(inputfile, delimiter="\t", header=0, index_col=0, comment="#")
    LLcontainer, loci_names, s1_list, s2_list, s_pairs, num_pairs = load_surface(inputfile)

    for i, snp in enumerate(loci_names):
        surface = np.array(LLcontainer[i,:].reshape( (len(s2_list)), len(s1_list)))
        # find neut
        neut_i, neut_j = list(s2_list).index(0), list(s1_list).index(0)
        snp_neut = surface[neut_i, neut_j] / np.log(10)
        max_idx = np.unravel_index(np.argmax(surface), surface.shape)
        snp_peak = surface[max_idx] / np.log(10)
        s1_0, s2_0 = s1_list[max_idx[1]], s2_list[max_idx[0]]
        print(snp, max_idx, s1_0, s2_0, snp_peak)
        # now let's plot
        figname = f'{outprefix}_{snp}_contour.png'
        fig, ax = plt.subplots(figsize=(4,4)) #
        ct = ax.contour(s1_list, s2_list, surface / np.log(10), levels=40, cmap='copper')
        ax.clabel(ct, ct.levels[1::2], inline=True)
        # ax.plot ([0,0], [-0.5,0.5], "--", c='gray', lw=0.5)
        # ax.plot ([-0.5,0.5], [-0,0], "--", c='gray', lw=0.5)
        # ax.set_axis_on()
        ax.set_xlabel('$s_{Aa}$')
        ax.set_ylabel('$s_{AA}$')
        mlr = (snp_peak-snp_neut) * 2
        if '_' in snp:
            ch, pos, rsid = snp.split("_")
            ax.set_title(f'SNP {rsid} on chr{ch}, {pos}\nLR={mlr : .3f}', size='small')
        else:
            ax.set_title(f'Locus {snp}\nLR={mlr : .3f}', size='small')
        # ax.set_title(f'%s on chr%s, %s\nLR=%.4f', size='small')
        ax.plot(s1_0, s2_0, 'x')
        ax.annotate('$log_{10}\mathcal{L}_{max}=$\n%.4f\n(%.3f, %.3f)' % (snp_peak, s1_0, s2_0),
                             (s1_0, s2_0), xytext=(s1_0, s2_0), ha='left', va='top', fontsize='small',
                             textcoords='offset points', annotation_clip=True)
        # add axis line
        ax.axvline(x=0, ymin=0, ymax=1, ls='--', color='#888888', lw=0.5)
        ax.axhline(y=0, xmin=0, xmax=1, ls='--', color='#888888', lw=0.5)
        # ax.set_xlim(-0.1, 0.102)
        fig.tight_layout()
        plt.savefig(figname, dpi=350)


if __name__ == '__main__':
    main()
