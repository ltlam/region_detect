#==============================================================================
# bind_dmr_col.py
# Created By: Larry Lam
# Date Created: 3/08/2016
#==============================================================================
import argparse
from scipy import stats
import numpy as np
import pandas as pd
import re

STR_SEP = ','
STR_FRAG_LBL = 'frag'
STR_CPGNUM_LBL = 'cpg_num'

def get_float(str_val):
    if str_val == 'NA':
        return float('nan')
    return float(str_val)

def get_fragrow(d_col_to_val,l_sam_ids):
    str_row = ''
    n_not_nan_sum = 0
    for i in sorted(d_col_to_val.keys()):
        col_nan = np.isnan(d_col_to_val[i])
        if i in l_sam_ids:
            col_not_nan = [not isnan for isnan in col_nan]
            n_not_nan_sum += sum(col_not_nan)
        if sum(col_nan) == len(d_col_to_val[i]):
            str_row += 'NA%s'%STR_SEP
        else:
            str_row += '%f%s'%(np.nanmean(d_col_to_val[i]),STR_SEP)
    return re.sub(',$','',str_row)

def collapse_data_frame(str_frag_matrix,str_out_file,l_metrics,l_sam_ids):
    fout = open(str_out_file,'w')
    str_cur_frag = ''
    b_header = True

    d_col_to_vals = {}
    for i in l_metrics:
        d_col_to_vals[i] = []

    n_frag_id = 0
    for str_line in open(str_frag_matrix,'r'):
        lcols = str_line.rstrip().split(',')

        #write column header
        if b_header:
            n_frag_id = lcols.index(STR_FRAG_LBL)
            str_row = lcols[n_frag_id]
            for i in l_metrics:
                str_row+='%s%s'%(STR_SEP,lcols[i])
            fout.write(str_row+'\n')
            b_header = False
            continue

        #reset new fragment
        if lcols[n_frag_id] != str_cur_frag and str_cur_frag != '':
            fout.write('%s%s%s\n'%(str_cur_frag,STR_SEP,get_fragrow(d_col_to_vals,l_sam_ids)))
            str_cur_frag = lcols[n_frag_id]
            for i in l_metrics:
                d_col_to_vals[i] = []

        for i in l_metrics:
            d_col_to_vals[i].append(get_float(lcols[i]))
        str_cur_frag = lcols[n_frag_id]

    if len(d_col_to_vals[l_metrics[0]]) > 0:
        fout.write('%s%s%s\n'%(str_cur_frag,STR_SEP,get_fragrow(d_col_to_vals,l_sam_ids)))
    fout.close()

#-----------------------------------------------------------------------------------------------------------------------
# extend_reset_frag_info
# [in/out] l_frag_lbls - list of fragments to bind to the data frame
# [in/out] l_n_cpgs - list of cpgs within the fragment to bind to the data frame
# [in] l_prev_pos - list of positions within the previous fragment
# [in] str_prev_chr - previous chromosome
#-----------------------------------------------------------------------------------------------------------------------
def extend_reset_frag_info(l_frag_lbls,l_n_cpgs,l_prev_pos,str_prev_chr):
    if len(l_prev_pos) > 0:
        l_pre_frags = ['%s:%d-%d'%(str_prev_chr,l_prev_pos[0],l_prev_pos[-1])] * len(l_prev_pos)
        l_frag_lbls.extend(l_pre_frags)
        l_n_cpgs.extend([len(l_prev_pos)]*len(l_prev_pos))

#-----------------------------------------------------------------------------------------------------------------------
# get_cpg_similarity - return pearson r sq, remove nan values
# [in] l_xval - list of numeric values
# [in] l_yval - list of numeric values
# [out] r_sq,p_val,max_delta - pearson rsq, and pearson p-value
#-----------------------------------------------------------------------------------------------------------------------
def get_cpg_similarity(l_xval,l_yval,b_ignorenan):
    if (np.isnan(l_xval).any() or np.isnan(l_yval).any()) and b_ignorenan:
        l_x_nan_id = [i for i,x in enumerate(l_xval) if np.isnan(x)]
        l_y_nan_id = [i for i,x in enumerate(l_yval) if np.isnan(x)]
        l_nan_ids = []
        for xid in l_x_nan_id:
            l_nan_ids.append(xid)
        for yid in l_y_nan_id:
            l_nan_ids.append(yid)
        nan_id = set(l_nan_ids)
        if len(nan_id) == len(l_xval):
            return float('nan'),float('nan'),float('nan')
        l_x_trim = [x for i,x in enumerate(l_xval) if not (i in nan_id)]
        l_y_trim = [x for i,x in enumerate(l_yval) if not (i in nan_id)]
        max_delta = max([abs(a - b) for a, b in zip(l_x_trim, l_y_trim)])
        r_sq,r_pval = stats.pearsonr(l_x_trim, l_y_trim)
    else:
        r_sq,r_pval = stats.pearsonr(l_xval, l_yval)
        max_delta = max([abs(a - b) for a, b in zip(l_xval, l_yval)])
    return r_sq,r_pval,max_delta

#-----------------------------------------------------------------------------------------------------------------------
# frag_bind
# [in] sort_df - sorted data frame
# [in] n_bp - bp to neighboring cpg
# [in] f_pearson_r_sq - correlation filter
# [in] f_pearson_pval - correlation pvfilter
# [in] n_rep_num - number of replicates
# [out] frag_df - sorted data frame with column fragment
#-----------------------------------------------------------------------------------------------------------------------
def frag_bind(sort_df,n_bp,f_pearson_r_sq,f_pearson_pval,l_ref_id,f_max_delta_thresh):
    #new col variables
    l_frag_count = []
    l_frags = []

    #temp frag var
    l_temp_frag_pos = []
    str_temp_prev_chr = ''
    l_temp_center_meth = []
    for lrow in sort_df.itertuples():
        #init row variables
        str_cur_chr = lrow[-2]
        n_cur_pos = lrow[-1]
        str_cpg_site = lrow[0]

        l_cur_rep_vals = [float(meth) for id, meth in enumerate(lrow) if id in l_ref_id]

        if len(l_temp_center_meth) == 0:
            l_temp_center_meth = l_cur_rep_vals

        #
        #if np.isnan(l_cur_rep_vals).any():
        #    print 'Has NAN'
        r_sq,r_pval,max_delta = get_cpg_similarity(l_cur_rep_vals, l_temp_center_meth,True)

        #extend reset for new file or new chr
        if (np.isnan(max_delta) or max_delta >= f_max_delta_thresh or len(l_temp_frag_pos) == 0 or str_temp_prev_chr != str_cur_chr) or (len(l_temp_frag_pos) > 0 and (n_cur_pos - l_temp_frag_pos[-1]) > n_bp) or (len(l_temp_frag_pos) > 0 and (r_sq < f_pearson_r_sq or r_pval > f_pearson_pval)):
            extend_reset_frag_info(l_frags,l_frag_count,l_temp_frag_pos,str_temp_prev_chr)
            l_temp_center_meth = l_cur_rep_vals
            str_temp_prev_chr = str_cur_chr
            l_temp_frag_pos = []

        #append position to frag temp
        l_temp_frag_pos.append(n_cur_pos)
        #l_temp_center_meth = [sum(x)/2.0 for x in zip(l_cur_rep_vals, l_temp_center_meth)]
        l_temp_center_meth = l_cur_rep_vals


    #extend if remaining temp pos
    extend_reset_frag_info(l_frags,l_frag_count,l_temp_frag_pos,str_temp_prev_chr)
    sort_df[STR_FRAG_LBL] = pd.Series(l_frags,index=sort_df.index)
    sort_df[STR_CPGNUM_LBL] = pd.Series(l_frag_count,index=sort_df.index)
    return sort_df

#-----------------------------------------------------------------------------------------------------------------------
# order_df
# [in] str_csv - path to pre processing csv
# [out] pre_meth_df - sorted data frame by primary and secondary col id chr and pos
#-----------------------------------------------------------------------------------------------------------------------
def order_df(str_csv):
    if str_csv.endswith('.csv'):
        pre_meth_df = pd.read_csv(str_csv,index_col=0)
    else:
        pre_meth_df = pd.read_csv(str_csv,sep='\t',index_col=0)
    l_str_chr = [ str_cpg.split('_')[0] for str_cpg in list(pre_meth_df.index)]
    l_n_pos = [ int(str_cpg.split('_')[1]) for str_cpg in list(pre_meth_df.index)]
    pre_meth_df['chr'] = pd.Series(l_str_chr,index=pre_meth_df.index)
    pre_meth_df['pos'] = pd.Series(l_n_pos,index=pre_meth_df.index)
    return pre_meth_df.sort_values(by=['chr','pos'],ascending=[1,1]).drop_duplicates()

def get_parser(parser=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-dmr_join','--dmr_join',type=str,help="path to sorted reference and sample methylation")
    parser.add_argument('-dmr_bpdist','--dmr_bpdist',type=int,help="distance to neighboring site",default=500)
    parser.add_argument('-dmr_rsq','--dmr_rsq',type=float,help="correlation minimum",default=-1.00)
    parser.add_argument('-dmr_pval','--dmr_pval',type=float,help="correlation pval max",default=1.00)
    parser.add_argument('-dmr_cpgmin','--dmr_cpgmin',type=int,help="minimum number of cpg for collapse",default=2)
    parser.add_argument('-dmr_maxdelta','--dmr_maxdelta',type=float,help="maximum delta between neighboring sites",default=0.15)

    return parser

def main(str_ref_sam_join,n_dist,f_rsq,f_pval,f_max_delta):
    #sort data frame by 2 col
    meth_df_sorted = order_df(str_ref_sam_join)
    print(list(meth_df_sorted.columns.values))

    #get col id of labels
    l_ref_id = range(1, len(meth_df_sorted.columns)-1)

    #bind fragment column
    frag_df = frag_bind(meth_df_sorted,n_dist,f_rsq,f_pval,l_ref_id,f_max_delta)
    str_frag_csv = re.sub('.csv$','_frag.csv',str_ref_sam_join)
    frag_df.to_csv(str_frag_csv,na_rep='NA')

    #collapse
    l_val_info = l_ref_id
    #l_val_info.append(len(frag_df.columns.values))
    l_val_info.append(len(frag_df.columns.values))
    str_collapse_csv = re.sub('.csv$','_collapse.csv',str_frag_csv)
    collapse_data_frame(str_frag_csv,str_collapse_csv,l_val_info,l_ref_id)

#-----------------------------------------------------------------------------------------------------------------------
# main
#-----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    main(args.dmr_join,args.dmr_bpdist,args.dmr_rsq,args.dmr_pval,args.dmr_maxdelta)