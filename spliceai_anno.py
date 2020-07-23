# -*- coding:utf-8 -*-
from multiprocessing import Pool, cpu_count
from functools import partial, reduce
from io import StringIO
import pandas as pd
import subprocess
import optparse
import yaml
import sys
import os


def print_usage(option, opt, value, parser):
    usage_message = """
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
  python3 spliceai_anno.py -i 20B2752186.out -o 20B2752186.out --process 4
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
    """
    print(usage_message)
    sys.exit()


def yaml_read(yaml_file):
    with open(yaml_file, 'r') as y:
        yaml_dic = yaml.load(y, Loader=yaml.FullLoader)
    return yaml_dic


def run_bedtools(bedtools, bed):
    command = bedtools + ' sort -i ' + bed + '|' + bedtools + ' merge -i -'
    output = subprocess.getoutput(command)
    return StringIO(output)


def run_tabix(tabix, vcf, bed):
    command = tabix + ' -h -R ' + bed + ' ' + vcf
    output = subprocess.getoutput(command)
    splice = pd.read_csv(StringIO(output), skiprows=range(28), sep='\t', dtype={'#CHROM': str})
    splice.rename(columns={'INFO': 'SpliceAI'}, inplace=True)
    splice.drop(columns=['ID', 'QUAL', 'FILTER'], inplace=True)
    splice_anno_format = vcf_format_2_bgi_anno(splice)
    return splice_anno_format


def spliceai_max_score(spliceai_res):
    spliceai_res_na = spliceai_res[(spliceai_res['SpliceAI'] == '.') | (spliceai_res['SpliceAI'].isna())].copy()
    spliceai_res_is = spliceai_res[(spliceai_res['SpliceAI'] != '.') & (~spliceai_res['SpliceAI'].isna())].copy()
    spliceai_res_na['SpliceAI Pred'] = '.'
    if not spliceai_res_is.empty:
        spliceai_res_is['SpliceAI Pred'] = spliceai_res_is['SpliceAI'].str.extract('SpliceAI=(.*?)$')
        spliceai_res_is['SpliceAI Pred'] = spliceai_res_is['SpliceAI Pred'].str.split(',').str[0]
        spliceai_res_is['SpliceAI Pred'] = spliceai_res_is['SpliceAI Pred'].str.split('|').str[2:6]
        spliceai_res_is['SpliceAI Pred'] = spliceai_res_is.apply(lambda x: max(x['SpliceAI Pred']), axis=1)
    spliceai_res_final = spliceai_res_is.append(spliceai_res_na, sort=False)
    return spliceai_res_final


def interpretation(spliceai_list, chrom, pos):
    threshold = 0.1
    ds = spliceai_list[0:4]
    dp = spliceai_list[4:8]
    ds = [float(i) for i in ds]
    dp = [int(i) for i in dp]
    ds_dic = {0: 'acceptor gain', 1: 'acceptor loss', 2: 'donor gain', 3: 'donor loss'}
    ds_index = [i for i, j in enumerate(ds) if j >= threshold]
    if len(ds_index) > 0:
        pred = []
        for k in ds_index:
            if dp[k] < 0:
                comment = str(chrom) + ':' + str(pos + dp[k]) + ' (=' + str(pos) + str(dp[k]) + ') ' + ds_dic[k] + ' ' + str(ds[k])
            else:
                comment = str(chrom) + ':' + str(pos + dp[k]) + ' (=' + str(pos) + '+' + str(dp[k]) + ') ' + ds_dic[k] + ' ' + str(ds[k])
            pred.append(comment)
        pred_content = ';'.join(pred)
        return pred_content
    else:
        return '.'


def spliceai_interpre(spliceai_res):
    spliceai_res_na = spliceai_res[(spliceai_res['SpliceAI Pred'] == '.') | (spliceai_res['SpliceAI Pred'].isna())].copy()
    spliceai_res_is = spliceai_res[(spliceai_res['SpliceAI Pred'] != '.') & (~spliceai_res['SpliceAI Pred'].isna())].copy()
    spliceai_res_na['SpliceAI Interpretation'] = '.'
    if not spliceai_res_is.empty:
        spliceai_res_is['SpliceAI Interpretation'] = spliceai_res_is['SpliceAI'].str.extract('SpliceAI=(.*?)$')
        spliceai_res_is['SpliceAI Interpretation'] = spliceai_res_is['SpliceAI Interpretation'].str.split(',').str[0]
        spliceai_res_is['SpliceAI Interpretation'] = spliceai_res_is['SpliceAI Interpretation'].str.split('|').str[2:10]
        spliceai_res_is['SpliceAI Interpretation'] = spliceai_res_is.apply(lambda x: interpretation(x['SpliceAI Interpretation'],
                                                                                                    x['#CHROM'], x['POS']), axis=1)
    spliceai_res_final = spliceai_res_is.append(spliceai_res_na, sort=False)
    return spliceai_res_final


def judge(pred):
    if pred != '.' and float(pred) >= 0.2:
        return 'D'
    elif pred != '.' and float(pred) < 0.2:
        return 'P'
    else:
        return '.'


def change_res(in_df):
    in_df['SpliceAI Pred'] = in_df.apply(lambda x: judge(x['SpliceAI Pred']), axis=1)
    return in_df


def vcf_format_2_bgi_anno(in_df):
    df = in_df.copy()
    df['MuType'] = '.'
    df.loc[(df['REF'].str.len() == 1) & (df['ALT'].str.len() == 1), 'MuType'] = 'snv'
    df.loc[(df['REF'].str.len() == 1) & (df['ALT'].str.len() > 1), 'MuType'] = 'ins'
    df.loc[(df['REF'].str.len() > 1) & (df['ALT'].str.len() == 1), 'MuType'] = 'del'
    df.loc[(df['REF'].str.len() > 1) & (df['ALT'].str.len() > 1) & (df['REF'].str[0] == df['ALT'].str[0]), 'MuType'] = 'delins_eq'
    df.loc[(df['REF'].str.len() > 1) & (df['ALT'].str.len() > 1) & (df['REF'].str[0] != df['ALT'].str[0]), 'MuType'] = 'delins_neq'
    df['#Chr'] = df['#CHROM']
    # df['#CHROM'] = df['#CHROM'].astype('str')
    # if len(df[df['#CHROM'].str.startswith('chr')]):
    #     df['#Chr'] = df['#CHROM']
    # else:
    #     df['#Chr'] = 'chr' + df['#CHROM']
    # df.loc[df['#CHROM'] == 'chrM_NC_012920.1', '#Chr'] = 'chrMT'
    df['Start'] = df['POS'].astype(int) - 1
    df['Stop'] = df['POS'].astype(int)
    df['Ref'] = df['REF']
    df['Call'] = df['ALT']
    df.loc[df['MuType'] == 'ins', 'Ref'] = '.'
    df.loc[df['MuType'] == 'ins', 'Call'] = df.loc[df['MuType'] == 'ins', 'ALT'].str[1:]
    df.loc[df['MuType'] == 'del', 'Call'] = '.'
    df.loc[df['MuType'] == 'del', 'Ref'] = df.loc[df['MuType'] == 'del', 'REF'].str[1:]
    df.loc[df['MuType'] == 'ins', 'Start'] = df.loc[df['MuType'] == 'ins', 'Stop']
    df.loc[df['MuType'] == 'del', 'Start'] = df.loc[df['MuType'] == 'del', 'Stop']
    df.loc[df['MuType'] == 'del', 'Stop'] = df.loc[df['MuType'] == 'del', 'Stop'] + df.loc[df['MuType'] == 'del', 'Ref'].str.len()
    df.loc[df['MuType'] == 'delins_eq', 'Stop'] = df.loc[df['MuType'] == 'delins_eq', 'Start'] + df.loc[df['MuType'] == 'delins_eq', 'Ref'].str.len()
    df.loc[df['MuType'] == 'delins_eq', 'Start'] = df.loc[df['MuType'] == 'delins_eq', 'POS']
    df.loc[df['MuType'] == 'delins_eq', 'Ref'] = df.loc[df['MuType'] == 'delins_eq', 'REF'].str[1:]
    df.loc[df['MuType'] == 'delins_eq', 'Call'] = df.loc[df['MuType'] == 'delins_eq', 'ALT'].str[1:]
    df.loc[df['MuType'] == 'delins_neq', 'Stop'] = df.loc[df['MuType'] == 'delins_neq', 'Start'] + df.loc[df['MuType'] == 'delins_neq', 'Ref'].str.len()
    # a = df[['#CHROM', 'POS', 'REF', 'ALT', '#Chr', 'Start', 'Stop', 'Ref', 'Call']].copy()
    df.drop(columns=['MuType'], inplace=True)
    return df


def split_snv_indel(in_df):
    df = in_df.copy()
    df['MuType'] = 'indel'
    df.loc[(df['Ref'].map(len) == 1) & (df['Call'].map(len) == 1) & (df['Ref'] != '.') & (df['Call'] != '.'), 'MuType'] = 'snv'
    df_snv = df[df['MuType'] == 'snv']
    df_indel = df[df['MuType'] == 'indel']
    return df_snv[['#Chr', 'Start', 'Stop']], df_indel[['#Chr', 'Start', 'Stop']]


def split_df(df, split_num):
    df.reset_index(drop=True, inplace=True)
    df_list = list()
    step = round(df.shape[0]/split_num)
    for i in range(split_num):
        if i == 0:
            df_list.append(df.loc[0: step-1])
        elif i == split_num-1:
            df_list.append(df.loc[step*i:])
        else:
            df_list.append(df.loc[step*i:step*(i+1)-1])
    return df_list


Path = os.path.split(os.path.realpath(__file__))[0]


def main():
    parser = optparse.OptionParser()
    parser.add_option('-u', '--usage', help='print more info on how to use this script', action="callback", callback=print_usage)
    parser.add_option('-b', '--bgi_anno', dest='bgi_anno', default=None, metavar='file')
    parser.add_option('-p', '--pwd', dest='pwd', default=Path, metavar='string')
    parser.add_option('-o', '--out', dest='out', default=None, metavar='string')
    parser.add_option('--skip_rows', dest='skip_rows', default=0, type=int)
    parser.add_option('--process', dest='process', default=cpu_count(), type=int)
    parser.add_option('-c', '--config', dest='config', default=os.path.join(Path, 'etc', 'spliceai.yaml'), metavar='file')
    (opts, args) = parser.parse_args()
    bgi_anno = opts.bgi_anno
    pwd = opts.pwd
    out = opts.out
    skip_rows = opts.skip_rows
    config = opts.config
    process_num = min(opts.process, cpu_count())
    if not os.path.exists(config):
        print(config + ' is not exist')
        sys.exit(1)
    config_dic = yaml_read(config)
    spliceai_snv = config_dic['spliceai_snv']
    if not os.path.exists(spliceai_snv):
        spliceai_snv = os.path.join(pwd, config_dic['spliceai_snv'])
        if not os.path.exists(spliceai_snv):
            print(spliceai_snv + ' is not exist')
            sys.exit(1)
    spliceai_indel = config_dic['spliceai_indel']
    if not os.path.exists(spliceai_indel):
        spliceai_indel = os.path.join(pwd, config_dic['spliceai_indel'])
        if not os.path.exists(spliceai_indel):
            print(spliceai_indel + ' is not exist')
            sys.exit(1)
    wes_spliceai_gene = config_dic['wes_spliceai_gene']
    if not os.path.exists(wes_spliceai_gene):
        wes_spliceai_gene = os.path.join(pwd, config_dic['wes_spliceai_gene'])
        if not os.path.exists(wes_spliceai_gene):
            print(wes_spliceai_gene + ' is not exist')
            sys.exit(1)
    wes_spliceai_gene_df = pd.read_csv(wes_spliceai_gene, sep='\t')
    bedtools = config_dic['bedtools']
    if os.path.exists(os.path.join(pwd, bedtools)):
        bedtools = os.path.join(pwd, bedtools)
    tabix = config_dic['tabix']
    if os.path.exists(os.path.join(pwd, tabix)):
        tabix = os.path.join(pwd, tabix)
    anno_df = pd.read_csv(bgi_anno, sep='\t', low_memory=False, dtype={'#Chr': str}, skiprows=range(skip_rows))
    if anno_df.empty:
        anno_df['SpliceAI'] = []
        anno_df['SpliceAI Pred'] = []
        anno_df['SpliceAI Interpretation'] = []
        anno_df.to_csv(out + '.spliceai.tsv', sep='\t', index=False)
        sys.exit(0)
    if len(anno_df[anno_df['#Chr'].str.startswith('chr')]):
        anno_df.replace({'#Chr': r'^chr'}, {'#Chr': ''}, inplace=True, regex=True)
    df_bed_snv, df_bed_indel = split_snv_indel(anno_df)
    df_bed_snv.to_csv(out + '.bgianno.snv.bed', sep='\t', index=False)
    df_bed_indel.to_csv(out + '.bgianno.indel.bed', sep='\t', index=False)
    splice_snv_all = pd.DataFrame(columns=['#CHROM', 'POS', 'REF', 'ALT', 'SpliceAI', '#Chr', 'Start', 'Stop', 'Ref', 'Call'])
    splice_indel_all = pd.DataFrame(columns=['#CHROM', 'POS', 'REF', 'ALT', 'SpliceAI', '#Chr', 'Start', 'Stop', 'Ref', 'Call'])
    if not df_bed_snv.empty:
        bed_snv = pd.read_csv(run_bedtools(bedtools, out + '.bgianno.snv.bed'), sep='\t', header=None)
        bed_snv_df_list = split_df(bed_snv, min(process_num, bed_snv.shape[0]))
        bed_snv_df_list_files = list()
        for i in range(len(bed_snv_df_list)):
            bed_snv_df_list[i].to_csv(out + '.bgianno.snv.sort_merge.bed.' + str(i), sep='\t', index=False, header=None)
            bed_snv_df_list_files.append(out + '.bgianno.snv.sort_merge.bed.' + str(i))
        partial_run_snv = partial(run_tabix, tabix, spliceai_snv)
        with Pool(process_num) as pool:
            splice_snv_list = pool.map(partial_run_snv, bed_snv_df_list_files)
        splice_snv_all = reduce(lambda x, y: x.append(y), splice_snv_list)
    if not df_bed_indel.empty:
        bed_indel = pd.read_csv(run_bedtools(bedtools, out + '.bgianno.indel.bed'), sep='\t', header=None)
        bed_indel_df_list = split_df(bed_indel, min(process_num, bed_indel.shape[0]))
        bed_indel_df_list_files = list()
        for i in range(len(bed_indel_df_list)):
            bed_indel_df_list[i].to_csv(out + '.bgianno.indel.sort_merge.bed.' + str(i), sep='\t', index=False, header=None)
            bed_indel_df_list_files.append(out + '.bgianno.indel.sort_merge.bed.' + str(i))
        partial_run_indel = partial(run_tabix, tabix, spliceai_indel)
        with Pool(process_num) as pool:
            splice_indel_list = pool.map(partial_run_indel, bed_indel_df_list_files)
        splice_indel_all = reduce(lambda x, y: x.append(y), splice_indel_list)
    splice_res = splice_snv_all.append(splice_indel_all, sort=False)
    splice_res['SpliceAI Gene'] = splice_res['SpliceAI'].str.extract(r'\|(.*?)\|')
    splice_res_gene = pd.merge(splice_res, wes_spliceai_gene_df, on=['SpliceAI Gene'])
    splice_res_gene.drop(columns=['SpliceAI Gene'], inplace=True)
    splice_res_max = spliceai_max_score(splice_res_gene)
    splice_res_max_inter = spliceai_interpre(splice_res_max)
    splice_res_max_inter_change = change_res(splice_res_max_inter)
    splice_res_max_inter_change.drop(columns=['#CHROM', 'POS', 'REF', 'ALT'], inplace=True)
    splice_res_merge = pd.merge(anno_df, splice_res_max_inter_change, on=['#Chr', 'Start', 'Stop', 'Ref', 'Call', 'Gene Symbol'], how='left')
    splice_res_merge.drop_duplicates(inplace=True)
    splice_res_merge.fillna(value={'SpliceAI': '.', 'SpliceAI Pred': '.', 'SpliceAI Interpretation': '.'}, inplace=True)
    splice_res_merge.to_csv(out + '.spliceai.tsv', sep='\t', index=False)


if __name__ == '__main__':
    main()
