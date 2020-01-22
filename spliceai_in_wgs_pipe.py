# -*- coding:utf-8 -*-
import pandas as pd
import optparse
import pyfaidx
import yaml
import sys
import os


def print_usage(option, opt, value, parser):
    usage_message = r"""
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
    python spliceai_in_wgs_pipe.py -b test.xlsx --sheet intron --out_format excel \
    -o prefix -p /path/to/work -c spliceai.yaml
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
    """
    print(usage_message)
    sys.exit()


def read_vcf_head(header):
    vh = open(header, 'r')
    fp_read = vh.read()
    vh.close()
    return fp_read


def yaml_read(yaml_file):
    with open(yaml_file, 'r') as y:
        yaml_dic = yaml.load(y, Loader=yaml.FullLoader)
    return yaml_dic


def bgi_anno_2_vcf_format(in_df, reference):
    df = in_df.copy()
    df['#Chr'] = df['#Chr'].astype('str')
    if len(df[df['#Chr'].str.startswith('chr')]):
        df['#CHROM'] = df['#Chr']
    else:
        df['#CHROM'] = 'chr' + df['#Chr']
    df.loc[df['#CHROM'] == 'chrMT', '#CHROM'] = 'chrM_NC_012920.1'
    df['ID'] = '.'
    df['QUAL'] = '.'
    df['FILTER'] = '.'
    df['INFO'] = '.'
    df['MuType'] = 'delins'
    df.loc[df['Ref'] == '.', 'MuType'] = 'ins'
    df.loc[df['Call'] == '.', 'MuType'] = 'del'
    df.loc[(df['Ref'].map(len) == 1) & (df['Call'].map(len) == 1) & (df['Ref'] != '.') & (df['Call'] != '.'), 'MuType'] = 'snp'
    df['POS'] = df['Stop']
    df.loc[df['MuType'] == 'del', 'POS'] = df.loc[df['MuType'] == 'del', 'Start']
    df.loc[df['MuType'] == 'delins', 'POS'] = df.loc[df['MuType'] == 'delins', 'Start'] + 1
    df['REF'] = df['Ref']
    df['ALT'] = df['Call']
    fa = pyfaidx.Fasta(reference)
    for i in range(df.shape[0]):
        if df.loc[i, 'MuType'] == 'ins':
            base = str(fa.get_seq(df.loc[i, '#CHROM'], df.loc[i, 'POS'], df.loc[i, 'POS'])).upper()
            df.loc[i, 'REF'] = base
            df.loc[i, 'ALT'] = base + df.loc[i, 'ALT']
        elif df.loc[i, 'MuType'] == 'del':
            base = str(fa.get_seq(df.loc[i, '#CHROM'], df.loc[i, 'POS'], df.loc[i, 'POS'])).upper()
            df.loc[i, 'ALT'] = base
            df.loc[i, 'REF'] = base + df.loc[i, 'REF']
        else:
            pass
    a = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']].copy()
    a.sort_values(by=['#CHROM', 'POS'], ascending=True, inplace=True)
    b = df[['#Chr', 'Start', 'Stop', 'Ref', 'Call', '#CHROM', 'POS', 'REF', 'ALT']].copy()
    df.drop(columns=['ID', 'QUAL', 'FILTER', 'INFO', 'MuType'], inplace=True)
    return a, b, df


def write_spliceai(yaml_dic, path, prefix):
    reference = yaml_dic['reference']
    spliceai = yaml_dic['spliceai']
    process = yaml_dic['process']
    with open(os.path.join(path, 'spliceai.' + prefix + '.sh'), 'w') as sp:
        shell = r'''#!/bin/bash
{spliceai} -I {path}/{prefix}.vcf -O {path}/{prefix}.spliceai.vcf -A grch37 -R {reference} -D 500 -P {process}
'''.format(**locals())
        sp.write(shell)


def spliceai_max_score(spliceai_res):
    spliceai_res_na = spliceai_res[(spliceai_res['SpliceAI'] == '.') | (spliceai_res['SpliceAI'].isna())].copy()
    spliceai_res_is = spliceai_res[(spliceai_res['SpliceAI'] != '.') & (~spliceai_res['SpliceAI'].isna())].copy()
    spliceai_res_na['SpliceAI Pred'] = '.'
    spliceai_res_is['SpliceAI Pred'] = spliceai_res_is['SpliceAI'].str.extract('SpliceAI=(.*?)$')
    spliceai_res_is['SpliceAI Pred'] = spliceai_res_is['SpliceAI Pred'].str.split(',').str[0]
    spliceai_res_is['SpliceAI Pred'] = spliceai_res_is['SpliceAI Pred'].str.split('|').str[2:6]
    spliceai_res_is['SpliceAI Pred'] = spliceai_res_is.apply(lambda x: max(x['SpliceAI Pred']), axis=1)
    spliceai_res_final = spliceai_res_is.append(spliceai_res_na)
    return spliceai_res_final


def merge_spliceai_res(path, prefix, skip_rows, input_bed):
    cols = ['#CHROM', 'POS', 'REF', 'ALT', 'SpliceAI']
    spliceai_df = pd.read_csv(path + '/' + prefix + '.spliceai.vcf', sep='\t', skiprows=range(skip_rows))
    spliceai_df.rename(columns={'INFO': 'SpliceAI'}, inplace=True)
    spliceai_res = spliceai_df[cols].copy()
    merged_spliceai = pd.merge(input_bed, spliceai_res, on=cols[:-1], how='left')
    merged_spliceai.drop_duplicates(inplace=True)
    merged_spliceai_final = spliceai_max_score(merged_spliceai)
    merged_spliceai_final.drop(columns=['#CHROM', 'POS', 'REF', 'ALT'], inplace=True)
    return merged_spliceai_final


if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('-u', '--usage', help='print more info on how to use this script', action="callback", callback=print_usage)
    parser.add_option('-b', '--bgi_anno', dest='bgi_anno', default=None, metavar='file')
    parser.add_option('-p', '--pwd', dest='pwd', default=None, metavar='string')
    parser.add_option('-o', '--out', dest='out', default=None, metavar='string')
    parser.add_option('-c', '--config', dest='config', default=None, metavar='file')
    parser.add_option('--sheet', dest='sheet', default='intron', metavar='string')
    parser.add_option('--out_format', dest='out_format', default='excel', type='string')
    (opts, args) = parser.parse_args()
    bgi_anno = opts.bgi_anno
    pwd = opts.pwd
    out = opts.out
    config = opts.config
    sheet = opts.sheet
    out_format = opts.out_format
    config_dic = yaml_read(config)
    anno_df = pd.read_excel(bgi_anno, sheet_name=sheet)
    vcf_df, bed, all_bed = bgi_anno_2_vcf_format(anno_df, config_dic['reference'])
    vcf_header = read_vcf_head(config_dic['vcf_header'])
    with open(os.path.join(pwd, out + '.vcf'), 'w') as f:
        f.write(vcf_header)
    vcf_df.to_csv(os.path.join(pwd, out + '.vcf'), sep='\t', index=False, mode='a')
    write_spliceai(config_dic, pwd, out)
    status = os.system('sh ' + os.path.join(pwd, 'spliceai.' + out + '.sh'))
    if status != 0:
        sys.exit(1)
    splice_pred = merge_spliceai_res(pwd, out, config_dic['skip_vcf_header_row'], all_bed)
    if out_format == 'excel':
        splice_pred.to_excel(pwd + '/' + out + '.spliceai.xlsx', index=False)
    else:
        splice_pred.to_csv(pwd + '/' + out + '.spliceai.tsv', index=False, sep='\t')
