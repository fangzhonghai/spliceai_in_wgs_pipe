# -*- coding:utf-8 -*-
import pandas as pd
import optparse
import pyfaidx
import yaml
import sys
import os


def print_usage(option, opt, value, parser):
    usage_message = """
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------

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


def bgi_anno_null_2_vcf(in_df):
    df = in_df.copy()
    df['#CHROM'] = []
    df['POS'] = []
    df['REF'] = []
    df['ALT'] = []
    df['ID'] = []
    df['QUAL'] = []
    df['FILTER'] = []
    df['INFO'] = []
    a = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']].copy()
    b = df[['#Chr', 'Start', 'Stop', 'Ref', 'Call', '#CHROM', 'POS', 'REF', 'ALT']].copy()
    return a, b, df


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


if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option('-u', '--usage', help='print more info on how to use this script', action="callback", callback=print_usage)
    parser.add_option('-b', '--bgi_anno', dest='bgi_anno', default=None, metavar='file')
    parser.add_option('-p', '--pwd', dest='pwd', default=None, metavar='string')
    parser.add_option('-o', '--out', dest='out', default=None, metavar='string')
    parser.add_option('-c', '--config', dest='config', default=None, metavar='file')
    parser.add_option('--in_format', dest='in_format', default='excel', metavar='string')
    parser.add_option('--sheet', dest='sheet', default='intron', metavar='string')
    parser.add_option('--out_format', dest='out_format', default='excel', type='string')
    (opts, args) = parser.parse_args()
    bgi_anno = opts.bgi_anno
    pwd = opts.pwd
    out = opts.out
    config = opts.config
    sheet = opts.sheet
    out_format = opts.out_format
    in_format = opts.in_format
    config_dic = yaml_read(config)
    vcf_header = read_vcf_head(config_dic['vcf_header'])
    with open(os.path.join(pwd, out + '.vcf'), 'w') as f:
        f.write(vcf_header)
    if in_format != 'excel':
        anno_df = pd.read_csv(bgi_anno, sep='\t')
    else:
        anno_df = pd.read_excel(bgi_anno, sheet_name=sheet)
    if anno_df.empty:
        vcf_df, bed, all_bed = bgi_anno_null_2_vcf(anno_df)
    else:
        vcf_df, bed, all_bed = bgi_anno_2_vcf_format(anno_df, config_dic['reference'])
    vcf_df.to_csv(os.path.join(pwd, out + '.vcf'), sep='\t', index=False, mode='a')
    bed.to_csv(os.path.join(pwd, out + '.bed2vcf.all.bed'), sep='\t', index=False)
    write_spliceai(config_dic, pwd, out)
