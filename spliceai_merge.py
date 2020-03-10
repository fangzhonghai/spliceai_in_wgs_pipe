# -*- coding:utf-8 -*-
import pandas as pd
import optparse
import yaml
import sys


def print_usage(option, opt, value, parser):
    usage_message = """
# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------
    """
    print(usage_message)
    sys.exit()


def yaml_read(yaml_file):
    with open(yaml_file, 'r') as y:
        yaml_dic = yaml.load(y, Loader=yaml.FullLoader)
    return yaml_dic


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


def merge_spliceai_res(path, prefix, skip_rows, input_bed, anno):
    cols = ['#CHROM', 'POS', 'REF', 'ALT', 'SpliceAI']
    spliceai_df = pd.read_csv(path + '/' + prefix + '.spliceai.vcf', sep='\t', skiprows=range(skip_rows))
    spliceai_df.rename(columns={'INFO': 'SpliceAI'}, inplace=True)
    spliceai_res = spliceai_df[cols].copy()
    spliceai_res.drop_duplicates(inplace=True)
    spliceai_res_pred = spliceai_max_score(spliceai_res)
    spliceai_res_pred_interpre = spliceai_interpre(spliceai_res_pred)
    merge1 = pd.merge(input_bed, spliceai_res_pred_interpre, on=cols[:-1], how='left')
    merge2 = pd.merge(anno, merge1, on=['#Chr', 'Start', 'Stop', 'Ref', 'Call'], how='left')
    merge2.drop(columns=['#CHROM', 'POS', 'REF', 'ALT'], inplace=True)
    return merge2


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
    if in_format != 'excel':
        anno_df = pd.read_csv(bgi_anno, sep='\t')
    else:
        anno_df = pd.read_excel(bgi_anno, sheet_name=sheet)
    anno_df['#Chr'] = anno_df['#Chr'].astype('str')
    all_bed = pd.read_csv(pwd + '/' + out + '.bed2vcf.all.bed', sep='\t', dtype={'#Chr': str})
    if anno_df.empty:
        splice_pred = anno_df.copy()
        splice_pred['SpliceAI'] = []
        splice_pred['SpliceAI Pred'] = []
        splice_pred['SpliceAI Interpretation'] = []
    else:
        splice_pred = merge_spliceai_res(pwd, out, config_dic['skip_vcf_header_row'], all_bed, anno_df)
    if out_format == 'excel':
        splice_pred.to_excel(pwd + '/' + out + '.spliceai.xlsx', index=False)
    else:
        splice_pred.to_csv(pwd + '/' + out + '.spliceai.txt', index=False, sep='\t')
