#!/bin/bash
/zfssz4/B2C_RD_P2/PMO/fangzhonghai/software/anaconda3/envs/python37/bin/python /zfssz4/B2C_RD_P2/PMO/fangzhonghai/work/lfx/spliceai/spliceai_in_wgs_pipe.py \
-b test.xlsx -p /zfssz4/B2C_RD_P2/PMO/fangzhonghai/work/lfx/spliceai -o test \
-c /zfssz4/B2C_RD_P2/PMO/fangzhonghai/work/lfx/spliceai/spliceai.yaml --sheet intron --out_format excel
