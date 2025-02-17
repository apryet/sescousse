#!/bin/zsh

python gen_hist.py --tpl_dir=drn0 --sim_dir=histo_drn0 > histo_drn0.out &
python gen_hist.py --tpl_dir=drn40 --sim_dir=histo_drn40 > histo_drn40.out &
python gen_hist.py --tpl_dir=drn110 --sim_dir=histo_drn110 > histo_drn110.out &
python gen_hist.py --tpl_dir=drn150 --sim_dir=histo_drn150 > histo_drn150.out &
python gen_hist.py --tpl_dir=drn200 --sim_dir=histo_drn200 > histo_drn200.out &

python gen_hist.py --tpl_dir=drn_alt_ld_40 --sim_dir=histo_drn_alt_ld_40 > histo_drn_alt_ld_40.out &
python gen_hist.py --tpl_dir=drn_alt_ld_110 --sim_dir=histo_drn_alt_ld_110 > histo_drn_alt_ld_110.out &
python gen_hist.py --tpl_dir=drn_alt_hd_40 --sim_dir=histo_drn_alt_hd_40 > histo_drn_alt_hd_40.out &
python gen_hist.py --tpl_dir=drn_alt_hd_110 --sim_dir=histo_drn_alt_hd_110 > histo_drn_alt_hd_110.out &


for i in {1..11}; do
    python gen_prosp.py --tpl_dir=drn110 --cm=${i} > prosp_$(printf "%02d" "$i").out &

