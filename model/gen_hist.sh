#!/bin/zsh
python gen_hist.py --tpl_dir=drn_alt_ld_40 --sim_dir=histo_drn_alt_ld_40 > histo_drn_alt_ld_40.out &
python gen_hist.py --tpl_dir=drn_alt_ld_110 --sim_dir=histo_drn_alt_ld_110 > histo_drn_alt_ld_110.out &
python gen_hist.py --tpl_dir=drn_alt_hd_40 --sim_dir=histo_drn_alt_hd_40 > histo_drn_alt_hd_40.out &
python gen_hist.py --tpl_dir=drn_alt_hd_110 --sim_dir=histo_drn_alt_hd_110 > histo_drn_alt_hd_110.out &


cd histo_drn_alt_dd_40 
python forward_run.py &
cd ..

cd histo_drn_alt_dd_110 
python forward_run.py &
cd ..

cd histo_drn_alt_hd_40 
python forward_run.py &
cd ..

cd histo_drn_alt_hd_110
python forward_run.py &
cd ..

/home/apryet/model/sescousse/sim/histo



for i in {1..1}; do
    python gen_prosp.py --tpl_dir=drn110 --cm=${i} > prosp_$(printf "%02d" "$i").out &

