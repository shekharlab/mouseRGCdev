#!/bin/bash

prefix=eps005
wot optimal_transport --matrix RGC_batch_uncorr_190910_small.mtx --cell_days cell_day.txt --growth_iters 3 --epsilon 0.005 --out tmaps/${prefix}
wot trajectory --tmap tmaps/${prefix} --cell_set cell_sets.gmt --day 0 --out tmaps/wot_${prefix}_E13_traj.txt --verbose
wot trajectory --tmap tmaps/${prefix} --cell_set cell_sets.gmt --day 1 --out tmaps/wot_${prefix}_E14_traj.txt --verbose
wot trajectory --tmap tmaps/${prefix} --cell_set cell_sets.gmt --day 3 --out tmaps/wot_${prefix}_E16_traj.txt --verbose
wot trajectory --tmap tmaps/${prefix} --cell_set cell_sets.gmt --day 6 --out tmaps/wot_${prefix}_P0_traj.txt --verbose
wot trajectory --tmap tmaps/${prefix} --cell_set cell_sets.gmt --day 11 --out tmaps/wot_${prefix}_P5_traj.txt --verbose
wot trajectory --tmap tmaps/${prefix} --cell_set cell_sets.gmt --day 16 --out tmaps/wot_${prefix}_P56_traj.txt --verbose

wot fates --tmap tmaps/${prefix} --cell_set cell_sets.gmt --day 0 --out tmaps/wot_${prefix}_E13_fate.txt --verbose
wot fates --tmap tmaps/${prefix} --cell_set cell_sets.gmt --day 1 --out tmaps/wot_${prefix}_E14_fate.txt --verbose
wot fates --tmap tmaps/${prefix} --cell_set cell_sets.gmt --day 3 --out tmaps/wot_${prefix}_E16_fate.txt --verbose 
wot fates --tmap tmaps/${prefix} --cell_set cell_sets.gmt --day 6 --out tmaps/wot_${prefix}_P0_fate.txt --verbose 
wot fates --tmap tmaps/${prefix} --cell_set cell_sets.gmt --day 11 --out tmaps/wot_${prefix}_P5_fate.txt --verbose
wot fates --tmap tmaps/${prefix} --cell_set cell_sets.gmt --day 16: --out tmaps/wot_${prefix}_P56_fate.txt --verbose


