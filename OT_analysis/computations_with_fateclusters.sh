#!/bin/bash

prefix=eps005
#wot optimal_transport --matrix RGC_batch_uncorr_190910_small.mtx --cell_days cell_day.txt --growth_iters 3 --epsilon 0.005 --out tmaps/${prefix}
wot trajectory --tmap tmaps/${prefix} --cell_set cell_sets_fate_clusters.gmt --day 0 --out tmaps/wot_${prefix}_E13_traj_fate_clusters.txt --verbose
wot trajectory --tmap tmaps/${prefix} --cell_set cell_sets_fate_clusters.gmt --day 1 --out tmaps/wot_${prefix}_E14_traj_fate_clusters.txt --verbose
wot trajectory --tmap tmaps/${prefix} --cell_set cell_sets_fate_clusters.gmt --day 3 --out tmaps/wot_${prefix}_E16_traj_fate_clusters.txt --verbose
wot trajectory --tmap tmaps/${prefix} --cell_set cell_sets_fate_clusters.gmt --day 6 --out tmaps/wot_${prefix}_P0_traj_fate_clusters.txt --verbose
wot trajectory --tmap tmaps/${prefix} --cell_set cell_sets_fate_clusters.gmt --day 11 --out tmaps/wot_${prefix}_P5_traj_fate_clusters.txt --verbose
wot trajectory --tmap tmaps/${prefix} --cell_set cell_sets_fate_clusters.gmt --day 16 --out tmaps/wot_${prefix}_P56_traj_fate_clusters.txt --verbose

wot fates --tmap tmaps/${prefix} --cell_set cell_sets_fate_clusters.gmt --day 0 --out tmaps/wot_${prefix}_E13_fate_fate_clusters.txt --verbose
wot fates --tmap tmaps/${prefix} --cell_set cell_sets_fate_clusters.gmt --day 1 --out tmaps/wot_${prefix}_E14_fate_fate_clusters.txt --verbose
wot fates --tmap tmaps/${prefix} --cell_set cell_sets_fate_clusters.gmt --day 3 --out tmaps/wot_${prefix}_E16_fate_fate_clusters.txt --verbose 
wot fates --tmap tmaps/${prefix} --cell_set cell_sets_fate_clusters.gmt --day 6 --out tmaps/wot_${prefix}_P0_fate_fate_clusters.txt --verbose 
wot fates --tmap tmaps/${prefix} --cell_set cell_sets_fate_clusters.gmt --day 11 --out tmaps/wot_${prefix}_P5_fate_fate_clusters.txt --verbose
wot fates --tmap tmaps/${prefix} --cell_set cell_sets_fate_clusters.gmt --day 16 --out tmaps/wot_${prefix}_P56_fate_fate_clusters.txt --verbose


