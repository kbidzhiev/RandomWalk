
###########################################
set terminal gif animate delay 5

dt_from_cpp = '0.5';




dt = dt_from_cpp + 0; ## this will convert string into a double

set output 'Sz_montecarlo_dt' .dt_from_cpp. '.gif' # . is a concatination of strings
FILE = '../Data/Sz_dt' . dt_from_cpp . '.dat'

stats FILE nooutput

set key c r

set xrange [-2:2]
set yrange [0:1]
set ylabel "Sz"
set xlabel "x"
percentile="P5 P10 P20 P25 P50 P75"
cir=7;

col1="black";
col2="red";
col3="orangel";
col4="#04B431";
col5="#3399FF";
col6="dark-violet";


f(x) = x;
do for [i=1:int(STATS_blocks)-1: floor(1.0/dt)] {
    plot FILE index (i) u ($1/sqrt(i)):($2*sqrt(i)) w p pt cir ps 1.5 lt rgb "red" title columnheader
}


