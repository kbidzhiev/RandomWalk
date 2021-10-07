
###########################################
set terminal gif animate delay 5

dt_from_cpp = '0.1';




dt = dt_from_cpp + 0; ## this will convert string into a double

#set output 'Sz_montecarlo_dt' .dt_from_cpp. '.gif' # . is a concatination of strings
#FILE = '../Data/Sz' . dt_from_cpp . '.dat'

set output 'Sz_montecarlo.gif' # . is a concatination of strings
FILE = '../Data/Sz.dat'

stats FILE nooutput

set key t r

set xrange [-4:4]
set yrange [0:0.5]
set ylabel "Sz * sqrt(2 prob * n)"
set xlabel "x/sqrt(2 prob * n)"
percentile="P5 P10 P20 P25 P50 P75"
cir=7;

col1="black";
col2="red";
col3="orangel";
col4="#04B431";
col5="#3399FF";
col6="dark-violet";



#floor(1.0/dt)
step = 1;


coeff= 2*dt;

do for [i=1:int(STATS_blocks)-1: step] {
    plot exp(-x**2/2)/sqrt(2*pi) lt rgb "black" ti "exp(-x^2/2)/sqrt(2*pi)", FILE index (i) u ($1/sqrt(coeff*i)):($2*sqrt(coeff*i)) w p pt cir ps 1.5 lt rgb "red" title columnheader
}

!eog Sz_montecarlo.gif &


