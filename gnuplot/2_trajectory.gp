#!/usr/bin/gnuplot -persist

#------------------------------------------------------------------------------
# header

set terminal epslatex color solid standalone 'ptm' 14 \
    header '\usepackage{xcolor, amsmath}'

set style line 1 lt 1 lc rgb "black" lw 3
set style line 2 lt 2 lc rgb "red" lw 3
set style line 3 lt 5 lc rgb "dark-blue" lw 6 
set style line 4 lt 2 lc rgb "orange" lw 6 
set style line 5 lt 3 lc rgb "dark-green" lw 3
set style line 6 lt 3 lc rgb "yellow" lw 3

set style line 10 lt 1 lc rgb "black" lw 6

set style line 11 lt 1 pt 5 lc rgb "black" lw 4.5
set style line 12 lt 2 pt 7 lc rgb "red" lw 4.5
set style line 13 lt 8 pt 9 lc rgb "dark-blue" lw 4.5
set style line 14 lt 8 pt 13 lc rgb "dark-orange" lw 4.5


set pointsize 2
set size 1,0.95
set origin 0,0
set size ratio 0.65
set lmargin 6.5
set rmargin 1.3
set tmargin 1.5
set bmargin 3
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# plot trajectories at raffinate port and extract port respectively

set output "trajectory_raffinate.tex"

set macros
set format y "$%2.1t$"

set xrange [0:60000]
set yrange [-0.0005:0.0015]

set ytics 0.001
set xtics 10000
set xtics add ("0" 0, "10" 10000, "20" 20000, "30" 30000, "40" 40000, "50" 50000, "60" 60000)

set xlabel "Iteration" offset 1.5
set ylabel "Concentration [mol]" offset 1.5

set grid 
set key nobox left 
set key width -5

plot "2_trajectory.dat" using 1 with lines ls 2 lw 3 title '$c_{\text{out},1}(t)$', \
     "2_trajectory.dat" using 2 with lines ls 5 lw 3 title '$c_{\text{out},2}(t)$'

set output

set output "trajectory_extract.tex"

set macros
set format y "$%2.1t$"

set xrange [0:60000]
set yrange [-0.0005:0.0015]

set xtics 10000
set ytics 0.001

set xtics add ("0" 0, "10" 10000, "20" 20000, "30" 30000, "40" 40000, "50" 50000, "60" 60000)

set xlabel "Iteration"
set ylabel "Concentration [mol]" offset 1.5

set grid 
set key nobox left
set key width -5

plot "2_trajectory.dat" using 3 with lines ls 2 lw 3 title '$c_{\text{out},1}(t)$', \
     "2_trajectory.dat" using 4 with lines ls 5 lw 3 title '$c_{\text{out},2}(t)$'

set output

