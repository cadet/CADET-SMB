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
set size 1,0.75
set origin 0,0
set size ratio 0.55
set lmargin 6.5
set rmargin 1.3
set tmargin 1.5
set bmargin 1.5
#------------------------------------------------------------------------------


set output "chromatogram.tex"

set macros
set format y "$%2.1t$"

set xrange [0:320]
set yrange [0:0.003]

set xtics 40
set xtics ("" 0, "Zone I" 40, "" 80, "II" 120, "" 160, "III" 200, "" 240, "IV" 280, "" 320)

set ytics 0.001

set grid

set ylabel "Concentration [mol]" offset 1.5

set key nobox left
set key width -2.8

plot "chromatogram.dat" using 1 with lines ls 2 lw 3.5 title '$c_1(\tau, z)$', \
     "chromatogram.dat" using 2 with lines ls 5 lw 3.5 title '$c_2(\tau, z)$'

set output

