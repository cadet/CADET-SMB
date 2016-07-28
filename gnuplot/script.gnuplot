#!/usr/bin/gnuplot -persist

set terminal postscript eps color solid linewidth 2 "Helvetica" 20

set output "trajectory_raffinate.eps"

set multiplot

set xrange [0:100000]
set yrange [-0.0005:0.003]

set xtics add ("0" 0, \
"20" 20000, \
"40" 40000, \
"60" 60000, \
"80" 80000, \
"100" 100000)

set xlabel "Switches"
set ylabel "Concentration [mMol]"

set grid 
set key box

plot "trajectory.dat" using 1 with lines title "comp 1", "trajectory.dat" using 2 with lines title "comp 2", "trajectory.dat" using 3 with lines title "comp 3"

set size 0.5,0.4
set origin 0.4,0.12

unset xlabel
unset ylabel

set xrange [96000:98000]
set yrange [0:0.0025]

set xtics 1000
set ytics 0.001

set xtics add ("96" 96000, \
"97" 97000, \
"98" 98000)

unset grid
unset key

plot "trajectory.dat" using 1 with lines, "trajectory.dat" using 2 with lines, "trajectory.dat" using 3 with lines

unset multiplot

set output


#=========================================================#
unset size
unset origin
unset xlabel
unset ylabel
unset xtics
unset ytics

set output "trajectory_extract.eps"

set multiplot

set xrange [0:100000]
set yrange [-0.001:0.0025]

set xtics 20000
set ytics 0.0005
set xtics add ("0" 0, \
"20" 20000, \
"40" 40000, \
"60" 60000, \
"80" 80000, \
"100" 100000)

set xlabel "Switches"
set ylabel "Concentration [mMol]"

set grid 
set key box

plot "trajectory.dat" using 4 with lines title "comp 1", "trajectory.dat" using 5 with lines title "comp 2", "trajectory.dat" using 6 with lines title "comp 3"

#---------------------------------------------------------#
set size 0.5,0.4
set origin 0.4,0.12

unset xlabel
unset ylabel

set xrange [96000:98000]
set yrange [0:0.0025]

set xtics 1000
set ytics 0.001

set xtics add ("96" 96000, \
"97" 97000, \
"98" 98000)

unset grid
unset key

plot "trajectory.dat" using 4 with lines, "trajectory.dat" using 5 with lines, "trajectory.dat" using 6 with lines

unset multiplot

set output


#=========================================================#
unset size
unset origin
unset xlabel
unset ylabel
unset xtics
unset ytics

set output "chromatogram.eps"

set xrange [0:8000]
set yrange [0:0.0025]

set xtics 2000
set xtics ("" 0, \
"Zone IV" 1000, \
"" 2000, \
"Zone III" 3000, \
"" 4000, \
"Zone II" 5000, \
"" 6000, \
"Zone I" 7000, \
"" 8000)

set ytics 0.0005

set grid

set ylabel "Concentration [mMol]"

set key box

plot "chromatogram.dat" using 1 with lines title "comp 1", "chromatogram.dat" using 2 with lines title "comp 2", "chromatogram.dat" using 3 with lines title "comp 3"

set output

set terminal wxt
