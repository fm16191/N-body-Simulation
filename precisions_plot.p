set style data histograms
set boxwidth 0.9
set style fill solid 1.00 border 0
set style histogram errorbars gap 2 lw 1 
#set style histogram errorbars linecolor rgb "#0d42ba" lw 1 lc
set style fill solid
set size ratio 0.8

set logscale y

set key left Right
set key samplen 2 spacing .8 height 3 font ",10"

set term svg
set output "avg_precision.svg"
set title "Average precision loss depending of the compiler"
plot "precisions.dat" every 3:3 using 5:4:3:xtic(1) title "gcc",\
"" skip 1 every 3:3  using 5:4:3:xtic(1) title "icc",\
"" skip 2 every 3:3  using 5:4:3:xtic(1) title "icx",\