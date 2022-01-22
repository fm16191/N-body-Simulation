# set title "US immigration from Northern Europe\nPlot selected data columns as histogram of clustered boxes"
# set auto x
# #set yrange [0:1]
# set style data histogram
# set style histogram cluster gap 1
# set style fill solid border -1
# set boxwidth 0.9
# set xtic rotate by -45 scale 0
# #set bmargin 10
# set style histogram gap 5

set output "bench.svg"

# binwidth=5
# set boxwidth binwidth
# set style fill solid 1.0 border lt -1
# set style data histograms
# #bin(x,width)=width*floor(x/width)
# bin(x,width)=width*floor(x/width) + width/2.0
# # (bin($4,binwidth)):(1.0)

set boxwidth 0.9 relative
set style histogram rows
#set style data histograms
#set style fill solid 1.0 border -1
plot 'precisions/avg_p_g.dat' using 4, 'precisions/avg_p_cl.dat' using 4, 'precisions/avg_p_amd.dat' using 4

plot 'precisions/avg_p_g.dat' using 4:1 smooth freq with boxes,\
'precisions/avg_p_cl.dat' using 4:1 smooth freq with boxes,\
'precisions/avg_p_amd.dat' using 4:1 smooth freq with boxes

# plot 'precisions/avg_p_g.dat' using 1:xtic(4) title columnheader(2)
# using 1:xtic(1) ti col, '' u 12 ti col, '' u 13 ti col, '' u 14 ti col






#
# set datafile fortran
# set term svg
# set size ratio 0.8
# set key right center Right
# set key samplen 2 spacing .8 height 3 font ",10"
# #set logscale y
# set xtics scale default 0,1
# #set xrange [0:]
# set xlabel offset 0,0.3
# set xlabel "size of the n matrix dimension"
# set ylabel offset 2 "Average execution time (sec)"
# set title font "Helvetica Bold, 10" #offset -5
# set style line 1 lt 1 linecolor rgb "#0d42ba" lw 2 pt 1
# set style line 2 lt 1 linecolor rgb "#0d98ba" lw 2 pt 1
# set style line 3 lt 1 linecolor rgb "#0dba86" lw 2 pt 1
# 
# set style line 4 lt 1 linecolor rgb "red" lw 2 pt 1
# set style line 5 lt 1 linecolor rgb "#00FF01" lw 2 pt 1
# 
# set style line 9 lt 1 linecolor rgb "red" lw 0.5 pt 1 # red dispersion
# 
# set ytics nomirror
# #set y2tics border scale default
# #set y2range [0:]
# #set logscale y2
# #set y2label offset -2 "Average speed difference"
# 
# set output "bench.svg"
# set title "Difference of execution time in logarithmic scale according to the size of the matrix"
# plot "precisions/avg_p_g.dat" using 1:4 ls 1 title "avg precision order 1" with lines, "precisions/avg_p_g.dat" using 1:4:3:2 ls 9 notitle with errorbars #,\
# # ".dat" using 1:8 ls 2 title "mCSRv" with lines, ".dat" using 1:8:10:9 ls 9 notitle with errorbars,\
# # ".dat" using 1:5 ls 3 title "csmtCSR" with lines, ".dat" using 1:5:7:6 ls 9 title "dispersion" with errorbars,\
# # ".dat" using 1:($8/$2) smooth bezier ls 4 title "speed mCSRv / ref" with lines axis x1y2
# 
# # set xtics scale default 0,1
# # set xlabel "sparse matrix density"
# # set output "_bench_density.svg"
# # set title "Difference of execution time in logarithmic scale according to the density of the matrix"
# # plot "data/CSR_density.dat" using 1:2 ls 1 title "reference mult" with lines, "data/CSR_density.dat" using 1:2:4:3 ls 9 notitle with errorbars,\
# # "data/CSR_density.dat" using 1:8 ls 2 title "mCSRv" with lines, "data/CSR_density.dat" using 1:8:10:9 ls 9 notitle with errorbars,\
# # "data/CSR_density.dat" using 1:5 ls 3 title "csmtCSR" with lines, "data/CSR_density.dat" using 1:5:7:6 ls 9 title "dispersion" with errorbars,\
# # "data/CSR_density.dat" using 1:($8/$2) smooth bezier ls 4 title "speed mCSRv / ref" with lines axis x1y2
# # 
# # set key right bottom Right
# # set xtics scale default 0,100
# # set xlabel "size of the n matrix dimension"
# # unset logscale y
# # set ytics mirror
# # unset y2tics
# # unset y2label
# # set ylabel "Average precision difference"
# # 
# # set output "data/CSR_bench_precision.svg"
# # set title "Difference of precision according to the size of the matrix"
# # plot "data/CSR.dat" using 1:11 ls 5 title "precision difference" with lines, "data/CSR.dat" using 1:11:13:12 ls 9 title "dispersion" with errorbars
# # 
# # #i, mean_t1, min(t(1,:)), max(t(1,:)), mean_t2, min(t(2,:)), max(t(2,:)), mean_norm, min(n(1,:)), max(n(1,:))