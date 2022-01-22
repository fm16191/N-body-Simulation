set style data histograms
set boxwidth 0.9
set style fill solid 1.00 border 0
set style histogram errorbars gap 2 
#set style histogram errorbars linecolor rgb "#0d42ba" lw 1 lc
set style fill solid
set size ratio 0.8
set term svg
set output "benchs/avg_perf.svg"
set title "Average performance in GFlops depending of the compiler"
stats [1:2] 'benchs/b_g.dat' using 4 prefix "A" nooutput
stats 'benchs/b_cl.dat' using 4 prefix "B" nooutput
stats 'benchs/b_amd.dat' using 4 prefix "C" nooutput
set yrange [(A_min-2):(A_max+2)]
plot "" u (A_mean):(A_min):(A_max):xtic(1) title 'GCC'
#,\
#"<(sed -n '1,1p' benchs/b_cl.dat)" using (B_mean):(B_min):(B_max):xtic(1) title 'Clang' ,\
#"<(sed -n '1,1p' benchs/b_amd.dat)" using (C_mean):(C_min):(C_max):xtic(1) title 'AOCC'