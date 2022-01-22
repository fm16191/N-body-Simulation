set style data histograms
set boxwidth 0.9
set style fill solid 1.00 border 0
set style histogram #errorbars gap 2 lw 1
set style fill solid
set size ratio 0.8
set term svg
set output "precisions/avg_precision_loss.svg"
set title "Average precision loss, depending of the compiler"
#plot 'precisions/avg_p.dat' using 4:2:3:xtic(1) title 'at order 1', '' using 5:2:3:xtic(1) title 'at order 2' linecolor 'blue'
plot 'precisions/avg_p.dat' using 4:xtic(1) title 'at order 1', '' using 5:xtic(1) title 'at order 2' linecolor 'blue'