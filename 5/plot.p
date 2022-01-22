set style data histograms
set output "test.png"
#plot './colors.dat' using 2:xtic(1)
plot 'colors.dat' using 2:xtic(1) title 'Values by Color'

plot 'precisions/avg_p.dat' using 4:xtic(1) title 'perf'