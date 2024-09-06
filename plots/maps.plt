
N = 6
V = 1

n = 21

plottype = "velocitiesU"

reset
datafile = sprintf("datafiles/%d_%d_%s.dat", N, V, plottype)

set terminal pngcairo size 580,500 enhanced

fps = 40


set palette defined (-4 "blue", 0 "white", 4 "red")


do for [v=0:n]{

	j = v * fps
	set output sprintf('%d_%d_plot_%s_%d.png', N, V, plottype, v)
	set view map # widok z gory
	set size ratio -1
	set cbr [-4:4]
	set xr [0:50] 
	set yr [0:50]
	unset border
	unset xtics
	unset ytics
	plot sprintf('< grep "^%d " %s', j, datafile) using 2:3:4  title "" with image
}
