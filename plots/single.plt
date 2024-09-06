# ------------------------------------------------------------------
# Rysowanie obrazka!
# ------------------------------------------------------------------

N = 6
V = 2

n = 5


reset
datafile = sprintf("datafiles/%d_%d_particles.dat", N, V)
map = sprintf("datafiles/%d_%d_map.dat", N, V)


set terminal pngcairo size 600,600 enhanced


fps = 40

set palette defined (0 "steelblue", 1 "forest-green", 2 "dark-pink", 3 "dark-yellow", 4 "orange", \
					5 "skyblue", 6 "dark-grey")
unset colorbox
set cbrange [0:6] 


do for [v=0:n]{

	j = v * fps
	set output sprintf('%d_%d_plot_%d.png', N, V, v)
	set view map # widok z gory
	set size ratio -1
	set cbr [0:]
	set xr [0:100] 
	set yr [0:100]
	unset border
	unset xtics
	unset ytics

	 
  	plot datafile i j u 2:3:4  title sprintf("") with points palette pointtype 7 pointsize 0.8, \
		map using (($1+$3)/2):(($2+$4)/2):1:3:2:4 title "" with boxxyerr fill fc "black"  fs solid 1.0

}






