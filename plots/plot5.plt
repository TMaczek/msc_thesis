# ------------------------------------------------------------------
# Rysowanie gifow
# ------------------------------------------------------------------

N = 8
V = 1
n= 565 #liczba klatek


reset
set term gif size 600, 300 animate delay 5 enhanced
set output sprintf("animations/%d_%d_anim.gif", N, V)
set view map # widok z gory
set size ratio -1
set cbr [0:]
set xr [0:200] 
set yr [0:100]
unset border
unset xtics
unset ytics

set palette defined (0 "steelblue", 1 "forest-green", 2 "dark-pink", 3 "dark-yellow", 4 "orange", \
					5 "skyblue", 6 "dark-grey")
unset colorbox
set cbrange [0:6] 

do for [j=0:n:2] {
	file = sprintf("datafiles/%d_%d_particles.dat", N, V)
	map = sprintf("datafiles/%d_%d_map.dat", N, V)
  
	plot file i j u 2:3:4  title sprintf("t=%i",j/40) at screen 0.9, 0.98 with points palette pointtype 7 pointsize 0.2 , \
		 map using (($1+$3)/2):(($2+$4)/2):1:3:2:4 title "" with boxxyerr fill fc "black"  fs solid 1.0

}

