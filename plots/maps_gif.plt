# ------------------------------------------------------------------
# Rysowanie gifow
# ------------------------------------------------------------------

reset

N = 6
V = 1

plottype = "velocitiesU"
datafile = sprintf("datafiles/%d_%d_%s.dat", N, V, plottype)

set term gif size 500,500 animate delay 5 enhanced 

set palette defined(-4 "midnight-blue", -3 "purple", -2 "skyblue", -1 "gray", 0 "black", 1 "goldenrod", 2 "forest-green", 3 "greenyellow", 4 "orange")

set output "anim2.gif"
n=866 #liczba klatek
set view map # widok z gory
set size ratio -1
set xr [0:50] 
set yr [0:50]
# unset colorbox
set cbrange [-4:4] 




do for [j=0:n:2] {
  plot sprintf('< grep "^%d " %s', j, datafile) using 2:3:4 title sprintf("t=%i",j/40) with image

}

