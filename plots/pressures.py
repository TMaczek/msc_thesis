import matplotlib.pyplot as plt
import numpy as np

steps = [0, 3, 6, 10, 15, 20]
fps = 40

for j in steps:
    timestep = j*fps
        
    P = np.empty([100, 100])
    
    with open('datafiles/pressures.dat') as f:
        lines = f.readlines()
           
        for i in lines:
            segments = i.split(" ")
            if len(segments) > 1 and int(segments[0]) == timestep:
                x = int(segments[1])
                y = int(segments[2])
                val = float(segments[3])
                P[x, y] = val
                
                    
    P = np.transpose(P)
        
    # mapki
    
    fig, ax = plt.subplots(figsize=(6, 5))
    plt.gca().invert_yaxis()
    im = ax.imshow(P, vmin=0, vmax=15, cmap='inferno')
    im.axes.invert_yaxis()
    plt.gca().axes.get_xaxis().set_visible(False)
    plt.gca().axes.get_yaxis().set_visible(False)
    cbar = fig.colorbar(im, ax=ax, ticks=[0, 5, 10, 15])
    cbar.ax.tick_params(labelsize=20) 
    fig.tight_layout()
    fig.savefig("3_2_pressures_" + str(j) + ".png")
    