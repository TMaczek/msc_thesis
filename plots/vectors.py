import matplotlib.pyplot as plt
import numpy as np


steps = [0, 2, 4, 8, 12, 16]
fps = 40

#j = 0

for j in steps:
    timestep = j*fps
    
    X, Y = np.meshgrid(np.arange(0, 50, 1), np.arange(0, 50, 1))
    
    U = np.empty([50, 50])
    V = np.empty([50, 50])
    
    with open('datafiles/6_2_velocitiesU.dat') as f:
        lines = f.readlines()
           
        for i in lines:
            segments = i.split(" ")
            if len(segments) > 1 and int(segments[0]) == timestep:
                x = int(segments[1])
                y = int(segments[2])
                val = float(segments[3])
                U[x, y] = val
                
    with open('datafiles/6_2_velocitiesV.dat') as f:
        lines = f.readlines()
           
        for i in lines:
            segments = i.split(" ")
            if len(segments) > 1 and int(segments[0]) == timestep:
                x = int(segments[1])
                y = int(segments[2])
                val = float(segments[3])
                V[x, y] = val
                    
    U = np.transpose(U)
    V = np.transpose(V)
    
    fig1, ax1 = plt.subplots(figsize=(10, 10))
    
    ax1.axis('equal')
    
    #ax1.set_title('Title')
    plt.gca().axes.get_xaxis().set_visible(False)
    plt.gca().axes.get_yaxis().set_visible(False)
    
    
    #Q = ax1.quiver(X, Y, U, V, angles='xy', scale_units='xy' , scale=3.5)
    Q = ax1.quiver(X[::3, ::3], Y[::3, ::3], U[::3, ::3], V[::3, ::3],
               pivot='mid', units='inches', scale=1.5)
    #strm = ax1.streamplot(X, Y, U, V, color=V, linewidth=4, cmap='autumn')

    fig1.tight_layout()
    
    fig1.savefig("6_2_vectors_" + str(j) + ".png")
    
    # mapki
    
    fig, ax = plt.subplots(figsize=(6, 5))
    plt.gca().invert_yaxis()
    im = ax.imshow(U, cmap='seismic', vmin=-4, vmax=4)
    im.axes.invert_yaxis()
    plt.gca().axes.get_xaxis().set_visible(False)
    plt.gca().axes.get_yaxis().set_visible(False)
    cbar = fig.colorbar(im, ax=ax, ticks=[-4, 0, 4])
    cbar.ax.tick_params(labelsize=20) 
    fig.tight_layout()
    fig.savefig("6_2_velocitiesU_" + str(j) + ".png")
    
    
    
    fig, ax = plt.subplots(figsize=(6, 5))
    plt.gca().invert_yaxis()
    im = ax.imshow(V, cmap='seismic', vmin=-4, vmax=4)
    im.axes.invert_yaxis()
    plt.gca().axes.get_xaxis().set_visible(False)
    plt.gca().axes.get_yaxis().set_visible(False)
    cbar = fig.colorbar(im, ax=ax, ticks=[-4, 0, 4])
    cbar.ax.tick_params(labelsize=20) 

    fig.tight_layout()
    fig.savefig("6_2_velocitiesV_" + str(j) + ".png")
