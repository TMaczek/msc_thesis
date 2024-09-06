from matplotlib import pyplot as plt
import numpy as np

if __name__ == '__main__':
  print('hello')

   # komorki
  with open('datafiles/komorki.dat') as f:
    lines = f.readlines()
    x = [int(line.split()[0]) for line in lines]
    y = [int(line.split()[1]) for line in lines]
    v = [(line.split()[2]) for line in lines]
  
  asp = (max(y) + 1) / (max(x) + 1)
  print(asp)
  plt.figure(figsize=(10, 10*asp))

  #scatter = plt.imshow(x, y, c=v, cmap='gnuplot', s=1000, marker='s',  alpha=0.8)
  t = [ [0]*200 for i in range(200)]

  for a, b, c in zip(x, y, v):
      t[199-b][a] = float(c)
      
  
  plt.imshow(t, cmap='gnuplot')

   # Adding color bar
  plt.colorbar()
  
   # Show plot
   # plt.grid(True)
  plt.title("Predkosc pozioma 1 iter")
  plt.axis('off')
  plt.show()
  
  # czastki
  # with open('datafiles/czastki.dat') as f:
  # lines = f.readlines()
  # x = [line.split()[0] for line in lines]
  # y = [line.split()[1] for line in lines]
  # data = np.loadtxt('datafiles/czastki.dat')
  # x = data[:, 0]
  # y = data[:, 1]
  # plt.figure()
  # plt.plot(x, y, 'b.', markersize=1)

  # Show plot
  # plt.grid(True)

  # # later: wziac rozmiary pudla z maxow w poprzednim lub przekazac w pliku w 1 linii
  # plt.xlim([0, 200])
  # plt.ylim([0, 200])

  # plt.gca().set_aspect('equal')
  # plt.title("Czastki znaczone")
  # # plt.tick_params(
  # # 	axis='x',  # changes apply to the x-axis
  # # 	which='both',  # both major and minor ticks are affected
  # # 	bottom=False,  # ticks along the bottom edge are off
  # # 	top=False,  # ticks along the top edge are off
  # # 	labelbottom=False)  # labels along the bottom edge are off

  # # plt.tick_params(
  # # 	axis='y',  # changes apply to the x-axis
  # # 	which='both',  # both major and minor ticks are affected
  # # 	left=False,  # ticks along the bottom edge are off
  # # 	right=False,  # ticks along the top edge are off
  # # 	labelleft=False)  # labels along the bottom edge are off

  # plt.show()
