import numpy as np
y=np.arange(39.5, 41.5, 0.5)
x=np.arange(80,87, 0.5)
mesh = np.meshgrid(x,y)
xx = mesh[0].flatten()
yy = mesh[1].flatten()
for i in np.arange(len(xx)):
    print(xx[i], yy[i], "0,")

