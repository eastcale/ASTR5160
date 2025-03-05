from HW.HW0 import *

x = np.linspace(0, 10)

y = line(x, 1, 0)

fig, ax = plt.subplots()

ax.scatter(x, y)

ticks(ax, x, y, np.size(x)*[0], np.size(y)*[0])

save_show(fig, 'testline.png')