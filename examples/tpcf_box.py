from contrast.box import tpcf_monopole, projected_tpcf
import numpy as np
import matplotlib.pyplot as plt

# read simulation data
positions = np.genfromtxt("data/R101_S15.csv",
    skip_header=1, delimiter=",", usecols=(0, 1, 2)
)

box_size = 2000
rbins = np.linspace(0, 10, 11)
rbins_c = 0.5*(rbins[1:] + rbins[:-1])

# delta_r = tpcf_monopole(positions1=positions, box_size=box_size, rbins=rbins)
sigma_r = projected_tpcf(positions1=positions, box_size=box_size, rbins=rbins)

# fig, ax = plt.subplots()

# ax.plot(rbins_c, rbins_c ** 2 * delta_r)
# plt.show()


