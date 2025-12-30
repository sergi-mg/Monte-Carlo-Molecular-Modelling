# -*- coding: utf-8 -*-
"""
@author: Sergi Martinez Galindo
"""

import matplotlib.pyplot as plt
import numpy as np


#%%
#reading the data
data=np.loadtxt("results_MC_2.10_100.txt")
time=np.arange(10,10**6+1,10)

# Figure 1: Energy 
fig1, ax1 = plt.subplots()

ax1.plot(time, data[1:, 0], color="#D4A017")
ax1.set_xscale("log")
ax1.grid(which="major", linestyle='-', linewidth=0.7)
ax1.grid(which="minor", linestyle='--', linewidth=0.3)
ax1.set_xlabel(r"$t$ (MCS)")
ax1.set_ylabel(r"$E$")
ax1.set_ylim([-1.75, -1.2])
ax1.set_xlim([10, 10**6])

fig1.savefig("time_series_E_2_1.pdf", bbox_inches="tight")


# Figure 2: Magnetization
fig2, ax2 = plt.subplots()

ax2.plot(time, data[1:, 1], color="#0A2A43")
ax2.set_xscale("log")
ax2.grid(which="major", linestyle='-', linewidth=0.7)
ax2.grid(which="minor", linestyle='--', linewidth=0.3)
ax2.set_xlabel(r"$t$ (MCS)")
ax2.set_ylabel(r"$M$")
ax2.set_ylim([-1, 0])
ax2.set_xlim([10, 10**6])

fig2.savefig("time_series_M_2_1.pdf", bbox_inches="tight")

plt.show()
