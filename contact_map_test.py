'''
Calculate contact map between two groups of atoms (residues)
based on a specified cut-off distance using multiple MD simulation
trajectories (plotting is not included).   
'''

import parmed as pmd
import MDAnalysis
from MDAnalysis.analysis import distances
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

topo = "data/amber.prmtop"
traj = "data/md.1.mdcrd.nc"
top = pmd.load_file(topo)
u = MDAnalysis.Universe(top, traj)

group_1 = u.select_atoms('resid 169:185 and name CA')
group_2 = u.select_atoms('resid 169:185 and name CA')
# group_1 = u.select_atoms('resid 169:185')
# group_2 = u.select_atoms('resid 169:185')
n1 = len(group_1)
n2 = len(group_2)
print("selected group N1, N2 =", n1, n2)

contact_sum = np.zeros((n1, n2))
max_distance = 10.0

n_frames = 0

for ts in u.trajectory:
    gp1 = group_1.positions
    gp2 = group_2.positions
    ts_dist = distances.distance_array(gp1, gp2, box=u.dimensions)
    ts_dist[ts_dist < max_distance] = 1
    ts_dist[ts_dist > max_distance] = 0
    contact_sum = ts_dist + contact_sum
    n_frames += 1

print("N_frames =", n_frames)
contact_ratio = contact_sum / n_frames
print("contact ratio =", contact_ratio)

plt.imshow(contact_ratio, cmap='hot', interpolation='nearest')
plt.colorbar()
plt.title('Contact Map Heatmap (original, all atoms)')
plt.show()

plt.savefig('heatmap_test_HVR_all_atoms.png')