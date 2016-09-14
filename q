[1mdiff --git a/structopt/cluster/individual/generators/fcc.py b/structopt/cluster/individual/generators/fcc.py[m
[1mindex dd6a931..5aa190c 100644[m
[1m--- a/structopt/cluster/individual/generators/fcc.py[m
[1m+++ b/structopt/cluster/individual/generators/fcc.py[m
[36m@@ -211,7 +211,7 @@[m [mdef get_vector_angle(orientation=None, v=None, angle=None):[m
         v = np.array([1, 0, 0])[m
     elif orientation == '110':[m
         angle = np.pi / 4[m
[31m-        v = np.array([1, 0, 0])[m
[32m+[m[32m        v = np.array([0, 1, 0])[m
     elif orientation == '111':[m
         angle = np.arcsin(1.0 / 3.0 ** 0.5)[m
         v = np.array([-1.0 / 2.0 ** 0.5, 1.0 / 2.0 ** 0.5, 0])[m
[1mdiff --git a/structopt/common/individual/relaxations/STEM.py b/structopt/common/individual/relaxations/STEM.py[m
[1mindex 608fe25..7fc9047 100644[m
[1m--- a/structopt/common/individual/relaxations/STEM.py[m
[1m+++ b/structopt/common/individual/relaxations/STEM.py[m
[36m@@ -36,6 +36,7 @@[m [mclass STEM(structopt.common.individual.fitnesses.STEM):[m
         parameters.setdefault('rotation_grid', 10)[m
         parameters.setdefault('rotation_iterations', 2)[m
         parameters.setdefault('surface_moves', 10)[m
[32m+[m[32m        parameters.setdefault('filter_size', 1)[m
 [m
         super().__init__(parameters)[m
 [m
[36m@@ -139,7 +140,7 @@[m [mclass STEM(structopt.common.individual.fitnesses.STEM):[m
 [m
         return bonds[m
 [m
[31m-    def get_STEM_projection(self, individual):[m
[32m+[m[32m    def get_STEM_projection(self, individual, test=False):[m
         """Gets the projection of bonds from the brighest STEM image"""[m
 [m
         if self.target is None:[m
[36m@@ -147,11 +148,12 @@[m [mclass STEM(structopt.common.individual.fitnesses.STEM):[m
 [m
         parameters = self.parameters[m
         target = self.target[m
[32m+[m[41m        [m
 [m
         # Get a cutoff between maximum points in the STEM image based[m
         # on nearest neighbor distances[m
         cutoff = get_avg_radii(individual) * 2 * 1.1[m
[31m-        size = cutoff / 8 * parameters['resolution'][m
[32m+[m[32m        size = cutoff * parameters['resolution'] * parameters['filter_size'][m
 [m
         # Get a list of xy positions from analyzing local maxima in STEM image[m
         # as well as the position of a spot near the center of mass[m
[36m@@ -167,4 +169,25 @@[m [mclass STEM(structopt.common.individual.fitnesses.STEM):[m
         dists = np.linalg.norm(vecs, axis=1)[m
         vecs = vecs[((dists < cutoff) & (dists > 0))][m
 [m
[32m+[m[32m        if test:[m
[32m+[m[32m            import matplotlib.pyplot as plt[m
[32m+[m[32m            import matplotlib.cm as cm[m
[32m+[m
[32m+[m[32m            fig, ax = plt.subplots(num=1)[m
[32m+[m[32m            fig.colorbar(ax.pcolormesh(target, cmap=cm.viridis, linewidths=0))[m
[32m+[m[32m            ax.set_xlim((0, parameters['dimensions'][0] * parameters['resolution']))[m
[32m+[m[32m            ax.set_ylim((0, parameters['dimensions'][1] * parameters['resolution']))[m
[32m+[m
[32m+[m[32m            fig, ax = plt.subplots(num=2)[m
[32m+[m[32m            fig.colorbar(ax.pcolormesh(data_max, cmap=cm.viridis, linewidths=0))[m
[32m+[m[32m            ax.set_xlim((0, parameters['dimensions'][0] * parameters['resolution']))[m
[32m+[m[32m            ax.set_ylim((0, parameters['dimensions'][1] * parameters['resolution']))[m
[32m+[m
[32m+[m[32m            fig, ax = plt.subplots(num=3)[m
[32m+[m[32m            fig.colorbar(ax.pcolormesh(maxima, cmap=cm.viridis, linewidths=0))[m
[32m+[m[32m            ax.set_xlim((0, parameters['dimensions'][0] * parameters['resolution']))[m
[32m+[m[32m            ax.set_ylim((0, parameters['dimensions'][1] * parameters['resolution']))[m
[32m+[m
[32m+[m[32m            plt.show();[m
[32m+[m
         return vecs[m
[1mdiff --git a/structopt/tools/lammps.py b/structopt/tools/lammps.py[m
[1mindex 198fb59..b2d0f61 100644[m
[1m--- a/structopt/tools/lammps.py[m
[1m+++ b/structopt/tools/lammps.py[m
[36m@@ -372,7 +372,8 @@[m [mclass LAMMPS(object):[m
         zhilo = (hi[2] - lo[2])[m
 [m
         cell = [[xhilo,0,0],[xy,yhilo,0],[xz,yz,zhilo]][m
[31m-        self.atoms.set_cell(cell)[m
[32m+[m[32m        if all(atoms.get_pbc()):[m
[32m+[m[32m            self.atoms.set_cell(cell)[m
                 [m
         return[m
         [m
