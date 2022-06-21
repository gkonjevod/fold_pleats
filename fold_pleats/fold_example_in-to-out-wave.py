import os

from fold_constants import timedelta
from fold_constants import k_edge, d_edge, k_crease, d_crease 
from fold_constants import k_flat, d_flat, k_axial
from fold_simulate import *
from fold_pleats import *
from fold_examples import *

def step(n, start, model_name):
    k = start
    for i in range(n):
        k += 1
        a_sim.step()
        if k % 100 == 0:
            a_sim.write_stl(model_name, k // 100)
            print('Snapshot')


# define the pattern and fold it --> FlatFoldedGrid
print('Folding the pattern.')
a = in_to_out_wave_origcoords(npleats = 8, cellsize = 1,
                              examplespath = '.')

print('Making a mesh.')
a_mesh = Mesh(a)

a_mesh.compute_free_complex()
# a_mesh.lock_bottom_edge()
# a_mesh.lock_left_edge()
# a_mesh.lock_top_edge()
# a_mesh.lock_right_edge()
a_mesh.update_node_references()


print('Triangulating the mesh and making springs.')
a_tri = SimpleTriMesh(a_mesh, k_edge, d_edge, 
                      k_crease, d_crease, 
                      k_flat, d_flat, k_axial)

print('Simulating the relaxation.')
a_sim = vectorized_Solver(a_tri, timedelta)

step(50000, 0, 'in_to_out_wave_08')
      
