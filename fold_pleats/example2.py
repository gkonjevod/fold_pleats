#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 00:51:13 2021

@author: goran
"""


import numpy as np

from fold_constants import constants
from fold_run_sim import Run_Simulation

from fold_simulate import Mesh, SimpleTriMesh, Solver

from make_sequence import fold_standard_sequence
from make_sequence import pyramid, in_to_out_wave, bowl, simple_tess, metaleaf
from make_sequence import random_seq, alternate_wave


def bound_z_array(a):
    minz = 0
    maxz = 24
    b = np.minimum(a[:, 2], maxz)
    a[:, 2] = np.maximum(b, minz)
    
    
if __name__ == '__main__':
    f = fold_standard_sequence(alternate_wave, 'alternate_wave_16_x16',
                               length = [16]*16, seq_mult = 1,
                               pleat_params = (3, 2, 1),
                               cp_scale = 10)

    # m = Mesh(f)
    # m.compute_free_complex()
    # m.attach_glued_nodes_by_glue_links()
    # t = SimpleTriMesh(m, 
    #                   constants['k_edge'], 
    #                   constants['d_edge'],
    #                   constants['k_crease'], 
    #                   constants['d_crease'],
    #                   constants['k_flat'], 
    #                   constants['d_flat'], 
    #                   constants['k_axial'],
    #                   constants['map_crease_targetTheta'])
    # s = Solver(t, constants['timedelta'],
    #            constraint_f = None)


    sim = Run_Simulation(f, name = 'alternate_wave_16_x16',
                         constraint_f = None,
                         constants = constants)



# from make_sequence import standard_grid_size
# from make_sequence import third_width_box_pleat_grid_size, quarter_width_box_pleat_grid_size
# from make_sequence import simple_interleaved_grid_size
# from make_sequence import skip_one_seq_wave, metaleaf
# from make_sequence import no_twist_without_rearrangements
# from make_sequence import simple_weave_without_rearrangements
# from make_sequence import confusion_without_rearrangements


# def bound_z(pos):
#     minz = 0
#     maxz = 6
#     x, y, z = pos
#     if z < minz:
#         z = minz
#     elif z > maxz: 
#         z = maxz
#     return (x, y, z)

#f = fold_standard_sequence(pyramid, 'pyramid',
#                           length = 16, seq_mult = 2,
#                           pleat_params = (3, 3, 1),
#                           cp_scale = 10)
#
#sim = Run_Simulation(f, name = 'pyramid_16_boundz_00-06', 
#                     constraint_f = bound_z_array,
#                     constants = constants)
#
#
#
#f = fold_standard_sequence(skip_one_seq_wave, 'skip_one_wave',
#                           length = 16, seq_mult = 2,
#                           pleat_params = (3, 3, 1),
#                           cp_scale = 10)
#sim = Run_Simulation(f, name = 'skip_one_wave_16', constants = constants)
#
#f = fold_standard_sequence(metaleaf, 'metaleaf',
#                           length = [5, 5, 5, 5, 5], seq_mult = 1,
#                           pleat_params = (3, 3, 1),
#                           cp_scale = 10)
#
#f = fold_standard_sequence(metaleaf, 'metaleaf',
#                           length = [6, 6, 6], seq_mult = 1,
#                           pleat_params = (3, 3, 1),
#                           cp_scale = 10)



#def flip_edges_basic_no_twist(ll = 8):
#    edges_to_flip = []
#    for i in range(ll):
#        for j in range(ll):
#            # vertical creases
#            xs = [offset + i * 16 for offset in [12, 13, 15, 16]]
#            ys = [(4 + j * 16, 5 + j * 16), (7 + j * 16, 8 + j * 16)]
#            for x in xs:
#                for y in ys:
#                    edges_to_flip.append( ((x, y[0]), (x, y[1])) )
#            # horizontal creases
#            ys = [offset + i * 16 for offset in [4, 5, 7, 8]]
#            xs = [(12 + j * 16, 13 + j * 16), (15 + j * 16, 16 + j * 16)]
#            for y in ys:
#                for x in xs:
#                    edges_to_flip.append( ((x[0], y), (x[1], y)) )
#    return edges_to_flip
#
#f = fold_standard_sequence(no_twist_without_rearrangements, 'no_twist_direct_rearr',
#                           length = ll, 
#                           compute_N = quarter_width_box_pleat_grid_size,
#                           seq_mult = 1,
#                           pleat_params = (4, 4, 1),
#                           cp_scale = 10,
#                           edges_to_flip = edges_to_flip)




#def flip_edges_simple_weave(ll = 8):
#    edges_to_flip = []
#    for i in range(ll):
#        for j in range(ll):
#            # vertical creases
#            xs = [offset + i * 6 for offset in [5, 6]]
#            ys = [(2 + j * 6, 3 + j * 6)]
#            for x in xs:
#                for y in ys:
#                    edges_to_flip.append( ((x, y[0]), (x, y[1])) )
#            # horizontal creases
#            ys = [offset + i * 6 for offset in [2, 3]]
#            xs = [(5 + j * 6, 6 + j * 6)]
#            for y in ys:
#                for x in xs:
#                    edges_to_flip.append( ((x[0], y), (x[1], y)) )
#    return edges_to_flip
#
#
#
#f = fold_standard_sequence(simple_weave_without_rearrangements, 'simple_weave_rearr',
#                           length = ll, 
#                           compute_N = simple_interleaved_grid_size,
#                           seq_mult = 1,
#                           pleat_params = (3, 2, 1),
#                           cp_scale = 10,
#                           edges_to_flip = flip_edges_simple_weave(ll))

# in origami simulator, there is a hysteresis effect between ~44% and ~62% (iirc)


# # the "1" pleats have abstract coordinates (3, 4), (10, 11), etc
# # the "-1" pleats have coordinates (6, 7), (13, 14), etc

# def flip_edges_confusion_base(ll = 8):
#     edges_to_flip = []
#     for i in range(ll):
#         for j in range(ll):
#             # vertical creases
#             xs = [offset + i * 7 for offset in [6, 7]]
#             ys = [(3 + j * 7, 4 + j * 7)]
#             for x in xs:
#                 for y in ys:
#                     edges_to_flip.append( ((x, y[0]), (x, y[1])) )
#             # horizontal creases
#             ys = [offset + i * 7 for offset in [3, 4]]
#             xs = [(6 + j * 7, 7 + j * 7)]
#             for y in ys:
#                 for x in xs:
#                     edges_to_flip.append( ((x[0], y), (x[1], y)) )
#     return edges_to_flip

# f = fold_standard_sequence(confusion_without_rearrangements, 'confusion_base',
#                            length = 6,
#                            compute_N = third_width_box_pleat_grid_size,
#                            seq_mult = 1,
#                            pleat_params = (1, 3, 1),
#                            cp_scale = 10,
#                            edges_to_flip = flip_edges_confusion_base(ll = 6))



#from fold_examples import in_to_out_wave_origcoords
#from fold_examples import metaleaf_pyramid_origcoords
#
#a = in_to_out_wave_origcoords(npleats = 8, cellsize = 1, examplespath = 'folded_examples')
#a.save_cp(fold_name = 'in_to_out_wave_08', outfilename = 'in_to_out_wave_08.fold', scale = 10 )
#sim = Run_Simulation(a, name = 'in_to_out_wave_08', constants = constants)
#
#a = in_to_out_wave_origcoords(npleats = 16, cellsize = 1, examplespath = 'folded_examples')
#a.save_cp(fold_name = 'in_to_out_wave_16', outfilename = 'in_to_out_wave_16.fold', scale = 6 )
#sim = Run_Simulation(a, name = 'in_to_out_wave_16', constants = constants)
#
#a = metaleaf_pyramid_origcoords(unit_list = [3, 5],
#                                center = 1, cellsize = 1, examplespath = 'folded_examples')
#a.save_cp(fold_name = 'metaleaf_pyramid_3-5', outfilename = 'metaleaf_pyramid_3-5.fold', scale = 1)
#sim = Run_Simulation(a, name = 'metaleaf_pyramid_3-5', constants = constants)
#
#a = in_to_out_wave_origcoords(npleats = 16, cellsize = 1, examplespath = 'folded_examples')
#a.save_cp(fold_name = 'in_to_out_wave_16', outfilename = 'in_to_out_wave_16.fold', scale = 6 )
#sim = Run_Simulation(a, name = 'in_to_out_wave_16', constants = constants)

#def standard_grid_size(length):
#    return 3 + 3 * length
#
#def fold_standard_sequence(sequence_generator, 
#                           sequence_name,
#                           length = 4, seq_mult = 1, 
#                           compute_N = standard_grid_size,
#                           pleat_params = (3, 2, 1)):
#    print('Making a ' + sequence_name + ' of length ', n)
#    N = compute_N(length * seq_mult)
#    g = Grid(N, N)
#    print('Made a grid of size', N)
#    f = FlatFoldedGrid(g)
#    print('made flat-folded grid')
#    abstract_pleat_sequence = sequence_generator(length)
#    scale, offset, width = pleat_params
#    for p in abstract_pleat_sequence:
#        f.pleat(p.make_pleat(scale=3, offset=2, width = 1), coords = 'original')
#    fold_name = sequence_name + '_' + str(length).zfill(3)
#    outfilename = fold_name + '.fold'
#    f.save_cp(fold_name, outfilename)
#    return f
