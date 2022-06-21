#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 00:51:13 2021

@author: goran
"""
import os
from datetime import datetime

from fold_simulate import Mesh, SimpleTriMesh, Solver, vectorized_Solver

class Run_Simulation(object):
    def __init__(self, folded_state, name = 'unnamed_simulation', constraint_f = None, constants = {}):
        self.start_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        self.name = name + '_' + self.start_time
        print('Initializing ' + self.name)

        self.f = folded_state              

        print('Making a mesh')
        self.define_constants(constants)      
        self.mesh = Mesh(self.f)
        self.mesh.compute_free_complex()
        self.mesh.lock_bottom_edge()
        self.mesh.lock_top_edge()
        self.mesh.lock_left_edge()
        self.mesh.lock_right_edge()

        #self.mesh.update_node_references()
        self.mesh.attach_glued_nodes_by_glue_links()
        
        print('Triangulating the mesh and making springs.')
        self.tri_mesh = SimpleTriMesh(self.mesh, 
                                      self.constants['k_edge'], 
                                      self.constants['d_edge'],
                                      self.constants['k_crease'], 
                                      self.constants['d_crease'],
                                      self.constants['k_flat'], 
                                      self.constants['d_flat'], 
                                      self.constants['k_axial'],
                                      self.constants['map_crease_targetTheta'],
                                      use_glue_links = True)
        print('Simulating the relaxation')
        self.sim = vectorized_Solver(self.tri_mesh, self.constants['timedelta'],
        #self.sim = Solver(self.tri_mesh, self.constants['timedelta'],
                          constraint_f = constraint_f)
        self.step(0)
        print('Done: ', self.constants['simulate_steps'], 'iterations.')
        print('Current time: ', datetime.now().strftime('%Y-%m-%d %H:%M:%S'))

    def define_constants(self, constants):
        self.constants = constants # a dictionary?
        
    def step(self, start):
        snapshot_freq = self.constants['snapshot_freq']
        k = start
        for i in range(self.constants['simulate_steps']):
            k += 1
            self.sim.step()
            (snap, rem) = divmod(k, snapshot_freq)
            if rem == 0:
                #a_sim.print_springs(edge = False, face = False)
                self.sim.write_ply(self.name, snap)
                print('Snapshot', snap, datetime.now())





