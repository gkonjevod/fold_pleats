import os
from collections import OrderedDict
from functools import cmp_to_key
import numpy as np
import math
import json

from fold_geometry import *

#cellsize = 1
EPS = 1e-2
#EPS = 0.1

# 20171014-- initial github upload
# 20170926--
#   Two things only:
# 1. Given a list of pleats, make the folds:
# map the squares of the unfolded sheet
# to axis-aligned squares in the z=0 plane
# 2. Shift the vertices of the folded mesh slightly 
# to place them in 3-d folded space.
# Their z coordinates should be consistent with layer ordering.
# Their (x, y) positions should be consistent with layer ordering.
# (That is, it should be possible to take a flat sheet and arrange it so
# that its vertex-corresponding points are exactly at these (x, y, z)
# positions, without any self-intersection, but with some small amount
# of stretching.
#
# simplification:
#  ---no parsing the input pleat list (define each pleat directly as a python 
#     object)
#  ---hardcode the basic sequence
#  ---hardcode the vertices before folding begins
#
# before we have a mesh folded in 3d, we need some more basic objects



class Node(object):
    # a node is just a thing with a name; I use its grid coordinates for its name
    def __init__(self, id):
        self.id = id

    def get_id(self):
        return self.id

class Edge(object):
    # an edge joins two nodes
    def __init__(self, nodeA, nodeB):
        self.nodeA = nodeA
        self.nodeB = nodeB
        self.faces = (None, None) # can't be set without knowing the whole grid

    def get_nodes(self):
        return (self.nodeA, self.nodeB)

    def set_faces(self, ff):
        self.faces = ff[:]

    def get_faces(self):
        return self.faces
       
class Face(object):
    # a face is bounded by edges that link the face's boundary nodes
    # the order of nodes in the list here is important
    def __init__(self, nodes):
        self.nodes = nodes[:]

    def get_nodes(self):
        return self.nodes

class Grid(object):
    # my grids are rectangular, gx by gy unit rectangles.
    def __init__(self, gx, gy):
        self.gx, self.gy = gx, gy
        # make nodes first
        self.nodes = OrderedDict()
        for i in range(gx + 1):
            for j in range(gy + 1):
                self.nodes[(i, j)] = Node((i, j))
        # then make edges
        self.edges = OrderedDict()
        for i in range(gx):
            for j in range(gy):
                node00 = (i, j)
                node01 = (i, j+1)
                node10 = (i+1, j)
                # node11 = (i+1, j+1)
                self.edges[(node00, node01)] = Edge(node00, node01)
                self.edges[(node00, node10)] = Edge(node00, node10)
        for i in range(gx):
            node0 = (i, gy)
            node1 = (i+1, gy)
            self.edges[(node0, node1)] = Edge(node0, node1)
        for i in range(gy):
            node0 = (gx, i)
            node1 = (gx, i+1)
            self.edges[(node0, node1)] = Edge(node0, node1)
        # finally, make faces
        self.faces = OrderedDict()
        for i in range(gx):
            for j in range(gy):
                nodeA = (i, j)
                nodeB = (i+1, j)
                nodeC = (i+1, j+1)
                nodeD = (i, j+1)
                self.faces[(nodeA, nodeB, nodeC, nodeD)] = Face([nodeA, nodeB, nodeC, nodeD])
        # the faces incident on an edge are fixed and known
        self.compute_edges_faces()

    def compute_edges_faces(self):
        edges_faces = {e: [None, None] for e in self.edges}
        for f in self.faces:
            for (v1, v2) in zip(f, f[1:] + f[0:1]):
                if (v1, v2) in edges_faces:
                    # if f is on the right side of v1 -> v2, it's the edge's f1
                    edges_faces[(v1, v2)][1] = f
                else:
                    edges_faces[(v2, v1)][0] = f
        for e_ind, e in self.edges.items():
            e.set_faces(edges_faces[e_ind])
                    
        
# a pleat is defined by the coordinates of its mountain and valley folds,
# together with the choice of direction ('h' or 'v')
class Pleat(object):
    def __init__(self, dir, xm, xv):
        self.dir = dir
        self.xm = xm
        self.xv = xv


# to think about folding a grid, we assign space coordinates to its nodes.
# in this case, we are working with flat-folded grids and folds along 
# axis-parallel lines, so for x and y, we use standard integer values;
# z is by definition 0.
# However, due to paper thickness, we may need to compute an offset for
# each node touching other nodes. The offset is a vector in 3d, with units
# on the order of paper thickness
class FlatFoldedGrid(object):
    def __init__(self, grid):
        # keep the grid for reference
        self.grid = grid
        # set node positions to their grid locations, all offsets 0
        self.pos = {}
        self.offset = {}
        self.zmin = 0
        self.zmax = 0
        for x, y in self.grid.nodes:
            self.pos[(x, y)] = (x, y)
            self.offset[(x, y)] = (0, 0, 0)
        # define offset units
        # x and y should be a few orders of magnitude below 1,
        # and z should be another few orders of magnitude smaller
        self.unit_offx = self.unit_offy = EPS
        self.unit_offz = EPS*EPS
        # Edge assignment: V (valley), M (mountain), B (boundary), F (flat)
        # We initialize an unfolded (flat) grid
        self.init_edge_assignment()

    def get_3d_coords(self, n):
        pos_2d = self.pos[n]
        offset = self.offset[n]
        x = pos_2d[0] + offset[0] * self.unit_offx
        y = pos_2d[1] + offset[1] * self.unit_offy
        z = offset[2] * self.unit_offz
        return (x, y, z)

    def set_edge_assignment(self, e, val):
        self.edge_assignment[e] = val 

    def init_edge_assignment(self):
        # at initialization, all edges are Flat unless they are Boundary
        self.edge_assignment = {}
        for e in self.grid.edges.keys():
            if None in self.grid.edges[e].get_faces():
                self.edge_assignment[e] = 'B'
            else:
                self.edge_assignment[e] = 'F'    
    
    # now, we need ways to fold the grid.
    # let's start with a horizontal valley fold.
    # this fold runs along y = yf,involves layers at zmin and higher offset,
    # and leaves stationary nodes with y < yf
    def horiz_valley(self, yf, zmin, side):
        # side=1: y>yf stationary; side=-1: y<yf stationary
       
        # then move the nodes
        newzmax = self.zmax
        for n in self.grid.nodes:
            x, y = self.pos[n]
            if (y - yf) * side > 0:
                continue
            off = self.offset[n]
            if off[2] < zmin:
                continue
            if (y - yf) * side < 0:
                newy = reflect(yf, y)
                newoffy = reflect(0, off[1])
                newoffz = reflect(self.zmax + 1, off[2])
                self.pos[n] = (x, newy)
                self.offset[n] = (off[0], newoffy, newoffz)
                if newoffz > newzmax:
                    newzmax = newoffz
            else: # y == yf
                newoffz = self.zmax + 1
                newoffy = off[1] - side * (self.zmax - off[2])
                self.offset[n] = (off[0], newoffy, newoffz)
        oldzmax = self.zmax + 1
        if newzmax == self.zmax:
            self.zmax += 1
        else:
            self.zmax = newzmax
        return oldzmax
    
    def vert_valley(self, xf, zmin, side):
        # side=1: x>xf stationary
        newzmax = self.zmax
        for n in self.grid.nodes:
            x, y = self.pos[n]
            if (x - xf) * side > 0:
                continue
            off = self.offset[n]
            if off[2] < zmin:
                continue
            if (x - xf) * side < 0:
                newx = reflect(xf, x)
                newoffx = reflect(0, off[0])
                newoffz = reflect(self.zmax + 1, off[2])
                self.pos[n] = (newx, y)
                self.offset[n] = (newoffx, off[1], newoffz)
                if newoffz > newzmax:
                    newzmax = newoffz
            else: # x == xf
                newoffz = self.zmax + 1
                newoffx = off[0] - side * (self.zmax - off[2])
                self.offset[n] = (newoffx, off[1], newoffz)
        oldzmax = self.zmax + 1
        if newzmax == self.zmax:
            self.zmax += 1
        else:
            self.zmax = newzmax
        return oldzmax

    def is_on_creaseline(self, e, dir, x):
        if self.get_edge_assignment(e) == 'B':
            return False
        nodeA, nodeB = e
        if dir == 'h':
            if self.pos[nodeA][1] != x or self.pos[nodeB][1] != x:
                return False
            return True
        else: # dir == 'v':
            if self.pos[nodeA][0] != x or self.pos[nodeB][0] != x:
                return False
            return True

    def update_edge_assignment(self, edgelist):
        for e in edgelist:
            self.compute_edge_assignment(e)

    def zsum(self, f):
        return sum([self.offset[n][2] for n in f])
  
    def xsum(self, f):
        return sum([self.offset[n][0] for n in f])

    def ysum(self, f):
        return sum([self.offset[n][1] for n in f])

    def sign_of_f(self, f):
        p1 = self.pos[f[0]]
        p2 = self.pos[f[1]]
        p3 = self.pos[f[2]]
        n1 = (p1[0], p1[1], 0)
        n2 = (p2[0], p2[1], 0)
        n3 = (p3[0], p3[1], 0)
        return np.sign(cross(vdir(n1, n2), vdir(n2, n3))[2])
    
    def face_pos(self, f):
        return set([self.pos[n] for n in f])
    
    def compute_edge_assignment(self, e):
        #print('Computing edge assignment for (', e.nodeA, ',', e.nodeB, ')' )
        f1, f2 = e.get_faces()
        #print('Faces:', f1, f2)
        if not f1 or not f2: # if only one face, it's a boundary edge
            self.set_edge_assignment((e.nodeA, e.nodeB), 'B')
        elif self.face_pos(f1) != self.face_pos(f2): # if not in same place the edge is flat
            self.set_edge_assignment((e.nodeA, e.nodeB), 'F')
        else:
            zdiff = self.zsum(f1) - self.zsum(f2)
            if zdiff != 0:
                signprod = (self.zsum(f1) - self.zsum(f2)) * self.sign_of_f(f1)
            else: # zdiff == 0:
                #print('f1=', f1, 'and f2=', f2, 'have all equal z offsets')
                #print('f1:', [self.get_3d_coords(n) for n in f1])
                #print('f2:', [self.get_3d_coords(n) for n in f2])
                horizdiff = self.xsum(f1) + self.ysum(f1) - (self.xsum(f2) + self.ysum(f2))
                #print('horizdiff(f1 - f2)', horizdiff)
                if horizdiff > 0:
                    #print('horizontal offsets of f1 are larger so f1 is below f2')
                    signprod = -self.sign_of_f(f1)
                elif horizdiff < 0:
                    #print('horizontal offsets of f1 are smaller so f2 is below f1')
                    signprod = self.sign_of_f(f1)
                else:
                    #print('this should not happen')
                    #print('falling back to old style crease assignment')
                    f1c = [self.get_3d_coords(v) for v in f1]
                    f2c = [self.get_3d_coords(v) for v in f2]
                    return crease_sense(f1c, f2c)

                #print('sign(f1) =', self.sign_of_f(f1), 'sign(f2)=', self.sign_of_f(f2))
            if signprod > 0:
                self.set_edge_assignment((e.nodeA, e.nodeB), 'M')
            else: #signprod < 0:
                self.set_edge_assignment((e.nodeA, e.nodeB), 'V')
         
    def valley(self, dir, x, zmin, side):
        # first move the grid
        if dir == 'h':
            to_return = self.horiz_valley(x, zmin, side)
        else:
            to_return = self.vert_valley(x, zmin, side)
        # then update the edge assignment
        crease_edges = [edgeitem
                        for edgeitem in self.grid.edges.items()
                        if self.is_on_creaseline(edgeitem[0], dir, x)] 
        #print('crease edges', [item[0] for item in crease_edges])
        self.update_edge_assignment([item[1] for item in crease_edges])
        return to_return
            
    def flatten_even(self):
        # Weaker version: make z-offset values a contiguous set of integers
        zs = set()
        off = self.offset
        for n in self.grid.nodes:
            zs.add(off[n][2])
        zl = sorted(zs)
        zm = {}
        for (i, z) in enumerate(zl):
            zm[z] = i
        for n in self.grid.nodes:
            prevoff = off[n]
            off[n] = (prevoff[0], prevoff[1], zm[prevoff[2]])
        self.zmax = len(zl)-1

    def flatten_total(self):
        # Stronger version: make z-offsets contiguous at every (x, y) position
        # Not clear if this is useful. It is definitely slower.
        #
        # It is useful for initializing position in A. Ghassaei's simulator.
        # Unfortunately, this version is not correct:
        # The ties in z-offset are ignored and points that should stay at same
        # level are arbitrarily assigned a total order.
        posmap = self.pos
        offsetmap = {n: self.offset[n][2] for n in self.grid.nodes}
        positions = set(posmap.values())
        stacks = {p: [] for p in positions}
        for n, pos in posmap.items():
            stacks[pos].append((n, offsetmap[n]))
        for pos in stacks:
            stacks[pos] = sorted(stacks[pos], key = lambda x: x[1])
        for pos in stacks:
            for i, (n, z) in enumerate(stacks[pos]):
                prevoff = self.offset[n]
                self.offset[n] = (prevoff[0], prevoff[1], i)
        self.zmax = max([len(a) for a in stacks.values()])-1

    def flatten_compressed(self):
        # Stronger version: make z-offsets contiguous at every (x, y) position.
        # Preserve partial order at each (x, y).
        # This should be useful for initializing simulations.
        posmap = self.pos
        positions = set(posmap.values())
        stacks = {p: [] for p in positions}
        offsetmap = {n: self.offset[n][2] for n in self.grid.nodes}
        for n, pos in posmap.items():
            stacks[pos].append((n, offsetmap[n]))
        for pos in stacks:
            stacks[pos] = sorted(stacks[pos], key = lambda x: x[1])
        self.zmax = 0
        for pos in stacks:
            last_seen_z = -1
            curr_z = -1
            for (n, z) in stacks[pos]:
                prevoff = self.offset[n]
                if prevoff[2] > last_seen_z:
                    curr_z += 1
                last_seen_z = prevoff[2]
                self.offset[n] = (prevoff[0], prevoff[1], curr_z)
            if curr_z > self.zmax:
                self.zmax = curr_z

    def stack_faces_with_normals(self, return_map = False):
        # Build a map faces_stack that associates with each grid face
        # a position in the stacking order of faces over each
        # position grid square.
        # This version uses 3d coordinates and "is_below" function
        # from fold_geometry.
        # UNTESTED
        self.f_stack = {f: 0 for f in self.grid.faces}
        coords_map = {n: self.get_3d_coords(n) for n in self.grid.nodes}
        posmap = self.pos
        f_positions = set()
        faces_positions = {}
        for f in self.grid.faces:
            fpos = frozenset(posmap[n] for n in f)
            faces_positions[f] = fpos
            f_positions.add(fpos)
        stacks = {p: [] for p in f_positions}
        for f, p in faces_positions.items():
            stacks[p].append((f, [coords_map[n] for n in f]))
        for pos in stacks:
            stacks[pos] = sorted(stacks[pos], key = cmp_to_key(f_lt_g))
        for pos in stacks:
            for i, (f, z) in enumerate(stacks[pos]):
                self.f_stack[f] = i
        self.fmax = max([len(a) for a in stacks.values()])-1
        if return_map:
            return faces_positions

        
    def stack_faces(self, return_map = False):
        # Build a map faces_stack that associates with each grid face
        # a position in the stacking order of faces over each
        # position grid square.
        self.f_stack = {f: 0 for f in self.grid.faces}
        offsetmap = {n: self.offset[n][2] for n in self.grid.nodes}
        posmap = self.pos
        f_positions = set()
        faces_positions = {}
        for f in self.grid.faces:
            fpos = frozenset(posmap[n] for n in f)
            faces_positions[f] = fpos
            f_positions.add(fpos)
        stacks = {p: [] for p in f_positions}
        for f, p in faces_positions.items():
            stacks[p].append((f, sorted([offsetmap[n] for n in f])))
        for pos in stacks:
            stacks[pos] = sorted(stacks[pos], key = lambda x: x[1])
        for pos in stacks:
            for i, (f, z) in enumerate(stacks[pos]):
                self.f_stack[f] = i
        self.fmax = max([len(a) for a in stacks.values()])-1
        if return_map:
            return faces_positions

        
    def stack_faces_by_z_offset_only(self, return_map = False):
        # Build a map faces_stack that associates with each grid face
        # a position in the stacking order of faces over each
        # position grid square.
        # This is the original version that uses only z offsets.
        self.f_stack = {f: 0 for f in self.grid.faces}
        offsetmap = {n: self.offset[n][2] for n in self.grid.nodes}
        posmap = self.pos
        f_positions = set()
        faces_positions = {}
        for f in self.grid.faces:
            fpos = frozenset(posmap[n] for n in f)
            faces_positions[f] = fpos
            f_positions.add(fpos)
        stacks = {p: [] for p in f_positions}
        for f, p in faces_positions.items():
            stacks[p].append((f, sorted([offsetmap[n] for n in f])))
        for pos in stacks:
            stacks[pos] = sorted(stacks[pos], key = lambda x: x[1])
        for pos in stacks:
            for i, (f, z) in enumerate(stacks[pos]):
                self.f_stack[f] = i
        self.fmax = max([len(a) for a in stacks.values()])-1
        if return_map:
            return faces_positions

        
    # we now have enough to make a pleat
    def pleat(self, p, coords, compress=True):
        axis = p.dir
        if coords == 'original':
            if axis == 'v':
                xv = self.pos[(p.xv, 0)][0]
                xm = self.pos[(p.xm, 0)][0]
                #print('pleating: v, xm = %2d, xv = %2d' % (p.xm, p.xv))
            else:
                xv = self.pos[(0, p.xv)][1]
                xm = self.pos[(0, p.xm)][1]
                #print('pleating: h, xm = %2d, xv = %2d' % (p.xm, p.xv))
        elif coords == 'folded':
            xv = p.xv
            xm = p.xm
        else:
            print('FlatFoldedGrid.pleat: coordinate system (' + coords +') unknown')
        sign = np.sign(xv - xm)
        pleatzmin = self.valley(axis, xv, self.zmin, sign)
        #print('pleatzmin', pleatzmin)
        self.valley(axis, xv*2 - xm, pleatzmin, -sign)
        if compress:
            self.flatten_compressed()

    def offset_swap(self, n1, n2):
        offset1 = self.offset[n1]
        offset2 = self.offset[n2]
        self.offset[n1] = (offset2[1], offset2[0], offset2[2])
        self.offset[n2] = (offset1[1], offset1[0], offset1[2])

    def offset_paws(self, n1, n2):
        offset1 = self.offset[n1]
        offset2 = self.offset[n2]
        self.offset[n1] = (- offset2[1], - offset2[0], offset2[2])
        self.offset[n2] = (- offset1[1], - offset1[0], offset1[2])
        
    def offset_swap_avg(self, n1, n2):
        offset1 = self.offset[n1]
        offset2 = self.offset[n2]
        self.offset[n1] = ((offset2[1] + offset1[0]) / 2.0,
                           (offset2[0] + offset1[1]) / 2.0,
                           (offset2[2] + offset1[2]) / 2.0)
        self.offset[n2] = ((offset1[1] + offset2[0]) / 2.0,
                           (offset1[0] + offset2[1]) / 2.0,
                           (offset1[2] + offset2[2]) / 2.0)

        
    def xy_swap(self, n):
        offset = self.offset[n]
        self.offset[n] = (offset[1], offset[0], offset[2])

    def xy_paws(self, n):
        offset = self.offset[n]
        self.offset[n] = (- offset[1], - offset[0], offset[2])

    def layer_parity(self, here):
        # Returns True (1) if the first two faces (lower left and the one above it)
        # differ by one in their layer order. In the basic sequence case,
        # this is equivalent to the h, v order
        # NOT TO BE USED in more general situations!!!
        layer_seq = [self.f_stack[a] for a in self.get_faces_at(here)]
        return abs(layer_seq[0] - layer_seq[1]) == 1
        
    def swap_along_auxi_diag(self, x_coords, y_coords):
        self.offset_swap((x_coords[1], y_coords[0]), (x_coords[0], y_coords[1]))
        self.offset_swap((x_coords[2], y_coords[0]), (x_coords[0], y_coords[2]))
        self.offset_swap((x_coords[3], y_coords[0]), (x_coords[0], y_coords[3]))
        self.offset_swap((x_coords[2], y_coords[1]), (x_coords[1], y_coords[2]))
        self.offset_swap((x_coords[3], y_coords[1]), (x_coords[1], y_coords[3]))
        self.offset_swap((x_coords[3], y_coords[2]), (x_coords[2], y_coords[3]))
        self.xy_swap((x_coords[1], y_coords[1]))
        self.xy_swap((x_coords[2], y_coords[2]))
    def swap_along_main_diag(self, x_coords, y_coords):
        self.offset_paws((x_coords[1], y_coords[3]), (x_coords[0], y_coords[2]))
        self.offset_paws((x_coords[2], y_coords[3]), (x_coords[0], y_coords[1]))
        self.offset_paws((x_coords[0], y_coords[0]), (x_coords[3], y_coords[3]))
        self.offset_paws((x_coords[1], y_coords[1]), (x_coords[2], y_coords[2]))
        self.offset_paws((x_coords[1], y_coords[0]), (x_coords[3], y_coords[2]))
        self.offset_paws((x_coords[2], y_coords[0]), (x_coords[3], y_coords[1]))
        self.xy_paws((x_coords[1], y_coords[2]))
        self.xy_paws((x_coords[2], y_coords[1]))

    def rearrange_pleats_at(self, here):
        ff = self.get_faces_at(here)
        # rearrange layers in square [x, x+1] x [y, y+1]
        x, y = here
        x_coords = sorted(list(set([n[0] for f in ff for n in f])))
        # print(x_coords)
        y_coords = sorted(list(set([n[1] for f in ff for n in f])))
        # print(y_coords)
        lower_left_layer = self.f_stack[self.get_faces_at(here)[0]]
        if lower_left_layer in [0, 8]:
            self.swap_along_auxi_diag(x_coords, y_coords)
            print(here, 'auxi')
        else:
            self.swap_along_main_diag(x_coords, y_coords)
            print(here, 'main')
            
    def get_edge_assignment(self, e):
        return self.edge_assignment[e]

    # this should be removed; the getter function should not be computing this
    # and anyway, we should use the folding process and/or combinatorics rather
    # than 3d geometry to determine edge assignments
    # def get_edge_assignment(self, e):
    #     # TODO: this should probably use the combinatorics
    #     # and not angle geometry
    #     ff = self.grid.edges[e].get_faces()
    #     if None in ff:
    #         return 'B'
    #     else:
    #         f1, f2 = tuple(ff)
    #         f1c = [self.get_3d_coords(v) for v in f1]
    #         f2c = [self.get_3d_coords(v) for v in f2]
    #         v1, v2 = [self.get_3d_coords(v) for v in self.grid.edges[e].get_nodes()]
    #         return crease_sense(vdir(v1, v2), f1c, f2c)
        
    # again, for now do this on the fly during the folding and if necessary
    # to add a from-scratch computation, do that later
    # def compute_edge_assignment(self):
    #     for e in self.grid.edges:
    #         self.edge_assignment[e] = self.get_edge_assignment(e)

    def get_faces_at(self, f_pos):
        x, y = f_pos
        pos_set = set([(x, y), (x+1, y), (x, y+1), (x+1, y+1)])
        return [face for face in list(self.grid.faces)
                if set([self.pos[n] for n in face]) == pos_set]

    def get_ll_at(self, f_pos):
        x, y = f_pos
        pos_set = set([(x, y), (x+1, y), (x, y+1), (x+1, y+1)])
        ff = [face for face in list(self.grid.faces)
              if set([self.pos[n] for n in face]) == pos_set]
        minx = min([f[0] for face in ff for f in face])
        miny = min([f[1] for face in ff for f in face])
        return (minx, miny)
    
    def where_is_face_mapped(self, hseg, vseg):
        corners = [self.pos[hseg[0], vseg[0]],
                   self.pos[hseg[0], vseg[1]],
                   self.pos[hseg[1], vseg[0]],
                   self.pos[hseg[1], vseg[1]]]
        x = min([c[0] for c in corners])
        y = min([c[1] for c in corners])
        return (x, y)
    
    def get_four_triangles_around_edge(self, e):
        f1, f2 = e.get_faces()
        if f1 is None or f2 is None:
            return None
        e_nodes = e.get_nodes()
        nodeA = self.get_3d_coords(e_nodes[0])
        nodeB = self.get_3d_coords(e_nodes[1])
        right_nodes = list(set(f1) - set(e_nodes))
        left_nodes = list(set(f2) - set(e_nodes))
        nodeR1 = self.get_3d_coords(right_nodes[0])
        nodeR2 = self.get_3d_coords(right_nodes[1])
        nodeL1 = self.get_3d_coords(left_nodes[0])
        nodeL2 = self.get_3d_coords(left_nodes[1])
        return [(nodeA, nodeB, nodeR, nodeL) 
                for nodeR in (nodeR1, nodeR2)
                for nodeL in (nodeL1, nodeL2)]
    
# should be obsolete now    
    # next: output a FlatFoldedGrid in Fold format, so we can visualize it
    # this does not triangulate the grid
    # def make_FOLD(self, file_spec=1, 
    #                     file_creator='fold_pleats.py',
    #                     file_author='Goran Konjevod',
    #                     frame_title='Flat-folded pleats',
    #                     file_classes='singleModel',
    #                     frame_classes='foldedForm',
    #                     frame_attributes=['3D']):
    #     F = { "file_spec": file_spec,
    #           "file_creator": file_creator,
    #           "file_author": file_author,
    #           "file_classes": file_classes,
    #           "frame_title": frame_title,
    #           "frame_classes": frame_classes,
    #           "frame_attributes": frame_attributes
    #         }
        
    #     # fix a node ordering to give them FOLD-appropriate names
    #     nodes = self.grid.nodes
    #     n_ind = {n: i for (i, n) in enumerate(nodes)}
    #     # compute node coordinates
    #     node_coords = [self.get_3d_coords(n) for n in nodes]
    #     F["vertices_coords"] = node_coords
    #     # make edge list
    #     e_v = []
    #     edges = self.grid.edges
    #     # unnecessary e_ind = {e: i for (i, e) in enumerate(edges)}
    #     for e in edges:
    #         n1, n2 = edges[e].get_nodes()
    #         e_v.append([n_ind[n1], n_ind[n2]])
    #     F["edges_vertices"] = e_v
    #     # make facelist
    #     faces = self.grid.faces
    #     f_v = []
    #     for f in faces:
    #         f_v.append([n_ind[n] for n in faces[f].nodes])
    #     F["faces_vertices"] = f_v
    #     # compute edge assignments
    #     #self.compute_edge_assignment()
    #     F["edges_assignment"] = [self.edge_assignment[e] for e in self.grid.edges]
    #     # faceOrders do this later        
    #     # now the F object should be ready to dump
    #     return json.dumps(F, indent=4)

    def make_FOLD_cp(self, file_spec=1,
                     frame_title = 'fold_pleats',
                     file_creator='fold_pleats.py',
                     file_author='Goran Konjevod',
                     file_classes=['singleModel'],
                     frame_classes=['creasePattern'],
                     scale = 1,
                     edges_to_flip = [],
                     merge_nodes = False,
                     fc = None):
        if frame_classes == ['creasePattern']:
            frame_title = 'fold_pleats_CP'
            frame_attributes = ['2D']
        if frame_classes == ['foldedForm']:
            frame_title = 'fold_pleats_foldedForm'
            frame_attributes = ['3D']
        if merge_nodes:
            frame_title += '_mergednodes'
        F = { "file_spec": file_spec,
              "file_creator": file_creator,
              "file_author": file_author,
              "file_classes": file_classes,
              "frame_title": frame_title,
              "frame_classes": frame_classes,
              "frame_attributes": frame_attributes
            }
        
        # fix a node ordering to give them FOLD-appropriate names
        print('Making the crease pattern')
        print('Scale = ', scale)
        #print('  Getting nodes')
        original_nodes = self.grid.nodes
        if merge_nodes:
            node_map = fc.eq_class.parents
        else:
            node_map = {n: n for n in original_nodes}
        nodeset = set(node_map.values())
        n_ind = {n: i for (i, n) in enumerate(nodeset)}
        # compute node coordinates
        #print('  Getting node coordinates')
        if frame_classes == ['foldedForm']:
            node_coords = [self.get_3d_coords(n) for n in nodeset]
        else:
            node_coords = [[a * scale for a in n] for n in nodeset]            
        F["vertices_coords"] = node_coords
        # make edge list
        #print('  Getting edges') 
        e_v = {}
        e_a = {}
        edges = self.grid.edges
        # unnecessary e_ind = {e: i for (i, e) in enumerate(edges)}
        #print('  Making list of unique edges')
        for e in edges:
            n1, n2 = edges[e].get_nodes()
            e_ind = (n_ind[node_map[n1]], n_ind[node_map[n2]])
            e_ind2 = (e_ind[1], e_ind[0])
            if e_ind not in e_v and e_ind2 not in e_v:
                e_v[e_ind] = len(e_v)
                e_a[e_ind] = self.edge_assignment[e]
            else:
                if e_ind in e_v:
                    e_a[e_ind] = 'F'
                    e_a[e_ind2] = 'F'
        print("The CP has", len(e_v), "edges and ", len(nodeset), "nodes.")
        F["edges_vertices"] = list(e_v.values())
        # make facelist
        faces = self.grid.faces
        f_v = []
        f_sets_v = set()
        for f in faces:
            f_ind = [n_ind[node_map[n]] for n in faces[f].nodes]
            if frozenset(f_ind) not in f_sets_v:
                f_v.append(f_ind)
                f_sets_v.add(frozenset(f_ind))
        F["faces_vertices"] = f_v
        # compute edge assignments
        #self.compute_edge_assignment()
        for e in edges_to_flip:
            print('checking edge', e)
            if self.edge_assignment[e] == 'M':
                self.edge_assignment[e] = 'V'
                print('flipping edge', e, 'to', self.edge_assignment[e])
            elif self.edge_assignment[e] == 'V':
                self.edge_assignment[e] = 'M'
                print('flipping edge', e, 'to', self.edge_assignment[e])
        if merge_nodes:
            F["edges_assignment"] = list(e_a.values())
        else:
            F["edges_assignment"] = [self.edge_assignment[e] for e in self.grid.edges]
        # faceOrders do this later        
        # now the F object should be ready to dump
        return json.dumps(F, indent=4)



    def save_folded_state(self, fold_name = 'some_fold', 
                          outfilename = 'out.fold'):
        fold_string = self.make_FOLD_cp(frame_title = fold_name,
                                        frame_classes = ['foldedForm'])
        outf = open(outfilename, 'w')
        outf.write(fold_string)
        outf.close()

    def save_cp(self, fold_name = 'some_fold', 
                outfilename = 'out_cp.fold',
                scale = 1,
                edges_to_flip = []):
        fold_string = self.make_FOLD_cp(frame_title = fold_name, scale = scale,
                                        frame_classes = ['creasePattern'],
                                        edges_to_flip = edges_to_flip)
        outf = open(outfilename, 'w')
        outf.write(fold_string)
        outf.close()

    def save_folded_state_merged(self, fc,
                                 fold_name = 'some_fold', 
                                 outfilename = 'out_merged.fold'):
        fold_string = self.make_FOLD_cp(frame_title = fold_name,
                                        frame_classes = ['foldedForm'],
                                        merge_nodes = True,
                                        fc = fc)
        outf = open(outfilename, 'w')
        outf.write(fold_string)
        outf.close()
