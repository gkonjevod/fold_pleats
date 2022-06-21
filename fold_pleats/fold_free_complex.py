# Maybe I need only the original complex and the node equivalence relation.
#
# Then:
# 1. Mesh is built from the FlatFoldedGrid
# 2. Compute the free complex (the node equivalence relation)
# 3. Build the triangulated mesh.
# 4. ``Decorate'' the functions that access and set node force and mass:
#  ---each access to node.f is passed to node_equiv_class.f
#  ---each access to node.mass is passed to node_equiv_class.mass
# 5. Run the solver

# 20171125: figured out how to modify MeshNode methods to reroute any
# accesses to Node.add_force, Node.force and Node.mass to the representative of
# the node's equivalence class
# 20171126: typed out the code using UnionFind to determine the equivalence
# classes of nodes.
# 20171126: first version seems to produce reasonable output (2-level simple bowl unlocked).

from fold_geometry import vdir, dot
from fold_pleats import *

from PADS.UnionFind import UnionFind

class FreeComplex(object):
    def __init__(self, fg):
        # fg is a FlatFoldedGrid
        # This partitions the set of nodes located at pos
        # into equivalence classes according to whether they can move
        # independently.
        # Two nodes n1, n2 are in the same class if either
        # n1 traps n2 or n2 traps n1.
        # To test if n1 traps n2, it is enough to find 
        #  two neighbors m1, m2 of n1 at pos + 1*w, such that
        #  angle between n1->n2 and n1->m1 is < Pi.
        nodes = fg.grid.nodes

        self.eq_class = UnionFind()
        for n in nodes:
            tmp = self.eq_class[n]

    def partition_grid(self, fg):
        # consider separately subsets of nodes at different pos's
        poss = {}
        nodes = fg.grid.nodes
        for n in nodes:
            p = fg.pos[n]
            if p not in poss:
                poss[p] = []
            poss[p].append(n)
        for pos, nn in poss.items():
            # partition nn into subsets of equivalent nodes
            self.partition_single_pos(pos, nn, fg)

    def partition_single_pos(self, pos, nodes_here, fg):
        # partition further according to z offset
        n_z_map = {n: fg.offset[n][2] for n in nodes_here}
        z_n_map = {}
        for n in nodes_here:
            z = n_z_map[n]
            if z not in z_n_map:
                z_n_map[z] = []
            z_n_map[z].append(n)
        # now find candidate equivalences
        for z, nn in z_n_map.items():
            for n1 in nn:
                for n2 in nn:
                    if n1 == n2:
                        continue
                    if self.traps(n1, n2, fg):
                        #print(str(n1) + ' traps ' + str(n2))
                        self.eq_class.union(n1, n2)

    def traps(self, n, n_p, fg):
        # returns True if n ``traps'' n_p:
        # there are edges n--m1 and n--m2
        # such that pos[m1]==pos[m2] 
        # and dot(n->n_p, n->n_m1) > 0
        # (since pos[n] == pos[n_p], this
        # is equivalent to checking that
        # n_p is on the ``inside'' of the
        # path m1--n--m2)
        #
        # check all neighbor candidates
        n_x, n_y = n
        # horizontal first
        if n_x > 0 and n_x < fg.grid.gx:
            m1 = (n_x - 1, n_y)
            m2 = (n_x + 1, n_y)
            if fg.pos[m1] == fg.pos[m2]:
                n_3d, n_p_3d, m13d = [fg.get_3d_coords(node) for node in [n, n_p, m1]]
                #if dot(vdir(n_3d, n_p_3d), vdir(n_3d, m13d)) > 1e-3:
                if (n_3d[0] - n_p_3d[0]) * (n_3d[0] - m13d[0]) > 0:
                    # print(str(n), str(n_p), str(m1))
                    return True
        # now vertical
        if n_y > 0 and n_y < fg.grid.gy:
            m1 = (n_x, n_y - 1)
            m2 = (n_x, n_y + 1)
            if fg.pos[m1] == fg.pos[m2]:
                n_3d, n_p_3d, m13d = [fg.get_3d_coords(node) for node in [n, n_p, m1]]
                #if dot(vdir(n_3d, n_p_3d), vdir(n_3d, m13d)) > 1e-3:
                if (n_3d[1] - n_p_3d[1]) * (n_3d[1] - m13d[1]) > 0:
                    # print(str(n), str(n_p), str(m1))
                    return True
        return False
                            
    def get_class(self, n):
        return self.eq_class[n]
