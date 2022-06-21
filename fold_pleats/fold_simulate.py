import types
from math import dist, sqrt, hypot, atan2, pi
from itertools import combinations
import numpy as np
from scipy.sparse import dok_array, eye

from fold_constants import targetThetaV, targetThetaM
from fold_geometry import cross, normal_np, dot, angle
from fold_geometry import np_len, np_unit, np_unit_normal, np_uangle

from fold_pleats import Node, Edge, Face
from fold_free_complex import FreeComplex

TWO_PI = pi * 2

map_face_edge_index = { frozenset([0, 1]): 0,
                        frozenset([1, 2]): 1,
                        frozenset([2, 3]): 2,
                        frozenset([3, 0]): 3 }

class MeshNode(Node):
    # a MeshNode is a Node in a Grid, but it also
    # belongs to a Mesh moving in space, so it knows its
    # position information (original and last two),
    # current velocity and mass
    # it also keeps track of the force acting on it
    def __init__(self, node, pos, mass = 1.0): # node is a Node from a grid
        self.id = node
        self.origPos = pos
        self.lastLastPos = pos
        self.lastPos = pos
        self.lastVel = np.zeros(3) #(0.0, 0.0, 0.0)
        self.mass = mass
        self.nullforce = (0.0, 0.0, 0.0)
        self.clear_force()
        self.pinned = False
        self.constrained = False
        
    def __str__(self):
        return str(self.id)
        
    def get_lastPos(self):
        return self.lastPos
    def get_origPos(self):
        return self.origPos
    def get_lastVel(self):
        return self.lastVel
    def clear_force(self):
        # self.f = self.nullforce
        self.f = np.array(self.nullforce)
    def add_force(self, f):
        self.f = self.f + f
    def get_fm_for_move(self):
        return self.f, self.mass
    
    def constrain(self, f_adjust):
        self.constrain_f = f_adjust
        self.constrained = True

    def pin_at(self, position):
        self.origPos = self.lastLastPos = self.lastPos = position
        self.lastVel = np.zeros(3) #(0.0, 0.0, 0.0)
        self.pinned = True
    def activate_gravity(self, g):
        self.nullforce = (0.0, 0.0, -g * self.mass)

    def move(self, timedelta):
        if self.pinned:
            return
        f, m = self.get_fm_for_move()
        nextPos = f * (timedelta * timedelta / m) + 2.0 * self.lastPos - self.lastLastPos
        self.lastLastPos = self.lastPos
        self.lastPos = nextPos
        if self.constrained:
            self.lastPos = self.constrain_f(self.lastPos)
            print('Constraint active: last pos', nextPos, 'adjusted to', self.lastPos)
        self.lastVel = (nextPos - self.lastLastPos) / timedelta
        # s = 'N ' + str(self) + ' from (%5.3f, %5.3f, %5.3f)' % self.lastLastPos
        # s += ' to (%5.3f, %5.3f, %5.3f)' % self.lastPos
        # s += ' under f = (%5.3f, %5.3f, %5.3f)' % self.f
        # print(s) 

    def make_node_refer_to_another(self, n_ref, mass):
        n_ref.mass = mass
        def add_ref_force(self, f):
            n_ref.add_force(f)
        self.add_force = types.MethodType(add_ref_force, self)
        def get_fm_for_move(self):
            return n_ref.f, n_ref.mass
        self.get_fm_for_move = types.MethodType(get_fm_for_move, self)
        
        
class MeshEdge(Edge):
    # like a MeshNode, MeshEdge is an Edge in a Grid moving in space
    def __init__(self, nodeA, nodeB, face1, face2):
        self.nodeA, self.nodeB = nodeA, nodeB
        self.face1, self.face2 = face1, face2
        origPosA = nodeA.get_origPos()
        origPosB = nodeB.get_origPos()
        self.nomLength = dist(origPosA, origPosB)
        self.assignment = 'F'
    def get_faces(self):
        return self.face1, self.face2
    def get_nomLength(self):
        return self.nomLength
    def set_assignment(self, a):
        self.assignment = a
    def get_assignment(self):
        return self.assignment

class MeshFace(Face):
    # ...
    def __init__(self, nodes): # these are MeshNodes
        self.nodes = nodes[:]
        self.update_geometry()

    def update_geometry(self):
        self.normal = normal_np(self.nodes[0].get_lastPos(),
                                self.nodes[1].get_lastPos(),
                                self.nodes[2].get_lastPos())
                                #[a.get_lastPos() for a in self.nodes[:3]])
    def get_normal(self):
        return self.normal

class Mesh(object):
    # a 2-complex (2-dimensional polyhedral complex) in R^3
    # --mesh nodes
    # --mesh edges
    # --mesh faces
    # Used as an intermediate step between a folded grid and a triangulated mesh
    # for the simulation.
    # Enables the computation and use of the free complex, edge gluing
    # and mesh draping.
    def __init__(self, folded_grid):
        self.fg = folded_grid
        # self.fg.compute_edge_assignment()
        self.edge_assignment = self.fg.edge_assignment

        g = self.fg.grid
        
        # first, copy over the nodes
        self.nodes = []
        self.node_map = {}
        for v in g.nodes:
            new_v = MeshNode(v, np.array(self.fg.get_3d_coords(v)))
            self.nodes.append(new_v)
            self.node_map[v] = new_v
            
        # then, copy over the faces
        self.faces = []
        face_map = {}
        face_map[None] = None
        for f_ind, f in g.faces.items():
            newNodes = [self.node_map[n] for n in f.get_nodes()]
            new_f = MeshFace(newNodes)
            self.faces.append(new_f)
            face_map[f_ind] = new_f

        # finally, copy over the edges
        self.edges = []
        for e_ind, e in g.edges.items():
            newA, newB = [self.node_map[n] for n in e.get_nodes()]
            new_f1, new_f2 = [face_map[f] for f in e.get_faces()]
            new_e = MeshEdge(newA, newB, new_f1, new_f2)
            new_e.set_assignment(self.edge_assignment[e_ind])
            self.edges.append(new_e)
            

    def set_constraint(self, f_adjust):
        for n in self.nodes:
            n.constrain(f_adjust)

    def pin(self, n, position = None):
        if not position:
            # pin n at n's current position
            position = np.array(self.fg.get_3d_coords(n))
            # position = (position[0], position[1], 0)
        self.node_map[n].pin_at(position)

    def activate_gravity(self):
        g = 1.0 / len(self.nodes)
        for n in self.nodes:
            n.activate_gravity(g)
        
    def compute_free_complex(self):
        self.fc = FreeComplex(self.fg)
        self.fc.partition_grid(self.fg)
        
    def update_node_references(self):
        # modify the behavior of nodes
        # to allow simulation code reuse with glued nodes
        mass = self.fc.eq_class.calculate_class_sizes()
        # print(mass)
        for n in self.nodes:
            n_ref_id = self.fc.eq_class[n.get_id()]
            n_ref = self.node_map[n_ref_id]
            if n == n_ref:
                continue
            n.make_node_refer_to_another(n_ref, mass[n_ref_id])
            
    def compute_node_eq_classes(self):
        classes = self.eq_classes = {}
        for n in self.fg.grid.nodes:
            n_ref = self.fc.eq_class[n]
            if n_ref not in classes:
                classes[n_ref] = set([n, n_ref])
            else:
                classes[n_ref].add(n)
            
    def attach_glued_nodes_by_glue_links(self):
        self.glue_links = []
        self.compute_node_eq_classes()
        for eq_class in self.eq_classes.values():
            if len(eq_class) < 2:
                continue
            for n1, n2 in combinations(eq_class, 2):
                self.glue_links.append((self.node_map[n1], self.node_map[n2]))
            
    def glue_through_at_pos_use_fc(self, pos):
        # glue together all nodes at pos.
        # Uses the free complex, so it should be
        # done AFTER the free complex calculation
        # and BEFORE node reference updates.
        nn = [n for n in self.fg.grid.nodes if self.fg.pos[n] == pos]
        if len(nn) > 1:
            #print('Gluing ' + str(len(nn)) + ' nodes.')
            n0 = nn[0]
            for n in nn[1:]:
                self.fc.eq_class.union(n, n0)
                #print('Gluing', n, 'and', n0)
                
    def lock_bottom_edge(self):
        fg = self.fg
        zerox, zeroy = fg.pos[(0, 0)]
        zerox = int(zerox)
        zeroy = int(zeroy)
        gx, gy = fg.grid.gx, fg.grid.gy
        print('zerox ' + str(zerox) + ' gx ' + str(gx))
        for i in range(zerox, gx + 1):
            self.glue_through_at_pos_use_fc((i, zeroy))
            self.glue_through_at_pos_use_fc((i, zeroy + 1))
        print('Bottom edge locked.')

    def lock_left_edge(self):
        fg = self.fg
        zerox, zeroy = fg.pos[(0, 0)]
        zerox = int(zerox)
        zeroy = int(zeroy)
        gx, gy = fg.grid.gx, fg.grid.gy 
        for i in range(zeroy, gy + 1):
            self.glue_through_at_pos_use_fc((zerox, i))
            self.glue_through_at_pos_use_fc((zerox + 1, i))
        print('Left edge locked.')

    def lock_top_edge(self):
        fg = self.fg
        zerox, zeroy = fg.pos[(0, 0)]
        zerox = int(zerox)
        zeroy = int(zeroy)
        gx, gy = fg.grid.gx, fg.grid.gy
        maxx, maxy = fg.pos[(gx, gy)]
        for i in range(zerox, maxx + 1):
            self.glue_through_at_pos_use_fc((i, maxy))
            self.glue_through_at_pos_use_fc((i, maxy - 1))
        print('Top edge locked.')

    def lock_right_edge(self):
        fg = self.fg
        zerox, zeroy = fg.pos[(0, 0)]
        zerox = int(zerox)
        zeroy = int(zeroy)
        gx, gy = fg.grid.gx, fg.grid.gy
        maxx, maxy = fg.pos[(gx, gy)]
        for i in range(zeroy, maxy + 1):
            self.glue_through_at_pos_use_fc((maxx, i))
            self.glue_through_at_pos_use_fc((maxx - 1, i))
        print('Right edge locked.')

    # 20171126
    # The first attempt at simulating the simple bowl
    # curves differently from what I usually make by folding.
    # It would be useful to have a way to distort the mesh slightly
    # so I can give it a different starting shape.
    # The argument surface is a function of two variables that
    # behaves reasonably on the unit square.  The mesh coordinates
    # get scaled to the unit square before draping.
    # There are no safety checks.
    def drape_over(self, surface):
        minx = min([n.get_origPos()[0] for n in self.nodes])
        miny = min([n.get_origPos()[1] for n in self.nodes])
        maxx = max([n.get_origPos()[0] for n in self.nodes])
        maxy = max([n.get_origPos()[1] for n in self.nodes])

        tr_x = lambda x: 2 * (0.0 + x - minx) / (maxx - minx) - 1
        tr_y = lambda y: 2 * (0.0 + y - miny) / (maxy - miny) - 1
        for n in self.nodes:
            (x, y, z) = n.get_origPos()
            # print('(tr_x, tr_y) = (%5.3f, %5.3f)' % (tr_x(x), tr_y(y)))
            pos = (x, y, - 10 * surface(tr_x(x), tr_y(y)) + z)
            # print('(%5.3f, %5.3f, %5.3f) moved to (%5.3f, %5.3f, %5.3f)' % (x, y, z, pos[0], pos[1], pos[2]))
            n.origPos = n.lastPos = n.lastLastPos = pos
        
class EdgeSpring(object):
    def __init__(self, nodeA, nodeB, k = 1.0, d = 1.0):
        self.nodeA, self.nodeB = nodeA, nodeB
        self.nomLength = dist(nodeA.get_origPos(), nodeB.get_origPos())
        self.curLength = dist(nodeA.get_lastPos(), nodeB.get_lastPos())
        self.k = k
        self.d = d
        self.f = np.zeros(3) #(0.0, 0.0, 0.0)
        
    def __str__(self):
        s = str(self.nodeA) + '---' + str(self.nodeB)
        s += ' nomL = %5.3f, curL = %5.3f' % (self.nomLength, self.curLength)
        s += ' from (%5.3f %5.3f, %5.3f) to' % self.nodeA.get_lastPos()
        s += ' (%5.3f %5.3f, %5.3f) ' % self.nodeB.get_lastPos()
        return s

    def get_nodes(self):
        return self.nodeA, self.nodeB
    
    def update_geometry(self):
        v, w = self.nodeA.get_lastPos(), self.nodeB.get_lastPos()
        vw = w - v
        self.curLength = ll = hypot(*vw)         
        stretch = ll - self.nomLength
        f = vw * self.k * stretch / ll
        self.f = f + self.d * (self.nodeB.get_lastVel() - self.nodeA.get_lastVel())
                
    def get_force(self, v):
        if v == self.nodeA:
            # print('Edge f (%5.3f %5.3f, %5.3f) on ' % self.f + str(v))
            return self.f
        else: 
            # print('Edge f (%5.3f %5.3f, %5.3f) on ' % scale(self.f, -1) + str(v))
            return -self.f #scale(self.f, -1.0)



class Crease(object):
    def __init__(self, e, k, d, targetTheta):
        self.length = e.get_nomLength()
        self.nodeA, self.nodeB = e.get_nodes() # correspond to nodes 3 and 4 in AG
        self.face1, self.face2 = e.get_faces()
        # cleanup
        self.nodeC = list(set(self.face1.get_nodes()).difference(set([self.nodeA, self.nodeB])))[0]
        self.nodeD = list(set(self.face2.get_nodes()).difference(set([self.nodeA, self.nodeB])))[0]
        self.k = k
        self.d = d
        self.targetTheta = targetTheta
        self.currentTheta = targetTheta
        self.angF = 0.0
        self.update_geometry()

    def __str__(self):
        s = str(self.nodeA) + '---' + str(self.nodeB)
        s += ' n1 = (%5.3f, %5.3f, %5.3f)' % self.n1
        s += ' n2 = (%5.3f, %5.3f, %5.3f)' % self.n2        
        s += ' Th = %5.3f, target = %5.3f' % (self.currentTheta, self.targetTheta)
        s += ' aF = %5.3f' % self.angF
        return s
        
    def get_nodes(self):
        return self.nodeA, self.nodeB, self.nodeC, self.nodeD
    
    def update_geometry(self):
        posA = self.nodeA.get_lastPos()
        posB = self.nodeB.get_lastPos()
        posC = self.nodeC.get_lastPos()
        posD = self.nodeD.get_lastPos()
        cv = posB - posA
        cL = hypot(*cv)
        v0 = cv / cL
        v1 = posC - posA
        v2 = posD - posA
        proj1L = np.dot(v0, v1)
        proj2L = np.dot(v0, v2)
        self.h1 = sqrt(np.dot(v1, v1) - proj1L * proj1L)
        self.h2 = sqrt(np.dot(v2, v2) - proj2L * proj2L)
        self.coef1 = proj1L / cL
        self.coef2 = proj2L / cL
        self.n1 = self.face1.get_normal()
        self.n2 = self.face2.get_normal()
        theta = atan2(np.dot(cross(cv, self.n1), self.n2), np.dot(self.n1, self.n2))
        diff = theta - self.currentTheta
        origDiff = diff
        if diff < -5:
            diff += TWO_PI
        elif diff > 5:
            diff -= TWO_PI
        self.currentTheta += diff
        # print(self)
        self.angF = self.k * (self.targetTheta - self.currentTheta)
        
    def get_force(self, v):
        if v == self.nodeA:
            c1 = 1.0 - self.coef1
            c2 = 1.0 - self.coef2
            term1 = self.n1 * (c1/self.h1)
            term2 = self.n2 * (c2/self.h2)
            f = -self.angF * (term2 + term1)
        elif v == self.nodeB:
            c1 = self.coef1
            c2 = self.coef2
            term1 = self.n1 * (c1/self.h1)
            term2 = self.n2 * (c2/self.h2)
            f = -self.angF * (term2 + term1)
        else:
            if v == self.nodeC:
                n = self.n1
                momentArm = self.h1
            else:
                n = self.n2
                momentArm = self.h2
            f = (self.angF / momentArm) * n
        # print('Crease force (%5.3f %5.3f, %5.3f) on ' % f + str(v))
        return f

class FaceSpring(object):
    def __init__(self, nodes, axialStiffness):
        self.axialStiffness = axialStiffness
        self.tol = 0.0000001
        self.collapsed = False
        
        self.nodeA, self.nodeB, self.nodeC = nodes
        self.posA = self.nodeA.get_origPos()
        self.posB = self.nodeB.get_origPos()
        self.posC = self.nodeC.get_origPos()
        self.nomAngles = self.get_angles()
        self.update_geometry()
        
    def __str__(self):
        s = str(self.nodeA) + '--' + str(self.nodeB) + '--' + str(self.nodeC)
        # s += ' anglesDiff = (%5.3f, %5.3f, %5.3f)' % self.anglesDiff
        return s
        
    def get_nodes(self):
        return self.nodeA, self.nodeB, self.nodeC

    def get_normal(self):
        return self.normal

    def get_angles(self):
        angA = angle(self.posB - self.posA, self.posC - self.posA)
        angB = angle(self.posC - self.posB, self.posA - self.posB)
        angC = angle(self.posA - self.posC, self.posB - self.posC)
        return np.array([angA, angB, angC])
        
    def update_geometry(self):
        self.posA = self.nodeA.get_lastPos()
        self.posB = self.nodeB.get_lastPos()
        self.posC = self.nodeC.get_lastPos()
        self.normal = normal_np(self.posA, self.posB, self.posC)

        AB = self.posB - self.posA
        BC = self.posC - self.posB
        CA = self.posA - self.posC
        len_AB = hypot(*AB)
        len_BC = hypot(*BC)
        len_CA = hypot(*CA)

        if len_AB < self.tol or len_BC < self.tol or len_CA < self.tol:
            self.collapsed = True
            return
        else:
            self.collapsed = False
            
        self.ab = AB / len_AB
        self.bc = BC / len_BC
        self.ca = CA / len_CA
        self.normCrossAB = cross(self.normal, self.ab) / len_AB
        self.normCrossBC = cross(self.normal, self.bc) / len_BC
        self.normCrossCA = cross(self.normal, self.ca) / len_CA

        self.anglesDiff = self.nomAngles - self.get_angles()
        self.anglesDiff = self.anglesDiff * self.axialStiffness/100.0
        
    def get_force(self, v): # from AG
        if self.collapsed:
            return (0.0, 0.0, 0.0)
        if v == self.nodeA:
            f = (self.normCrossCA + self.normCrossAB) * self.anglesDiff[0] - \
                self.normCrossAB * self.anglesDiff[1] - \
                self.normCrossCA * self.anglesDiff[2]
        elif v == self.nodeB: 
            f = (self.normCrossAB + self.normCrossBC) * self.anglesDiff[1] - \
                self.normCrossBC * self.anglesDiff[2] - \
                self.normCrossAB * self.anglesDiff[0]   
        else: # v == self.nodeC
            f = (self.normCrossCA + self.normCrossBC) * self.anglesDiff[2] - \
                self.normCrossCA * self.anglesDiff[0] - \
                self.normCrossBC * self.anglesDiff[1]
        # print('Face force (%5.3f %5.3f, %5.3f) on ' % f + str(v))
        return f


class SimpleTriMesh(Mesh):
    # Triangulated by adding a single diagonal to each square.
    # 
    def __init__(self, mesh,
                 k_edge, d_edge,
                 k_crease, d_crease,
                 k_flat, d_flat,
                 k_axial,
                 map_crease_targetTheta = {'V': pi * (1 - 1/4),
                                           'M': -pi * (1 - 1/4),
                                           'F': 0.0 },
                 use_glue_links = True):
        self.nodes = mesh.nodes # no need to copy; just use the existing set

        self.k_edge = k_edge
        self.d_edge = d_edge
        self.use_glue_links = use_glue_links
        print('Use glue links:', use_glue_links)
        
        if 'eq_classes' not in mesh.__dict__:
            mesh.compute_node_eq_classes() # as long as the nodes stay the same
        
        # make a node-to-node equivalence array to capture the
        # free complex equivalence relation

        mesh_node_map = mesh.node_map
        node_index_map = {n: i for (i, n) in enumerate(self.nodes)}
        #self.node_to_node_eq = eye(len(self.nodes), dtype = np.int8,
        #                           format = 'csr')
        if not self.use_glue_links:
            for n in mesh.nodes:
                n1 = n.get_id()
                if n1 not in mesh.eq_classes:
                    continue
                equivs = mesh.eq_classes[n1]
                for n2 in equivs:
                    if n2 == n1:
                        continue
                    self.node_to_node_eq[node_index_map[mesh_node_map[n1]],
                                         node_index_map[mesh_node_map[n2]]] == 1
                    #print(n1, 'is equivalent to', n2)
                    
        # faces will be all new (each existing face gets split into two)
        self.faces = []
        newfaces = {}

        # we'll copy the original edges later when we know their faces
        self.edges = []
        
        # for each mesh face,
        for f in mesh.faces:
            nn = f.get_nodes()

            # make two new faces
            f1nodes = [nn[0], nn[1], nn[2]]
            newf1 = MeshFace(f1nodes)
            self.faces.append(newf1)
            
            f2nodes = [nn[2], nn[3], nn[0]]
            newf2 = MeshFace(f2nodes)
            self.faces.append(newf2)
                       
            # keep an access map
            newfaces[(f, 0)] = newf1
            newfaces[(f, 1)] = newf2
            
            # finally, make the new edge
            self.edges.append(MeshEdge(nn[0], nn[2], newf1, newf2))

        # copy the old edges now 
        for e in mesh.edges:
            nodeA, nodeB = e.get_nodes()
            f1, f2 = e.get_faces()
            if f1 is None:
                f1new = None
            else:
                s = set([f1.get_nodes().index(nodeA), f1.get_nodes().index(nodeB)])
                if 1 in s:
                    f1new = newfaces[(f1, 0)]
                else:
                    f1new = newfaces[(f1, 1)]
            if f2 is None:
                f2new = None
            else:
                s = set([f2.get_nodes().index(nodeA), f2.get_nodes().index(nodeB)])
                if 1 in s:
                    f2new = newfaces[(f2, 0)]
                else:
                    f2new = newfaces[(f2, 1)]
            newEdge = MeshEdge(nodeA, nodeB, f1new, f2new)
            newEdge.set_assignment(e.get_assignment())
            self.edges.append(newEdge)
        
        # make EdgeSprings (for all edges)
        self.edge_springs = []
        for e in self.edges:
            nodeA, nodeB = e.get_nodes()
            # print('New edge spring from '+str(nodeA)+' to '+str(nodeB))
            e_s = EdgeSpring(nodeA, nodeB, k_edge, d_edge)
            # print('EdgeSpring ', str(e_s))
            self.edge_springs.append(e_s)
        
        # if there are glue links, add those too
        if self.use_glue_links:
            for n1, n2 in mesh.glue_links:
                d = dist(n1.origPos, n2.origPos)
                k = k_edge/d
                e_s = EdgeSpring(n1, n2, k, d_edge)
                #print('Added glue link between', n1.id, 'and',
                #      n2.id, 'with k=', k)
                self.edge_springs.append(e_s)
                
            
        # make Creases (for non-boundary edges)
        self.creases = []
        for e in self.edges:
            if e.get_assignment() == 'B':
                continue
            #f1, f2 = e.get_faces()
            #nodeA, nodeB = e.get_nodes()
            #if nodeA in newnodes or nodeB in newnodes:
            #    continue
            # I also need the target angle
            # this depends on the crease sense.
            sense = e.get_assignment()
            targetTheta = map_crease_targetTheta[sense]
            if sense == 'F':
                k = k_flat
                d = d_flat
            else:
                k = k_crease
                d = d_crease
            c = Crease(e, k, d, targetTheta)
            # print('Crease', str(c))
            self.creases.append(c)
            
        # make FaceSprings (for all faces)
        # right now not using these
        self.face_springs = []
        #for f in self.faces:
        #    f_s = FaceSpring(f.get_nodes(), k_axial)
        #    # print('FaceSpring', str(f_s))
        #    self.face_springs.append(f_s)


class Solver(object):
    """Simulate the relaxation of a folded model.

    trimesh is a triangulated mesh; it has nodes, edges, faces and springs.
    In each step, the solver
      1. updates geometry of springs and faces
      2. calculates external and internal forces
      3. moves points
    I need to
       ---make sure these steps can be run separately like that
       ---make sure the geometry update can be done in a generic loop
       ---make sure the data to compute all the forces is available
    Then I may be able to
       ---add code for creating a scene from trimeshes and associated info
    In any case, I need code for saving the state of the simulation
    (for individual meshes, this may be possible using an extension of the fold format;
    for scenes, I don't know)
    """
    def __init__(self, trimesh, timedelta, constraint_f = None):
        self.nodes = trimesh.nodes
        self.edges = trimesh.edges
        self.faces = trimesh.faces
        self.edge_springs = trimesh.edge_springs
        self.creases = trimesh.creases
        self.face_springs = trimesh.face_springs
        self.timedelta = timedelta
        if constraint_f:
            self.constrained = True
            self.constraint_f = constraint_f
        
    def step(self):
        self.update_geometry() 
        self.update_external_forces() 
        self.calculate_internal_forces()
        self.move_points()

    def update_geometry(self):
        # recalculate lengths, angles, normals, projections
        for e in self.edge_springs:
            e.update_geometry()
        for f in self.faces: # before creases because of face normals
            f.update_geometry()
        for c in self.creases:
            c.update_geometry()
        for f_s in self.face_springs: # empty
            f_s.update_geometry()
    
    def update_external_forces(self):
        # nothing for now;
        # this will help make changes during the simulation
        pass

    def print_springs(self, edge = True, crease = True, face = True):
        if edge:
            for e in self.edge_springs:
                print(e)
        if crease:
            for c in self.creases:
                print(c)
        if face:
            for f in self.face_springs:
                print(f)
    
    def calculate_internal_forces(self):
        self.clear_internal_forces()
        for e in self.edge_springs:
            self.add_edge_forces(e)
        for c in self.creases:
            self.add_crease_forces(c)
        for f in self.face_springs:
            self.add_face_forces(f)
            
    def clear_internal_forces(self):
        for v in self.nodes:
            v.clear_force()
            
    def add_edge_forces(self, e):
        # compute the forces exerted by e on its endpoints and
        # add them to the current force totals
        for v in e.get_nodes():
            # print('Edge spring ', end='')
            v.add_force(e.get_force(v))
        
    def add_crease_forces(self, c):
        # compute the forces exerted by c on its four vertices and
        # add them to the current force totals
        # c is constructed from e: e's endpoints, as well as
        # the other vertices of both faces e belongs to are the four
        # vertices that feel the force of c's flexing
        for v in c.get_nodes():
            # print('Crease spring ', end='')
            v.add_force(c.get_force(v))
        
    def add_face_forces(self, f):
        for v in f.get_nodes():
            v.add_force(f.get_force(v))

    def move_points(self):
        for v in self.nodes:
            v.move(self.timedelta)

    # def write_stl(self, model, i):
    #     filename = model + str(i).zfill(4) + '.stl'
    #     outf = open(filename, 'w')
    #     outf.write('solid ' + model + '\n')
    #     for f in self.faces:
    #         normal = f.get_normal()
    #         pos = [n.get_lastPos() for n in f.nodes]
    #         outf.write('  facet normal %8.6f %8.6f %8.6f\n' % tuple(normal))
    #         outf.write('    outer loop\n')
    #         outf.write('      vertex %8.6f %8.6f %8.6f\n' % tuple(pos[0]))
    #         outf.write('      vertex %8.6f %8.6f %8.6f\n' % tuple(pos[1]))
    #         outf.write('      vertex %8.6f %8.6f %8.6f\n' % tuple(pos[2]))
    #         outf.write('    endloop\n')
    #         outf.write('  endfacet\n')
    #     outf.write('endsolid\n')
    #     outf.close()
        
    def write_ply(self, model, i):
        filename = model + '_' + str(i).zfill(4) + '.ply'
        outf = open(filename, 'w')
        outf.write('ply \n')
        outf.write('format ascii 1.0\n')
        outf.write('comment ' + model + '\n')
        outf.write('element vertex ' + str(len(self.nodes)) + '\n')
        outf.write('property float x\n')
        outf.write('property float y\n')
        outf.write('property float z\n')
        outf.write('element face ' + str(len(self.faces)) + '\n')
        outf.write('property list uchar int vertex_index\n')
        outf.write('end_header\n')
        for n in range(len(self.nodes)):
            pos = self.node_lastpos[n]
            outf.write('%8.6f %8.6f %8.6f\n' % tuple(pos))
        for f in self.faces: 
            nodes = [self.node_index_map[i] for i in f.get_nodes()]
            outf.write('3 %d %d %d\n' % tuple(nodes))
        outf.close()


    def write_node_coords(self, filename):
        outf = open(filename, 'w')
        outf.write('nodex nodey x y z x0 y0 z0\n')
        for n in self.nodes:
            grid_x, grid_y = n.id
            x, y, z = n.get_lastPos()
            x0, y0, z0 = n.get_origPos()
            outf.write('%d %d %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f\n' %
                       (grid_x, grid_y, x, y, z, x0, y0, z0))
        outf.close()


# this takes too much memory for more than 32 by 32 pleats
# def make_incidence_matrix(rows, cols, incidence_list):
#     mm = np.zeros((len(rows), len(cols)))
#     for i, inc in enumerate(incidence_list):
#         mm[inc, i] = 1
#     return mm
def make_incidence_matrix(rows, cols, incidence_list):
    mm = dok_array((len(rows), len(cols)), dtype = np.int8)
    for i, inc in enumerate(incidence_list):
        mm[inc, i] = 1
    return mm.tocsr()
        
# def make_plus_minus_incidence_matrix(rows, cols, incidence_list_pos, incidence_list_neg):
#     mm = np.zeros((len(rows), len(cols)), dtype = np.int8)
#     for i, inc in enumerate(incidence_list_pos):
#         mm[inc, i] = 1
#     for i, inc in enumerate(incidence_list_neg):
#         mm[inc, i] = -1
#     return mm
def make_plus_minus_incidence_matrix(rows, cols, incidence_list_pos, incidence_list_neg):
    mm = dok_array((len(rows), len(cols)), dtype = np.int8)
    for i, inc in enumerate(incidence_list_pos):
        mm[inc, i] = 1
    for i, inc in enumerate(incidence_list_neg):
        mm[inc, i] = -1
    return mm.tocsr()




class vectorized_Solver(object):
    """Vectorize the simulation by storing the data (positions, forces) in
    numpy arrays.
    """
    def __init__(self, trimesh, timedelta, constraint_f = None):

        print("initializing vectorized_Solver")
        
        # Basic triangulated mesh data
        self.nodes = trimesh.nodes
        self.edges = trimesh.edges
        self.faces = trimesh.faces
        print("Copied nodes, edges and faces")

        self.es_k = trimesh.k_edge
        self.es_d = trimesh.d_edge
        
        ### Node position, velocity and force
        
        # initiate node positions to those in trimesh
        self.node_origpos = np.array([n.get_origPos() for n in self.nodes])
        self.node_lastpos = self.node_origpos.copy()
        
        # the simulation is just starting, so previous position same, velocity 0
        self.node_lastlastpos = self.node_lastpos.copy()
        self.node_lastvel = np.zeros((len(self.nodes), 3))

        # by default, there are no external forces (e.g. gravity)
        self.node_nullF = np.zeros((len(self.nodes), 3))

        # by default, all nodes are free to move
        self.node_pinned = np.full(len(self.nodes), False)

        print("Initialized node positions, velocities, forces.")

        # In order to create reference arrays that will help us
        # retrieve appropriate node/face information, we need mappings
        # from nodes/faces to their respective indices in node/face lists
        node_index_map = {n: i for (i, n) in enumerate(self.nodes)}
        self.node_index_map = node_index_map # I use this in a print function
        face_index_map = {f: i for (i, f) in enumerate(self.faces)}
        #self.node_to_node_eq = trimesh.node_to_node_eq
        print("Defined node and face index maps.")
        
        # initialize all nodes to unit mass, multiplied by the node-to-node equivalence
        # matrix, to distribute correct mass value to each node
        #self.node_mass = self.node_to_node_eq @ np.ones((len(self.nodes), 1))
        self.node_mass = np.ones((len(self.nodes), 1))
        print("Initialized node masses.")
        
        ### Edge structural info (node and face references)
        self.edge_nodeA = [node_index_map[e.get_nodes()[0]] for e in self.edges]
        self.edge_nodeB = [node_index_map[e.get_nodes()[1]] for e in self.edges]
 
        # Face structural info (node and edge references; normals)
        self.face_nodeA = [node_index_map[f.get_nodes()[0]] for f in self.faces]
        self.face_nodeB = [node_index_map[f.get_nodes()[1]] for f in self.faces]
        self.face_nodeC = [node_index_map[f.get_nodes()[2]] for f in self.faces]
        print("Defined edge_node, face_node maps.")
        
        # Edge springs
        self.edge_springs = trimesh.edge_springs
        print("Copied edge springs.")
        
        # we store references to edge endpoints in these two arrays
        self.es_nodeA = [node_index_map[e.get_nodes()[0]] for e in self.edge_springs]
        self.es_nodeB = [node_index_map[e.get_nodes()[1]] for e in self.edge_springs]

        self.es_node_inc = make_plus_minus_incidence_matrix(self.nodes, self.edge_springs,
                                                            self.es_nodeA, self.es_nodeB)

        # self.es_nodeA_inc = make_incidence_matrix(self.nodes, self.edge_springs, self.es_nodeA)
        # self.es_nodeB_inc = make_incidence_matrix(self.nodes, self.edge_springs, self.es_nodeB)

        # the nominal length is fixed
        self.es_nomL = np.linalg.norm(self.node_origpos[self.es_nodeA] - 
                                      self.node_origpos[self.es_nodeB], axis = 1)

        print("Initialized edge springs.")

        self.update_edge_springs()
        print("Updated edge springs.")
        
        # Crease fixed info (four nodes, two faces, constants)
        self.creases = trimesh.creases
        self.cs_nodeA = [node_index_map[e.get_nodes()[0]] for e in self.creases]
        self.cs_nodeB = [node_index_map[e.get_nodes()[1]] for e in self.creases]
        self.cs_nodeC = [node_index_map[e.get_nodes()[2]] for e in self.creases]
        self.cs_nodeD = [node_index_map[e.get_nodes()[3]] for e in self.creases]
        
        self.cs_nodeA_inc = make_incidence_matrix(self.nodes, self.creases, self.cs_nodeA)
        self.cs_nodeB_inc = make_incidence_matrix(self.nodes, self.creases, self.cs_nodeB)
        self.cs_nodeC_inc = make_incidence_matrix(self.nodes, self.creases, self.cs_nodeC)
        self.cs_nodeD_inc = make_incidence_matrix(self.nodes, self.creases, self.cs_nodeD)
        print("Defined crease maps.")
        
        self.crease_face1 = [face_index_map[c.face1] for c in self.creases]
        self.crease_face2 = [face_index_map[c.face2] for c in self.creases]
        self.crease_k = np.array(list(c.k for c in self.creases))
        self.crease_d = np.array(list(c.d for c in self.creases))
        self.crease_curTheta = [c.targetTheta for c in self.creases]
        self.crease_tgtTheta = [c.targetTheta for c in self.creases]
        print("Initialized creases.")
        
        # Face info (still used for the computation of normals)
        self.faces = trimesh.faces
        self.fs_nodeA = np.array([node_index_map[f.get_nodes()[0]] for f in self.faces])
        self.fs_nodeB = np.array([node_index_map[f.get_nodes()[1]] for f in self.faces])
        self.fs_nodeC = np.array([node_index_map[f.get_nodes()[2]] for f in self.faces])
        print("Defined faces.")

        # this is what things should look like if we used face springs
        # self.fs_nodeA_inc = make_incidence_matrix(self.nodes, self.face_springs, self.fs_nodeA)
        # self.fs_nodeB_inc = make_incidence_matrix(self.nodes, self.face_springs, self.fs_nodeB)
        # self.fs_nodeC_inc = make_incidence_matrix(self.nodes, self.face_springs, self.fs_nodeC)

        # apos = self.node_lastpos[self.fs_nodeA]
        # bpos = self.node_lastpos[self.fs_nodeB]
        # cpos = self.node_lastpos[self.fs_nodeC]
        # self.fs_nomAngles = np_get_angles(apos, bpos, cpos)
        # self.face_axStiff = [f.axialStiffness for f in self.face_springs]
        # print(len(self.face_springs))
        # self.face_tol = self.face_springs[0].tol #[f.tol for f in self.face_springs]
        # self.face_collapsed = np.full(len(self.face_springs), False)

        self.update_face_normals()
        print("Updated face normals.")

        self.update_crease_springs()
        print("Updated crease springs.")

        # Define the step duration
        self.timedelta = timedelta
        
        if constraint_f:
            self.constrained = True
            self.constraint_f = constraint_f
        else:
            self.constrained = False

    # update geometry for the various pieces
    def update_edge_springs(self):
        diff = self.node_lastpos[self.es_nodeB] - self.node_lastpos[self.es_nodeA]
        self.es_curL = np.linalg.norm(diff, axis = 1)
        f_dirAB = ((1/self.es_curL) * diff.T).T
        stretch = self.es_curL - self.es_nomL
        v_diff = self.node_lastvel[self.es_nodeB] - self.node_lastvel[self.es_nodeA]

        self.es_F = ((self.es_k * stretch) * f_dirAB.T).T + (self.es_d * v_diff.T).T
        

    def update_crease_springs(self):
        posA = self.node_lastpos[self.cs_nodeA]
        posB = self.node_lastpos[self.cs_nodeB]
        
        v0 = np_unit(posB - posA)

        v1 = self.node_lastpos[self.cs_nodeC] - posA
        v2 = self.node_lastpos[self.cs_nodeD] - posA

        proj1L = np.sum(v0 * v1, axis=1)
        proj2L = np.sum(v0 * v2, axis=1)

        self.crease_h1 = np.sqrt(np.sum(v1 * v1, axis=1) - proj1L * proj1L)
        self.crease_h2 = np.sqrt(np.sum(v2 * v2, axis=1) - proj2L * proj2L)

        cL = np_len(posB - posA)
        
        self.crease_coef1 = proj1L / cL
        self.crease_coef2 = proj2L / cL
        
        self.crease_n1 = self.face_normal[self.crease_face1]
        self.crease_n2 = self.face_normal[self.crease_face2]
        
        theta = np.arctan2(np.sum(np.cross(v0, self.crease_n1) * self.crease_n2, axis=1),
                           np.sum(self.crease_n1 * self.crease_n2, axis=1))
        diff = theta - self.crease_curTheta
        diff[diff < -5] += TWO_PI
        diff[diff > 5] -= TWO_PI
        self.crease_curTheta += diff

        self.crease_angF = self.crease_k * (self.crease_tgtTheta - self.crease_curTheta)

        c1h1 = ((1 - self.crease_coef1)/self.crease_h1)[:,np.newaxis]
        c2h2 = ((1 - self.crease_coef2)/self.crease_h2)[:,np.newaxis]

        term1 = self.crease_n1 * c1h1
        term2 = self.crease_n2 * c2h2
        self.cs_FA = - ((term2 + term1) * self.crease_angF[:,np.newaxis])

        c1h1 = (self.crease_coef1/self.crease_h1)[:,np.newaxis]
        c2h2 = (self.crease_coef2/self.crease_h2)[:,np.newaxis]

        term1 = self.crease_n1 * c1h1
        term2 = self.crease_n2 * c2h2
        self.cs_FB = - ((term2 + term1) * self.crease_angF[:,np.newaxis])

        self.cs_FC = (self.crease_angF / self.crease_h1)[:,np.newaxis] * self.crease_n1        
        self.cs_FD = (self.crease_angF / self.crease_h2)[:,np.newaxis] * self.crease_n2

    def update_face_normals(self):
        self.ab = self.node_lastpos[self.fs_nodeB] - self.node_lastpos[self.fs_nodeA]
        self.bc = self.node_lastpos[self.fs_nodeC] - self.node_lastpos[self.fs_nodeB]
        self.face_normal = np_unit_normal(self.ab, self.bc)
        
    # def update_face_springs(self):
    #     #ab = self.node_lastpos[self.fs_nodeB] - self.node_lastpos[self.fs_nodeA]
    #     #bc = self.node_lastpos[self.fs_nodeC] - self.node_lastpos[self.fs_nodeB]
    #     ab = self.ab
    #     bc = self.bc
    #     ca = self.node_lastpos[self.fs_nodeA] - self.node_lastpos[self.fs_nodeC]
    #     len_ab = np_len(ab)
    #     len_bc = np_len(bc)
    #     len_ca = np_len(ca)

    #     # check that lengths > tol
    #     # self.face_collapsed = [(len_ab < self.face_tol) |
    #     #                        (len_bc < self.face_tol) |
    #     #                        (len_ca < self.face_tol)]
    #     # no real point in doing this unless we later act on it...
                
    #     u_ab = np_unit(ab)
    #     u_bc = np_unit(bc)
    #     u_ca = np_unit(ca)
                
    #     # self.face_normal = np_unit_normal(ab, bc)
        
    #     angA = np_uangle(u_ab, -u_ca)
    #     angB = np_uangle(u_bc, -u_ab)
    #     angC = np_uangle(u_ca, -u_bc)
    #     self.fs_angles = np.column_stack((angA, angB, angC))

    #     fs_nxAB = (1/len_ab)[:,np.newaxis] * np.cross(self.face_normal, ab)
    #     fs_nxBC = (1/len_bc)[:,np.newaxis] * np.cross(self.face_normal, bc)
    #     fs_nxCA = (1/len_ca)[:,np.newaxis] * np.cross(self.face_normal, ca)

    #     fs_angleDiff = self.fs_angles - self.fs_nomAngles

    #     self.fs_FA = fs_angleDiff[:,0][:,np.newaxis] * (fs_nxCA + fs_nxAB)
    #     self.fs_FA += -fs_angleDiff[:,1][:,np.newaxis] * fs_nxAB
    #     self.fs_FA -= fs_angleDiff[:,2][:,np.newaxis] * fs_nxCA
    
    #     self.fs_FB = fs_angleDiff[:,1][:,np.newaxis] * (fs_nxAB + fs_nxBC)
    #     self.fs_FB += -fs_angleDiff[:,2][:,np.newaxis] * fs_nxBC
    #     self.fs_FB -= fs_angleDiff[:,0][:,np.newaxis] * fs_nxAB
        
    #     self.fs_FC = fs_angleDiff[:,2][:,np.newaxis] * (fs_nxCA + fs_nxBC)
    #     self.fs_FC += -fs_angleDiff[:,0][:,np.newaxis] * fs_nxCA
    #     self.fs_FC -= fs_angleDiff[:,1][:,np.newaxis] * fs_nxBC

    def step(self):
        self.update_geometry()
        self.update_external_forces()
        self.move_points()
        
    def update_geometry(self):
        self.update_edge_springs() 
        self.update_crease_springs()
        self.update_face_normals()
        # self.update_face_springs()

    def move_points(self):
        # reset force
        #node_F = self.node_nullF.copy()

        # add up all internal forces on the nodes
        # (this could be as bad as N^2 time because we instantiate
        # and use whole node x spring matrices instead of 
        # just the nonzero parts, but even so, the simulation
        # is faster than the non-numpy version
        # Current implementation uses scipy sparse arrays
        # to avoid running out of memory
        node_F = self.es_node_inc @ self.es_F + \
                 self.cs_nodeA_inc @ self.cs_FA + \
                 self.cs_nodeB_inc @ self.cs_FB + \
                 self.cs_nodeC_inc @ self.cs_FC + \
                 self.cs_nodeD_inc @ self.cs_FD

        # add forces over the nodes' equivalence classes
        #node_F = self.node_to_node_eq @ node_F

        # keep pinned nodes pinned
        node_F[self.node_pinned] = 0
         
        # compute the new positions
        # node_nextpos = (self.timedelta **2) * node_F / self.node_mass + \
        #                 2.0 * self.node_lastpos - self.node_lastlastpos
        node_nextpos = (self.timedelta **2) * node_F + \
                        2.0 * self.node_lastpos - self.node_lastlastpos
                          
        
        if self.constrained:
            self.constraint_f(node_nextpos)           
        
        # update previous position and velocity values
        np.copyto(self.node_lastlastpos, self.node_lastpos)
        self.node_lastpos = node_nextpos
        self.node_lastvel = (node_nextpos - self.node_lastlastpos) / self.timedelta
        
    def update_external_forces(self):
        pass
    
    def print_springs(self, edge = True, crease = True, face = True):
        if edge:
            for e in self.edge_springs:
                print(e)
        if crease:
            for c in self.creases:
                print(c)
        if face:
            for f in self.face_springs:
                print(f)
                
    def write_stl(self, model, i):
        filename = model + '_' + str(i).zfill(4) + '.stl'
        outf = open(filename, 'w')
        outf.write('solid ' + model + '\n')
        for i in range(len(self.faces)):
            normal = self.face_normal[i,:]
            pos = [self.node_lastpos[self.face_nodeA[i]],
                    self.node_lastpos[self.face_nodeB[i]],
                    self.node_lastpos[self.face_nodeC[i]]]
            outf.write('  facet normal %8.6f %8.6f %8.6f\n' % tuple(normal))
            outf.write('    outer loop\n')
            outf.write('      vertex %8.6f %8.6f %8.6f\n' % tuple(pos[0]))
            outf.write('      vertex %8.6f %8.6f %8.6f\n' % tuple(pos[1]))
            outf.write('      vertex %8.6f %8.6f %8.6f\n' % tuple(pos[2]))
            outf.write('    endloop\n')
            outf.write('  endfacet\n')
        outf.write('endsolid\n')
        outf.close()

    def write_ply(self, model, i):
        filename = model + '_' + str(i).zfill(4) + '.ply'
        outf = open(filename, 'w')
        outf.write('ply \n')
        outf.write('format ascii 1.0\n')
        outf.write('comment ' + model + '\n')
        outf.write('element vertex ' + str(len(self.nodes)) + '\n')
        outf.write('property float x\n')
        outf.write('property float y\n')
        outf.write('property float z\n')
        outf.write('element face ' + str(len(self.faces)) + '\n')
        outf.write('property list uchar int vertex_index\n')
        outf.write('end_header\n')
        for n in range(len(self.nodes)):
            pos = self.node_lastpos[n]
            outf.write('%8.6f %8.6f %8.6f\n' % tuple(pos))
        for f in self.faces: 
            nodes = [self.node_index_map[i] for i in f.get_nodes()]
            outf.write('3 %d %d %d\n' % tuple(nodes))
        outf.close()

    def write_node_coords(self, filename):
        outf = open(filename, 'w')
        outf.write('nodex nodey x y z x0 y0 z0\n')
        for n in self.nodes:
            grid_x, grid_y = n.id
            x, y, z = self.node_lastpos[self.node_index_map[n],:]
            x0, y0, z0 = n.get_origPos()
            outf.write('%d %d %8.6f %8.6f %8.6f %8.6f %8.6f %8.6f\n' %
                       (grid_x, grid_y, x, y, z, x0, y0, z0))
        outf.close()


        

        

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
    #     self.compute_edge_assignment()
    #     F["edges_assignment"] = [self.edge_assignment[e] for e in self.grid.edges]
    #     # faceOrders do this later        
    #     # now the F object should be ready to dump
    #     return json.dumps(F, indent=4)

    # def save_folded_state(self, fold_name = 'some_fold', 
    #                       outfilename = 'out.fold'):
    #     fold_string = self.make_FOLD(frame_title = fold_name)
    #     outf = open(outfilename, 'w')
    #     outf.write(fold_string)
    #     outf.close()
            
