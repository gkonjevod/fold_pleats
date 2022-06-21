# Basic geometry functions useful for folding.
# Most work with arbitrary-dimension vectors,
# as long as iterating over a vector lists its coordinates.
# Not much checking is done.

from math import sqrt, atan2, acos, hypot
import numpy as np
linalg_norm = np.linalg.norm

from fold_constants import EPS, pihalf

# the vector pointing from p1 to p2
def vdir(p1, p2):
    return tuple(x2 - x1 for x1, x2 in zip(p1, p2))

# translation of point p by vector v
def shift(p, v):
    return tuple(x1 + x2 for x1, x2 in zip(p, v))
def neg_shift(p, v):
    return tuple(x1 - x2 for x1, x2 in zip(p, v))

# vector v scaled by ll
def scale(v, ll):
    return tuple(ll*x for x in v)

# the dot product of two vectors
def dot(v1, v2):
    return sum(x1 * x2 for x1, x2 in zip(v1, v2))

# unit vector in the direction of v
def unit(v, eps = 1e-6):
    #TODO: this is clearly not safe
    #l = dot(v, v)
    #if l < EPS:
    #    return 0
    return scale(v, 1.0/hypot(*v))

# cross product of two vectors (assumption: three-dimensional)
def cross(v1, v2):
    return np.array((v1[1]*v2[2] - v1[2]*v2[1], 
                     v1[2]*v2[0] - v1[0]*v2[2],
                     v1[0]*v2[1] - v1[1]*v2[0]))

# angle between two vectors
def angle(v1, v2):
    uv1 = unit(v1)
    uv2 = unit(v2)
    return acos(1.0*dot(uv1, uv2))

# reflection about an axis
# one-dimensional; axis-aligned reflection
def reflect(about, x):
    return 2*about - x

# normal of a face
# assumption: face is a triangle or at least flat in three-dimensional space
def normal(f):
    return unit(cross(vdir(f[0], f[1]), vdir(f[1], f[2])))

def normal_np(f0, f1, f2):
    cr = cross(f1 - f0, f2 - f1)
    return cr / hypot(*cr)

def dihedral_angle(n1, n2):
    crease = cross(n1, n2)
    return atan2(dot(cross(crease, n1), n2), dot(n1, n2))
    
# sense of the crease between two faces
# assumption (transitive from normal()) that the faces are flat
# used only as a last resort if the combinatorial edge assignment fails
def crease_sense(f1, f2):
    n1 = normal(f1)
    n2 = normal(f2)

    check_value = dihedral_angle(n1, n2)
    if abs(check_value) < 0.1:
        return 'F' # Flat
    elif check_value > 0:
        return 'V' # Valley
    else:
        return 'M' # Mountain

# the following are for the vectorized solver

def np_len(v):
    return linalg_norm(v, axis = 1)
    
def np_unit(v):
    v_len = linalg_norm(v, axis = 1)
    return (1/v_len)[:, np.newaxis] * v

def np_unit_normal(v0, v1):
    return np_unit(np.cross(v0, v1))

def np_uangle(u1, u2):
    return np.arccos(u1 * u2)

def np_get_angles(posA, posB, posC):
    ab = posB - posA
    bc = posC - posB
    ca = posA - posC

    u_ab = np_unit(ab)
    u_bc = np_unit(bc)
    u_ca = np_unit(ca)
    angA = np_uangle(u_ab, -u_ca)
    angB = np_uangle(u_bc, -u_ab)
    angC = np_uangle(u_ca, -u_bc)

    return np.column_stack((angA, angB, angC))
