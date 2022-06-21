from math import pi

# geometry

# pleat fold offset (standard pleat width = 1)
EPS = 1e-2
# this is on the low side--some of my small EH pieces are
# closer to the ratio 1/30 (~.1  mm paper thickness, ~3mm pleat widths)


pihalf = pi/2

# fold angle
targetOff = 1.0 / 4
targetThetaV = pi * (1 - targetOff)
targetThetaM = -pi * (1 - targetOff)

constants = { 'timedelta': 0.008,
             # spring parameters
             'k_edge' : 0.7,
             'd_edge' : 0.85,
             'k_flat' : 0.7,
             'd_flat' : 0.85,
             'k_crease' : 0.7,
             'd_crease' : 0.85,
             'k_axial' : 20,
             # fold resistance
             'map_crease_targetTheta' : { 'V': targetThetaV,
                                         'M': targetThetaM,
                                         'F': 0.0 },
             'snapshot_freq': 100,
             'simulate_steps': 30000}

#solver_type = 'vectorized'


