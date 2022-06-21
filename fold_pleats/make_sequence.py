
# coding: utf-8

# Code and examples to help me understand and develop the sequence generator.
# ===========================================
# 
# I'd like human-writable formulas that generate basic pleat sequences used 
# in my tension-folded models.  At this point, I'm not trying to represent 
# locks, but only alternating sequences in the square grid.  There may be 
# ways to relax alternation, but I'll think about those a bit later.
# 
# Examples:
# ---------
# 
#     1. Basic fold  
#     2. Simple bowl
#     3. Simple pyramid
#     4. Wave
#     5. Double wave
#     6. Fancy bowl
#     7. Metaleaf
#     8. Metaleaf bowl
#     9. Metaleaf pyramid
#     10. Interleaved sequences
#     
#     
# The basic pleat sequence is a sequence of pleats.  
# The basic sequence is defined by parameters $x_0$, $x_1$, $d$, $s$, $inc$.  
# The parameters $x_0$ and $x_1$ define the start pleat and the end pleat of 
# the sequence.  The direction of each pleat in the sequence is defined by 
# $d\in\{$'v', 'h'$\}$, and its sense by $s\in\{-1,1\}$.  (Without loss of 
# generality---thanks to pleat combinations, as explained later---we can 
# assign the same direction and sense to all pleats in a basic sequence.)  
# The final parameter $inc$ defines the increment between the consecutive 
# pleats in the sequence.  Its default value is 1.  
# 
# Two pleat sequences can be combined either by concatenation or by 
# alternation.
# 
# Concatenation means literally sequence concatenation: list the pleats from 
# the two sequences in order, all pleats from the first sequence, then all 
# pleats from the second sequence.  Alternation means to construct the new 
# sequence by alternating between pleats from the first and the second
# sequence.
# 
# For example, the basic fold consists of an alternation between a horizontal- 
# and vertical-orientation basic sequences.  Sometimes we'll call such a
# sequence an arrow.  A simple bowl consists of an alternation of two arrows. 
# Varying the relative positions of start and end pleats and their sense 
# in an alternation of two arrows can generate the simple bowl, pyramid or 
# wave, depending on the sense of the pleats in the two arrows. 
# (You may notice there are four possible combinations of $s$ values and 
# I've mentioned only three different outcomes.  The fourth one is a version 
# of the basic fold constructed by an alternating pair of basic arrows.)
# 
# For a combined sequence to be consistent, the pleat sets used by its 
# constituent basic sequences must be disjoint.
# 
# One question I'm still not completely clear on is how much to work on 
# cleaning up incomplete or inconsistent specifications.  For example, 
# should I allow two basic sequences in a formula to use the same start or 
# end point, or to be defined so that when expanded on their own, the 
# pleat sets they use overlap?  There are at least two fairly difficult, 
# but somewhat complementary issues with trying to do this in an aesthetically 
# pleasing way: resolving conflicts between requested pleat sets, and using up 
# all the pleats in a portion of the sheet.  Perhaps the best course for now 
# is to leave these questions of fitting and packing for the next version and 
# focus only on sequences that have been fully and consistently specified.
# 


from itertools import chain, cycle, islice, zip_longest
from numpy.random import permutation as random_perm
from numpy import cumsum
from fold_pleats import Pleat, Grid, FlatFoldedGrid

def roundrobin(*iterables):
    "roundrobin('ABC', 'D', 'EF') --> A D E B F C"
    # Recipe credited to George Sakkis;
    # This code taken from official python docs
    pending = len(iterables)
    nexts = cycle(iter(it).__next__ for it in iterables)
    while pending:
        try:
            for next in nexts:
                yield next()
        except StopIteration:
            pending -= 1
            nexts = cycle(islice(nexts, pending))

def interleave(*iterables):
    return chain.from_iterable(zip(*iterables))

def grouper(iterable, n, fillvalue = None):
    args = [iter(iterable)] * n
    return zip_longest(*args, fillvalue = fillvalue)

SignValueSet = set([-1, 1])
DirValueSet = set(['h', 'v'])

class Abs_Pleat(object):
    ''' Abstract pleat class used to generate pleat sequences.
    A pleat is defined by its coordinate in a grid, its direction (horizontal 
    or vertical) and its sign (+1 or -1).  To align with the definitions in 
    fold_pleats, a pleat is positive if its valley fold is at a higher x value 
    than its mountain fold: self.sign = np.sign(xv - xm).  The main difference 
    between an Abs_Pleat and a Pleat is that the x coordinate of an Abs_Pleat
    is NOT the exact coordinate of either the mountain or the valley fold in 
    the crease pattern.  Instead, it is a coordinate in a grid where each 
    square corresponds to one pleat (and not just a third of it).'''
 
    def __init__(self, x, d, s):
        self.x = x
        if s not in SignValueSet:
            raise ValueError(str(s) + ' not a valid pleat sign.')
        self.sign = s
        if d not in DirValueSet:
            raise ValueError(str(d) + ' not a valid pleat direction.')
        self.dir = d
    
    def make_pleat(self, scale = 3, offset = 2, width = 1):
        ''' To make a pleat in an actual grid, we need to translate between the
        two coordinate systems. This can be a difficult problem in general,
        for several different reasons. One, adjacent +1 and -1 pleats may 
        require more spacing than the uniform unit grid. Two, we may want
        to fit the sequence into a grid of slightly inappropriate size, in
        which case it may be necessary to trim some subsequences. Three,
        the design may specify additional spacing. All of these cases should
        probably be handled by a separate class. In the simplest case, there
        is simple scale parameter that specifies the number of pleat widths
        that fit into the unit of the abstract grid and perhaps an additional
        shift at the "beginning" of the grid (at two adjacent square boundary
        edges).
        
        Let's consider for now just this setting: the folded pleat locations
        are defined by a scaling parameter (how many units of the folded grid 
        coordinate system is each abstract pleat unit?) and a shift (apply an 
        offset at the beginning of the folded grid). The two parameters are 
        called "scale" and "offset" and default to 3 and 2, respectively.
        
        This should suffice to specify any design that uses only the creases
        in a regular grid as mountain folds of the pleats, but it may require
        us to be more explicit in defining the actual pleats. To make some
        of these a little easier, I use instead three parameters: offset,
        scale and width. One abstract grid unit is equivalent to scale times
        the folded grid units. Width specifies the pleat width (the distance
        between the mountain and valley folds of each pleat).
        '''
        print('Abstract pleat', self.__str__())
        xm = (self.x - 1) * scale + offset
        xv = xm + self.sign * width
        p = Pleat(self.dir, xm, xv)
        print('Becomes Pleat(', p.dir, ',', str(p.xm), ', ', str(p.xv), ').')
        return p
    
    def __str__(self):
        return '<Abs_Pleat x:%d d:%s s:%d>' % (self.x, self.dir, self.sign)

    def __repr__(self):
        return "Abs_Pleat(%d, '%s', %d)" % (self.x, self.dir, self.sign)

    
class Sequence(object):
    '''Generic abstract sequence of pleats.  Base class for Simple Sequence,
    as well as for Combined Sequence.'''

    def __iter__(self):
        return self.make_seq()
    def __repr__(self):
        return "\n".join((repr(a) for a in self.make_seq()))

class Simple_Sequence(Sequence):
    '''A sequence of pleats with reference to an abstract grid. A simple
    sequence is defined by its 
     1. range (span: (from, to)),
     2. direction (d: 'h' or 'v'),
     3. sense (s: 1 or -1), and
     4. step size (defaults to 1).'''

    def __init__(self, span, d, s, step = 1):
        self.x0 = span[0]
        self.x1 = span[1]
        if self.x0 <= self.x1:
            stepsign = 1
        else:
            stepsign = -1
        self.sign = s
        self.dir = d
        self.step = stepsign * step

    def make_seq(self):
        #print('make_seq:', list(range(self.x0, self.x1 + self.step, self.step)))
        return (Abs_Pleat(x, self.dir, self.sign)
                for x in range(self.x0, self.x1 + self.step, self.step))

def empty_sequence():
    return iter([])    

class Combined_Sequence(Sequence):
    '''Generic class for combining sequences.  Inherited by concatenation
    and alternation.'''

    def __init__(self, *iterables):
        self.iterables = iterables

class Concat_Sequence(Combined_Sequence):
    def make_seq(self):
        return chain(*self.iterables)
        
class Altern_Sequence(Combined_Sequence):
    def make_seq(self):
        return interleave(*self.iterables)

def sequence_from_units(it, d):
    return Concat_Sequence(*[Simple_Sequence((i, i), d, 1) for i in it])

def standard_grid_size(length):
    return 3 + 3 * length

def simple_interleaved_grid_size(length):
    return 3 + 2 * 3 * length - 1

def third_width_box_pleat_grid_size(length):
    return 3 + length * 7

def quarter_width_box_pleat_grid_size(length):
    return 4 + 2 * length * (4+4)

def fold_standard_sequence(sequence_generator, 
                           sequence_name,
                           length = 4, seq_mult = 1, 
                           compute_N = standard_grid_size,
                           pleat_params = (3, 2, 1),
                           cp_scale = 1,
                           edges_to_flip = []):
    print('Making a ' + sequence_name + ' of length ', length)
    if isinstance(length, list):
        ll = sum(length)
    else:
        ll = length
    N = compute_N(ll * seq_mult)
    g = Grid(N, N)
    print('Made a grid of size', N)
    f = FlatFoldedGrid(g)
    print('made flat-folded grid')
    abstract_pleat_sequence = sequence_generator(length)
    scale, offset, width = pleat_params
    for p in abstract_pleat_sequence:
        f.pleat(p.make_pleat(scale=scale, offset=offset, width = width), coords = 'original')
    fold_name = sequence_name + '_' + str(length).zfill(3)
    outfilename = fold_name + '.fold'
    f.save_cp(fold_name, outfilename, scale = cp_scale, edges_to_flip = edges_to_flip)
    return f


# example sequence generators
    
def simple_example():
    a = Simple_Sequence((1, 5), 'v', 1)
    b = Simple_Sequence((1, 5), 'h', 1)
    c = Concat_Sequence(a, b)
    d = Altern_Sequence(a, b)
    return a, b, c, d

def bowl_example_1():
    a = Simple_Sequence((1, 7), 'v', 1)
    b = Simple_Sequence((1, 7), 'h', 1)
    c = Simple_Sequence((15, 8), 'v', -1)
    d = Simple_Sequence((15, 8), 'h', -1)
    return Altern_Sequence(Altern_Sequence(a, b),
                           Altern_Sequence(c, d))

def bowl_example_2():
    a = Simple_Sequence((1, 7), 'v', 1)
    b = Simple_Sequence((1, 7), 'h', 1)
    c = Simple_Sequence((15, 8), 'v', -1)
    d = Simple_Sequence((15, 8), 'h', -1)
    return Altern_Sequence(a, b, c, d)

def bowl_example_3():
    a = Simple_Sequence((1, 7), 'v', 1)
    b = Simple_Sequence((1, 7), 'h', 1)
    c = Simple_Sequence((15, 8), 'v', -1)
    d = Simple_Sequence((15, 8), 'h', -1)
    return Altern_Sequence(a, c, b, d)

def simple_skip2():
    a = Simple_Sequence((1, 7), 'v', 1, step = 2)
    b = Simple_Sequence((1, 7), 'h', 1, step = 2)
    c = Simple_Sequence((2, 6), 'v', 1, step = 2)
    d = Simple_Sequence((2, 6), 'h', 1, step = 2)
    return Concat_Sequence(Altern_Sequence(a, b), Altern_Sequence(c, d))

def simple_tess_0(length = [8, 8]):
    size = length[0]
    repeats = len(length)
    units = []
    print('Pleating simple tess with', repeats, 'copies of size', size)
    for i in range(size):
        hs = Simple_Sequence((1 + i, 1 + i + size * (repeats - 1)), 'h', 1, step = size)
        vs = Simple_Sequence((1 + i, 1 + i + size * (repeats - 1)), 'v', 1, step = size)
        units.append(Concat_Sequence(hs, vs))
    return Concat_Sequence(*units)
    
def simple_tess_1(length = [8, 8]):
    size = length[0]
    repeats = len(length)
    units = []
    print('Pleating simple tess with', repeats, 'copies of size', size)
    for i in range(repeats):
        hs = Simple_Sequence((1 + i * size, 1 + (i+1) * size - 1), 'h', 1, step = 1)
        vs = Simple_Sequence((1 + i * size, 1 + (i+1) * size - 1), 'v', 1, step = 1)
        units.append(Altern_Sequence(hs, vs))
    return Altern_Sequence(*units)

def simple_tess(length = [8, 8]):
    size = length[0]
    repeats = len(length)
    units = []
    print('Pleating simple tess with', repeats, 'copies of size', size)
    for i in range(repeats):
        units.append(Simple_Sequence((1 + i * size, 1 + (i+1) * size - 1), 
                                     'h', 1, step = 1))
    for i in range(repeats):
        units.append(Simple_Sequence((1 + i * size, 1 + (i+1) * size - 1), 
                                     'v', 1, step = 1))
    return Altern_Sequence(*units)


def boxpleat_spiral_pyramid(length = 8):
    print('make_sequence: boxpleat_spiral_pyramid with length = ', length)
    center_x = (2 + 3 * length)
    center_y = (2 + 3 * length)
    print('Abstract center coordinate', center_x)
    # first pair h:
    first_h = Concat_Sequence(Simple_Sequence((center_y-1, center_y-1), 'h', 1),
                              Simple_Sequence((center_y+1, center_y+1), 'h', -1))
    first_y = Concat_Sequence(Simple_Sequence((center_x-1, center_x-1), 'v', 1),
                              Simple_Sequence((center_x+1, center_x+1), 'v', -1))
    first_square = Concat_Sequence(first_h, first_y)
    # center--> bottom edge:
    # TODO: define BoxPleat_Sequence to simplify this construction
    if length == 0:
        return first_square
    seq1 = Altern_Sequence(Simple_Sequence((center_y-2, 3), 'h', -1, step = 3),
                           Simple_Sequence((center_y-4, 1), 'h', 1, step = 3))
    boxed_seq1 = grouper(seq1, 2) # this groups the two parts of the box pleat
    seq2 = Altern_Sequence(Simple_Sequence((center_x-2, 3), 'v', -1, step = 3),
                           Simple_Sequence((center_x-4, 1), 'v', 1, step = 3))
    boxed_seq2 = grouper(seq2, 2)
    seq3 = Altern_Sequence(Simple_Sequence((center_y+2, 2*center_y - 3), 'h', 1, step = 3),
                           Simple_Sequence((center_y+4, 2*center_y - 1), 'h', -1, step = 3))
    boxed_seq3 = grouper(seq3, 2)
    seq4 = Altern_Sequence(Simple_Sequence((center_x+2, 2*center_x - 3), 'v', 1, step = 3),
                           Simple_Sequence((center_x+4, 2*center_x - 1), 'v', -1, step = 3))
    boxed_seq4 = grouper(seq4, 2)
    return Concat_Sequence(first_square,
                           chain.from_iterable(Altern_Sequence(boxed_seq1, boxed_seq2, 
                                                               boxed_seq3, boxed_seq4)))

def bowl(length = 8):
    print('make sequence: bowl, length = ', length)
    return Altern_Sequence(Simple_Sequence((1, length), 'h', -1),
                           Simple_Sequence((2 * length, length + 1), 'h', 1),
                           Simple_Sequence((1, length), 'v', -1),
                           Simple_Sequence((2 * length, length + 1), 'v', 1))
    
def pyramid(length = 8):
    print('make sequence: pyramid, length = ', length)
    return Altern_Sequence(Simple_Sequence((length, 1), 'h', -1),
                           Simple_Sequence((length + 1, 2 * length), 'h', 1),
                           Simple_Sequence((length, 1), 'v', -1),
                           Simple_Sequence((length + 1, 2 * length), 'v', 1))


def in_to_out_wave(length = 8):
    print('make sequence: in to out wave, length = ', length)
    return Altern_Sequence(Simple_Sequence((length, 1), 'h', 1),
                           Simple_Sequence((length, 1), 'v', 1),
                           Simple_Sequence((length+1, 2*length), 'h', 1),
                           Simple_Sequence((length+1, 2*length), 'v', 1))
    
def skip_one_seq_wave(length = 8):
    print('make sequence: skip-one sequence based wave, length = ', length)
    return Altern_Sequence(Simple_Sequence((1, 2*length - 1), 'h', 1, step = 2),
                           Simple_Sequence((1, 2*length - 1), 'v', 1, step = 2),
                           Simple_Sequence((2*length, 2), 'h', 1, step = 2),
                           Simple_Sequence((2*length, 2), 'v', 1, step = 2))
    
def random_seq(length = 16):
    print('make sequence: random, length', length)
    v_seq = sequence_from_units(random_perm(range(1, length+1)), 'v')
    h_seq = sequence_from_units(random_perm(range(1, length+1)), 'h')
    return Altern_Sequence(v_seq, h_seq)
       
def metaleaf(length = [3, 4, 5, 6, 7]):
    print('make sequence: metaleaf with length', length)
    s = empty_sequence()
    current_loc = 1
    for i in length:
        next_part = Altern_Sequence(Simple_Sequence((current_loc, current_loc + i - 1), 'h', 1),
                                    Simple_Sequence((current_loc + i - 1, current_loc), 'v', 1))
        s = Concat_Sequence(s, next_part)
        current_loc += i 
    return s

def alternate_wave(length = [8, 8]):
    print('make sequence: alternate wave, length', length)
    v_seq = Simple_Sequence((1, sum(length)), 'v', 1)
    starts = cumsum([1] + length[:-1])
    h_pairs = [(starts[i], starts[i] + length[i] - 1) for i in range(len(length))]
    h_ordered = []
    for (i, p) in enumerate(h_pairs):
        if i % 2:
            h_ordered.append(p)
        else:
            h_ordered.append((p[1], p[0]))
    h_seq = Concat_Sequence(*[Simple_Sequence(h_ordered[i], 'h', 1) for i in range(len(length))])
    return Altern_Sequence(v_seq, h_seq)
                            
                              
def no_twist_without_rearrangements(length = 2):
    print('make sequence: basic no twist without actual rearrangements, grid size', length)
    odds_h = Altern_Sequence(Simple_Sequence((1, 4 * length - 4), 'h', 1, step = 4),
                             Simple_Sequence((2, 4 * length - 4), 'h', -1, step = 4))
    odds_v = Altern_Sequence(Simple_Sequence((1, 4 * length - 4), 'v', 1, step = 4),
                             Simple_Sequence((2, 4 * length - 4), 'v', -1, step = 4))
    evens_h = Altern_Sequence(Simple_Sequence((3, 4 * length - 2), 'h', 1, step = 4),
                              Simple_Sequence((4, 4 * length - 2), 'h', -1, step = 4))
    evens_v = Altern_Sequence(Simple_Sequence((3, 4 * length - 2), 'v', 1, step = 4),
                              Simple_Sequence((4, 4 * length - 2), 'v', -1, step = 4))
    return Concat_Sequence(odds_h, odds_v, evens_h, evens_v)
        
def simple_weave_without_rearrangements(length = 2):
    print('make sequence: simple weave without actual rearrangements, grid size', length)
    odds_h = Simple_Sequence((1, 2 * length - 2), 'h', 1, step = 2)
    odds_v = Simple_Sequence((1, 2 * length - 2), 'v', 1, step = 2)
    evens_h = Simple_Sequence((2, 2 * length - 1), 'h', 1, step = 2)
    evens_v = Simple_Sequence((2, 2 * length - 1), 'v', 1, step = 2)
    return Concat_Sequence(odds_h, odds_v, evens_h, evens_v)
         
def confusion_without_rearrangements(length = 2):
    print('make sequence: simple weave without actual rearrangements, grid size', length)
    odds_h = Simple_Sequence((1, 7 * length - 7), 'h', 1, step = 7)
    odds_v = Simple_Sequence((1, 7 * length - 7), 'v', 1, step = 7)
    evens_h = Simple_Sequence((5, 7 * length - 3), 'h', -1, step = 7)
    evens_v = Simple_Sequence((5, 7 * length - 3), 'v', -1, step = 7)
    return Concat_Sequence(odds_h, odds_v, evens_h, evens_v)
 
