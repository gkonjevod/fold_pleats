# fold_pleats
Code for simulating folded pleat tessellations

Originally developed in 2017 for a 7OSME paper. The physics model is based 
directly on Amanda Ghassaei's Origami Simulator. Written in pure python 
with only a few dependencies.  Uses numpy and scipy to speed up low-level 
computation, and Eppstein's implementation of union-find from PADS (this
last bit is overkill for most of my examples because the equivalence
classes are small.)

Unlike Origami Simulator, this has no GUI. The pattern, with creases
and their target angles needs to be specified ahead of time and is static.

The first component (fold_pleats) lets us model and fold axis-aligned pleats.
This is done symbolically, keeping the topology of the folded sheet and the
physical locations in 3d space in separate data structures. 

After a sequence of pleats has been folded, the resulting model
can be saved into a file using the FOLD format.  This allows visualization
and importing into tools such as Amanda Ghassaei's interactive Origami
Simulator.

The next major component is a physical point-mass-spring system based on
a python reimplementation of Amanda Ghassaei's javascript code, with a few
changes, mostly to the organization of the calculations.  This part includes
code to build a mesh from the folded sheet (currently not fully general,
it makes some assumptions that are valid as long as we're dealing with
square grids and simple axis-aligned pleats), as well as to triangulate
it.

Finally, as an intermediate step, we can calculate the free complex of
a sequence of pleats.  Intuitively, this is a 2-complex consisting of nodes,
edges and faces of the folded model, where sets of nodes held together by
folds are merged.  Each such supernode is incident on all the edges and faces
that its member nodes were incident in the original folded sheet.  The forces
on each original node are still calculated exactly as they would be if all nodes
were free to move independently.  However, before actually moving the nodes,
one of two tricks are applied to prevent independent movement of nodes that
are trapped together by the folded pleats. In the original version of the code,
the trick was to simply add up the forces acting on all members of each supernode
and then displace all the nodes that comprise a supernode by exactly the same vector.
This doesn't do quite the right thing when the layers bend and the axis containing the
parts of a super node rotates. In the current version, there is an alternative
(which seems to work better and hasn't caused any problems in the examples I've run
so far): the individual members of each supernode are joined by a complete graph of
stiff springs, so that they are hopefully kept close together in roughly the same
relative spatial arrangement as the simulation progresses. The spring constant is
set as the reciprocal of the distance between each pair of nodes. Both versions
can then reuse the simulation code without any changes, the original by redefining
the functions that access force setting and retrieval from nodes so that each node
forwards its calls to a representative of the equivalence class (supernode) it
belongs to, and the current by simply including the new springs in the list
of structural springs that roughly preserve internode distances in the mesh.

In addition to the calculation of the free complex, there is also (very simple)
code that allows us to glue together subsets of selected nodes and run the
simulation with the resulting constraint.  Specific implementations of the idea
currently included are the gluing of boundary edges of the folded sheet (similar
to what I do with real paper when I need to lock a folded edge by folding it over
or by applying actual glue).

A few years after the original code was put together, I spent some time on a more
intuitive approach to generating the crease patterns. In the case of pleat tessellations,
this means mostly the specification of sequences of pleats that the code then
simulates. It is not quite as simple as I'd like it to be, but the examples included
in make_sequence do go quite a bit further than was easy with the original
method of specifying the pleat coordinates by hand.

Overall, the code is still messy.  It appears to work, but can clearly be improved both
for readability and for efficiency. However, it seems to be relatively efficient
for pure python and even though several data structures are built on top of each
other to represent different levels of detail about the crease pattern, folded
state and the simulation mesh, it doesn't require too much memory. I've been able
to run examples based on 256x256 pleats, resulting in a mesh with almost 600k nodes, 
which takes up about 20GB of memory and runs at betwen 1 and 2 iterations a second.
