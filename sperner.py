import itertools

import numpy as np
import matplotlib.pyplot as plt
import argparse
from subdivider_hardly import SimplexSubdivision
import time

PARSER = argparse.ArgumentParser()
PARSER.add_argument('--dim', type=int, required=False, default=2,
                    help='dimension of simplex')
PARSER.add_argument('--subdivisions', type=int, required=False, default=2,
                    help='barycentric subdivisions of simplex')
PARSER.add_argument('--invis-pts', action='store_true', required=False,
                    help='dont show points of graph')
args = PARSER.parse_args()

dim = args.dim
sub_divisions = args.subdivisions
timothy = time.time()
quite_simplex = SimplexSubdivision(
    dim=dim,
    subdivisions=sub_divisions,
    track_colors=True,
    track_adjacency=True,
)
gen_time = time.time()

for d, faces in enumerate(quite_simplex.faces_by_dim):
    print('# of', str(d) + '-dimensional faces:', len(faces))

rainbow_faces_straight_count = quite_simplex.num_rainbow_faces
print("RAINBOW FACES COUNT:", rainbow_faces_straight_count, )
print("RAINBOW FACES PROPORTION:", rainbow_faces_straight_count/len(quite_simplex.faces_by_dim[-1]))
print("RAINBOW FACES EXPECTED PROPORTION:", np.prod(np.arange(dim) + 1)/(dim + 1)**dim)

# make graph
crossing_colors = set(quite_simplex.colors[:dim])
graph_time = time.time()
# find edges in the graph connecting max dimensional faces based on if their edges are correct
collor_connections = quite_simplex.make_silly_graph(crossing_colors=crossing_colors)

incidence = np.zeros((len(quite_simplex.faces_by_dim[-1]) + 1,
                      len(quite_simplex.faces_by_dim[-1]) + 1))
face_to_idx = {face: i for i, face in enumerate(quite_simplex.faces_by_dim[-1])}
face_to_idx[quite_simplex.inf_face_idx] = len(quite_simplex.faces_by_dim[-1])
for face, facep in collor_connections:
    incidence[face_to_idx[face], face_to_idx[facep]] = 1
    incidence[face_to_idx[facep], face_to_idx[face]] = 1
incidence_time = time.time()

parity = (np.sum(incidence, axis=0)%2)[:-1]  # ignore the infinity vertex
rainbow_faces = [quite_simplex.faces_by_dim[-1][idx] for idx in np.where(parity)[0]]
print("RAINBOW FACES:")
for face in rainbow_faces:
    print(' ', face)
print('NUMBER OF RAINBOW FACES:', len(rainbow_faces))

print('subdivisions took', gen_time - timothy, 's')
print('finding edges took', graph_time - gen_time, 's')

if dim <= 2:  # geometic realization is plottable
    print('plotting graph')
    if dim == 2:
        inf_pos = np.array([0., -1.69])
    if dim == 1:
        inf_pos = np.array([1.25, -.25])
    quite_simplex.plot_coloring(
        plotter=plt,
        lines_kwargs=dict(color='black', linewidth=1, alpha=.5),
        pts_kwargs=dict(zorder=420, )
    )
    quite_simplex.plot_silly_graph(quite_simplex.make_silly_graph(crossing_colors=crossing_colors),
                                   inf_vertex_pos=inf_pos,
                                   show_pts=not args.invis_pts,
                                   pts_kwargs=dict(c='black',
                                                   zorder=420,
                                                   marker='x',
                                                   ),
                                   edge_kwargs=dict(color='orange',
                                                    linewidth=2,
                                                    alpha=1,
                                                    zorder=419,
                                                    ),
                                   )
    quite_simplex.plot_rainbow_faces(dim2kwargs=dict(color='purple',
                                                     alpha=.5,
                                                     zorder=-1, ),
                                     dim1kwargs=dict(color='purple',
                                                     linewidth=2,
                                                     alpha=1,
                                                     zorder=419,
                                                     linestyle='dotted'))
    plt.show()
