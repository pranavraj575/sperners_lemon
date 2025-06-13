import itertools

import numpy as np
import matplotlib.pyplot as plt
import argparse

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

# first make the simplex
if dim == 2:
    V = np.array([
        [-np.sqrt(3), -1],
        [np.sqrt(3), -1],
        [0, 2],
    ])
    colors = 'red', 'blue', 'purple'
elif dim==1:
    colors = 'red', 'blue'
else:
    V = np.identity(dim + 1)
    colors = list(range(dim + 1))
V = list(V)
faces_by_dim = []
for d in range(dim + 1):
    d_dim_faces = list(itertools.combinations(range(len(V)), d + 1))
    faces_by_dim.append(d_dim_faces)

# barycentric subdivision
color_choices = [{c} for c in colors]
levels = [0 for _ in V]
for lev in range(sub_divisions):
    Vp = []
    color_choicesP = []
    faces_by_dimP = [list() for _ in faces_by_dim]
    face_to_generated_faces = {(): [()]}
    idx = 0
    levelsP = []
    for faces in faces_by_dim:
        for face in faces:
            face_to_generated_faces[face] = []
            Vp.append(sum([V[i] for i in face])/len(face))
            if len(face) == 1:
                levelsP.append(levels[face[0]])
            else:
                levelsP.append(lev + 1)
            choices = set.union(*[color_choices[i] for i in face])
            color_choicesP.append(choices)
            for faceface in itertools.combinations(face, len(face) - 1):
                for faceface_generated in face_to_generated_faces[faceface]:
                    generated = tuple(list(faceface_generated) + [idx])

                    face_to_generated_faces[face].append(generated)
                    for ss in range(1, 1 + len(generated)):
                        for thingy in itertools.combinations(generated, ss):
                            if thingy not in faces_by_dimP[len(thingy) - 1]:
                                faces_by_dimP[len(thingy) - 1].append(thingy)
            idx += 1
    V = Vp
    faces_by_dim = faces_by_dimP
    color_choices = color_choicesP
    levels = levelsP

for d, faces in enumerate(faces_by_dim):
    print('# of', str(d) + '-dimensional faces:', len(faces))

colors = [
    np.random.choice(list(ch))
    for ch in color_choices
]
crossing_colors = set(colors[:dim])

INF_FACE_idx = 'inf'
# find edges in the graph connecting max dimensional faces based on if their edges are correct
collor_connections = []
incidence = np.zeros((len(faces_by_dim[-1]) + 1, len(faces_by_dim[-1]) + 1))
face_to_idx = {face: i for i, face in enumerate(faces_by_dim[-1])}
face_to_idx[INF_FACE_idx] = len(faces_by_dim[-1])

for face, facep in itertools.combinations(faces_by_dim[-1], 2):
    shared_vertices = set.intersection(set(face), set(facep))
    if len(shared_vertices) == len(face) - 1:
        # then they are friends :)
        if crossing_colors == {colors[i] for i in shared_vertices}:
            collor_connections.append((face, facep))
            incidence[face_to_idx[face], face_to_idx[facep]] = 1
            incidence[face_to_idx[facep], face_to_idx[face]] = 1

# special faces that border the special face of triangle
for face in faces_by_dim[-1]:
    for subface in itertools.combinations(face, len(face) - 1):
        if all(color_choices[i].issubset(crossing_colors) for i in subface):
            if crossing_colors == {colors[i] for i in subface}:
                collor_connections.append((face, INF_FACE_idx))
                incidence[face_to_idx[face], face_to_idx[INF_FACE_idx]] = 1
                incidence[face_to_idx[INF_FACE_idx], face_to_idx[face]] = 1

parity = (np.sum(incidence, axis=0)%2)[:-1]  # ignore the infinity vertex
rainbow_faces = [faces_by_dim[-1][idx] for idx in np.where(parity)[0]]
print("RAINBOW FACES:")
for face in rainbow_faces:
    print(' ', face)
print('NUMBER OF RAINBOW FACES FOUND:', len(rainbow_faces))

if dim <= 2:  # geometic realization is plottable
    if dim == 2:
        inf_pos = np.array([0., -1.69])
    if dim == 1:
        inf_pos = np.array([1.25, -.25])
    E = [[V[i], V[j]] for i, j in faces_by_dim[1]]  # 1 dim faces
    plt.scatter([v[0] for v in V],
                [v[1] for v in V],
                c=colors,
                zorder=420,
                s=[40/(1.6**(lev)) for lev in levels],
                )
    for p, q in E:
        plt.plot([p[0], q[0]], [p[1], q[1]], color='black', linewidth=1, alpha=.5)

    hi_dim_face_pts = {face: sum([V[i] for i in face])/len(face) for face in faces_by_dim[-1]}
    hi_dim_face_pts[INF_FACE_idx] = inf_pos
    if not args.invis_pts:
        plt.scatter([hi_dim_face_pts[k][0] for k in hi_dim_face_pts],
                [hi_dim_face_pts[k][1] for k in hi_dim_face_pts],
                c='black',
                zorder=420,
                marker='x',
                # s=[40/(2**(lev)) for lev in levels],
                )

    for face, facep in collor_connections:
        p = hi_dim_face_pts[face]
        q = hi_dim_face_pts[facep]
        plt.plot([p[0], q[0]], [p[1], q[1]],
                 color='orange',
                 linewidth=2,
                 alpha=1,
                 zorder=419,
                 )
    for rface in rainbow_faces:
        if dim==2:
            plt.fill([V[i][0] for i in rface],
                     [V[i][1] for i in rface],
                     color='purple',
                     alpha=.5,
                     zorder=-1,
                     )
        if dim==1:
            p,q=[V[i] for i in rface]
            plt.plot([p[0], q[0]], [p[1], q[1]],
                     color='purple',
                     linewidth=2,
                     alpha=1,
                     zorder=419,
                    linestyle='dotted'
                     )
    plt.show()
