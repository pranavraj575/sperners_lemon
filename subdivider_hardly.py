import numpy as np
import itertools
from collections import defaultdict
import matplotlib.pyplot as plt


class SimplexSubdivision:
    def __init__(self, dim=2, subdivisions=0, track_colors=True, track_adjacency=True):
        """
        :param dim:
        :param subdivisions:
        :param track_colors:
        :param track_adjacency: tracks (dim-1)-dimensional intersections betweeen dim-dimensional faces
            will be able to quickly get a list of all dim-dimensional faces that intersect a particular dim-dimensional face
            at a (dim-1)-dimensional face
        """
        # make the dim-simplex, subdivide it if specified
        # vertex stuff
        self.dim = dim
        if self.dim == 2:
            V = np.array([
                [-np.sqrt(3), -1],
                [np.sqrt(3), -1],
                [0, 2],
            ])
            all_colors = 'red', 'blue', 'purple'
        elif self.dim == 1:
            V = np.identity(self.dim + 1)
            all_colors = 'red', 'blue'
        else:
            V = np.identity(self.dim + 1)
            all_colors = list(range(self.dim + 1))
        self.V = list(V)
        # level of each vertex, (# of subdivisions before it appears)
        self.lev = 0
        self.levels = [self.lev for _ in self.V]
        # assign a color for each vertex
        self.track_colors = track_colors
        if self.track_colors:
            self.all_colors = all_colors
            self.color_choices = [{c} for c in all_colors]

        # simplex stuff
        self.faces_by_dim = []
        for d in range(self.dim + 1):
            d_dim_faces = list(itertools.combinations(range(len(V)), d + 1))
            self.faces_by_dim.append(d_dim_faces)
        self.track_adjacency = track_adjacency
        self.inf_face_idx = 'inf'
        if self.track_adjacency:
            self.adjacency_tracker = {self.faces_by_dim[-1][0]: {self.inf_face_idx},
                                      self.inf_face_idx: {self.faces_by_dim[-1][0]},
                                      }
        # there is one max dimensional face, connected to the outside

        for _ in range(subdivisions):
            self.barycentric_subdivision()

        if self.track_colors:
            self.sample_coloring()

    def barycentric_subdivision(self):
        self.lev += 1
        Vp = []
        color_choicesP = []
        faces_by_dimP = [list() for _ in self.faces_by_dim]
        faces_by_dim_set = [set() for _ in self.faces_by_dim]
        # all new subfaces that a specific face generates
        # includes the empty face so that we do not have to explicitly case a bunch of times
        face_to_generated_subfaces = {(): {()}}
        idx = 0
        levelsP = []
        # for adjacency, this is a dict of (non neighbored subface -> face)
        unmatched_boundaries = dict()

        adjacency_tracker = defaultdict(lambda: set())
        for face_dim, faces in enumerate(self.faces_by_dim):
            for face in faces:
                face_to_generated_subfaces[face] = set()
                # to generate this, we consider all generated subfaces of any k-1 dim face of face (dim k)
                #  add these subfaces, along with the cone on the vertex we just added for k
                Vp.append(sum([self.V[i] for i in face])/len(face))
                # ad vertex idx to the simplex, average of all vertices on a face
                if len(face) == 1:
                    # if face is a vertex, the level stays the same
                    levelsP.append(self.levels[face[0]])
                else:
                    # otherwise, this a newly generated vertex
                    levelsP.append(self.lev)
                if self.track_colors:
                    choices = set.union(*[self.color_choices[i] for i in face])
                    color_choicesP.append(choices)
                for faceface in itertools.combinations(face, len(face) - 1):
                    # faceface is a (face dim - 1)-dimensional subface of face
                    for faceface_generated in face_to_generated_subfaces[faceface]:
                        # iterate over all generated faces that faceface makes
                        # take the cone over this generated face and the new vertex at center of face
                        generated = tuple(list(faceface_generated) + [idx])
                        face_to_generated_subfaces[face].add(generated)
                        # the since faceface is in face, so is faceface generated
                        face_to_generated_subfaces[face].add(faceface_generated)
                        if generated not in faces_by_dim_set[len(generated) - 1]:
                            faces_by_dimP[len(generated) - 1].append(generated)
                            faces_by_dim_set[len(generated) - 1].add(generated)

                        if self.track_adjacency and (len(generated) - 1 == self.dim):
                            # tracking adjacency of this dimension
                            for gen_subface in itertools.combinations(generated, len(generated) - 1):
                                # since this is a subdivision, each generated maximal dimensional face
                                #   shares a specific face with at most ONE other generated face
                                if gen_subface in unmatched_boundaries:
                                    generated_p = unmatched_boundaries.pop(gen_subface)
                                    adjacency_tracker[generated].add((generated_p, gen_subface))
                                    adjacency_tracker[generated_p].add((generated, gen_subface))
                                else:
                                    unmatched_boundaries[gen_subface] = generated
                idx += 1
        self.V = Vp
        self.faces_by_dim = faces_by_dimP
        if self.track_colors:
            self.color_choices = color_choicesP
        if self.track_adjacency:
            for boundary_subface in unmatched_boundaries:
                outside_face = unmatched_boundaries[boundary_subface]
                adjacency_tracker[outside_face].add((self.inf_face_idx, boundary_subface))
                adjacency_tracker[self.inf_face_idx].add((outside_face, boundary_subface))
            self.adjacency_tracker = adjacency_tracker
        self.levels = levelsP

    def sample_coloring(self):
        self.colors = [
            np.random.choice(sorted(ch)) # sort so that random seed selection works
            for ch in self.color_choices
        ]

    @property
    def num_rainbow_faces(self):
        return sum(
            set(self.colors[i] for i in face) == set(self.all_colors)
            for face in self.faces_by_dim[-1]
        )

    @property
    def rainbow_faces(self):
        return [
            face
            for face in self.faces_by_dim[-1]
            if set(self.colors[i] for i in face) == set(self.all_colors)
        ]

    def make_silly_graph(self, crossing_colors=None):
        """
        Args:
            crossing_colors: list of self.dim colors to use in proof
        Returns: edge list of sperner lemma graph
        """
        if crossing_colors is None:
            crossing_colors = set(self.colors[:self.dim])
        else:
            crossing_colors = set(crossing_colors)
        collor_connections = []
        # incidence = np.zeros((len(self.faces_by_dim[-1]) + 1, len(self.faces_by_dim[-1]) + 1))
        face_to_idx = {face: i for i, face in enumerate(self.faces_by_dim[-1])}
        face_to_idx[self.inf_face_idx] = len(self.faces_by_dim[-1])
        for face in self.faces_by_dim[-1]:
            if self.track_adjacency:
                neigh = self.adjacency_tracker[face]
            else:
                neigh = [(facep, set.intersection(set(face), set(facep)))
                         for facep in self.faces_by_dim[-1]
                         if len(set.intersection(set(face), set(facep))) == len(face) - 1]
            for facep, shared_vertices in neigh:
                if len(shared_vertices) == len(face) - 1:
                    # then they are friends :)
                    if crossing_colors == {self.colors[i] for i in shared_vertices}:
                        collor_connections.append((face, facep))
                        # incidence[face_to_idx[face], face_to_idx[facep]] = 1
                        # incidence[face_to_idx[facep], face_to_idx[face]] = 1
        if self.track_adjacency:
            neigh = self.adjacency_tracker[self.inf_face_idx]
        else:
            # just check everything
            neigh = sum([[
                (face, subface)
                for subface in itertools.combinations(face, len(face) - 1)
            ]
                for face in self.faces_by_dim[-1]],
                []
            )
            # sum(list of lists,[]) concatenates lists
        # special faces that border the special face of triangle
        for outside_face, shared_vertices in neigh:
            # check that each vertex on the outside is on the 'crossing colors' face
            if all(self.color_choices[i].issubset(crossing_colors) for i in shared_vertices):
                if crossing_colors == {self.colors[i] for i in shared_vertices}:
                    # they are friends :)
                    collor_connections.append((outside_face, self.inf_face_idx))
                    # incidence[face_to_idx[outside_face], face_to_idx[self.inf_face_idx]] = 1
                    # incidence[face_to_idx[self.inf_face_idx], face_to_idx[outside_face]] = 1

        # parity = (np.sum(incidence, axis=0)%2)[:-1]  # ignore the infinity vertex

        return collor_connections  # , incidence

    def plot_coloring(self,
                      plotter=None,
                      lines_kwargs=dict(color='black',
                                        linewidth=1,
                                        alpha=.5,
                                        ),
                      pts_kwargs=dict(
                          zorder=420,
                      ),
                      ):
        assert self.dim <= 2
        if plotter is None:
            plotter = plt
        plotter.scatter([v[0] for v in self.V],
                        [v[1] for v in self.V],
                        c=self.colors,
                        s=[40/(1.6**(lev)) for lev in self.levels],
                        **pts_kwargs,
                        )
        E = [[self.V[i], self.V[j]]
             for i, j in self.faces_by_dim[1]]  # 1 dim faces
        for p, q in E:
            plotter.plot([p[0], q[0]],
                         [p[1], q[1]],
                         **lines_kwargs,
                         )

    def plot_silly_graph(self,
                         color_connections,
                         inf_vertex_pos=None,
                         show_pts=True,
                         pts_kwargs=dict(
                             c='black',
                             zorder=420,
                             marker='x',
                         ),
                         edge_kwargs=dict(
                             color='orange',
                             linewidth=2,
                             alpha=1,
                             zorder=419,
                         ),
                         ):
        assert self.dim <= 2
        if inf_vertex_pos is None:
            if self.dim == 1:
                inf_vertex_pos = np.array([1.25, -.25])
            if self.dim == 2:
                inf_vertex_pos = np.array([0., -1.69])

        hi_dim_face_pts = {face: sum([self.V[i] for i in face])/len(face)
                           for face in self.faces_by_dim[-1]}
        hi_dim_face_pts[self.inf_face_idx] = inf_vertex_pos
        if show_pts:
            plt.scatter([hi_dim_face_pts[k][0] for k in hi_dim_face_pts],
                        [hi_dim_face_pts[k][1] for k in hi_dim_face_pts],
                        **pts_kwargs,
                        )
        for face, facep in color_connections:
            p = hi_dim_face_pts[face]
            q = hi_dim_face_pts[facep]
            plt.plot([p[0], q[0]],
                     [p[1], q[1]],
                     **edge_kwargs,
                     )

    def plot_rainbow_faces(self,
                           dim2kwargs=dict(
                               color='purple',
                               alpha=.5,
                               zorder=-1,
                           ),
                           dim1kwargs=dict(
                               color='purple',
                               linewidth=2,
                               alpha=1,
                               zorder=419,
                               linestyle='dotted'
                           ),
                           ):
        assert self.dim <= 2
        for rface in self.rainbow_faces:

            if self.dim == 2:
                plt.fill([self.V[i][0] for i in rface],
                         [self.V[i][1] for i in rface],
                         **dim2kwargs,
                         )
            if self.dim == 1:
                p, q = [self.V[i] for i in rface]
                plt.plot([p[0], q[0]], [p[1], q[1]],
                         **dim1kwargs,
                         )


if __name__ == '__main__':
    import time

    for track_adjacency in True, False:
        timothy = time.time()
        stuff = SimplexSubdivision(dim=3, subdivisions=2, track_adjacency=track_adjacency, track_colors=True)
        mid_timothy = time.time()
        print('adj', track_adjacency)
        print('generation  time:', mid_timothy - timothy)
        stuff.make_silly_graph()
        print('graph makin time:', time.time() - mid_timothy)
        print('total       time:', time.time() - timothy)
        print()

    stuff = SimplexSubdivision(dim=2, subdivisions=3, track_adjacency=True, track_colors=True)
    for face, _ in stuff.adjacency_tracker[stuff.inf_face_idx]:
        plt.fill([stuff.V[i][0] for i in face],
                 [stuff.V[i][1] for i in face],
                 color='blue',
                 alpha=.25,
                 zorder=-1,
                 )
    stuff.plot_coloring(plotter=plt)
    stuff.plot_silly_graph(stuff.make_silly_graph())
    stuff.plot_rainbow_faces()
    plt.show()
