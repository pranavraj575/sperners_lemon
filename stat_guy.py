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

PARSER.add_argument('--n', type=int, required=False, default=30,
                    help='times to rechoose colors')
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
n_faces = len(quite_simplex.faces_by_dim[-1])
gen_time = time.time()
limiting_p = np.prod(np.arange(dim) + 1)/((dim + 1)**dim)
limiting_mean = limiting_p*n_faces
limiting_variance = (limiting_p*(1 - limiting_p))/n_faces
limiting_stdev = np.sqrt(limiting_variance)

for d, faces in enumerate(quite_simplex.faces_by_dim):
    print('# of', str(d) + '-dimensional faces:', len(faces))
rainbow_face_count = []
for _ in range(args.n):
    rainbow_face_count.append(quite_simplex.num_rainbow_faces)
    quite_simplex.sample_coloring()
rainbow_face_count = np.array(rainbow_face_count)
rainbow_face_proportions = rainbow_face_count/n_faces

print('found mean', np.mean(rainbow_face_proportions))
print('limiting mean', limiting_p)
print('found stdev', np.std(rainbow_face_proportions))
print('limiting stdev', limiting_stdev)
plt.hist(rainbow_face_count)
# max_cnt = np.max(rainbow_face_count)
# plt.bar(np.arange(max_cnt) + 1, [np.sum(rainbow_face_count == i) for i in range(1, 1 + max_cnt)])
plt.show()
