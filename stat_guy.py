import numpy as np
import matplotlib.pyplot as plt
import argparse
from subdivider_hardly import SimplexSubdivision
import time

PARSER = argparse.ArgumentParser()
PARSER.add_argument("--dim", type=int, required=False, default=2, help="dimension of simplex")
PARSER.add_argument(
    "--subdivisions",
    type=int,
    required=False,
    default=2,
    help="barycentric subdivisions of simplex",
)
PARSER.add_argument(
    "--dont_show",
    action="store_true",
    required=False,
    help="dont show histogram",
)

PARSER.add_argument(
    "--save",
    type=str,
    required=False,
    default=None,
    help="file to save histogram",
)
PARSER.add_argument("--n", type=int, required=False, default=30, help="times to rechoose colors")
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

for d, faces in enumerate(quite_simplex.faces_by_dim):
    print("# of", str(d) + "-dimensional faces:", len(faces))
expected_mean = quite_simplex.expected_num_rainbow_faces
print(f"expected # of rainbow faces: {expected_mean}")
rainbow_face_count = []
for i in range(args.n):
    rainbow_face_count.append(quite_simplex.num_rainbow_faces)
    quite_simplex.sample_coloring()
    print(f"sampled {i + 1}/{args.n}", end="\r")
rainbow_face_count = np.array(rainbow_face_count)
rainbow_face_proportions = rainbow_face_count / n_faces

print("mean number of rainbow faces:", np.mean(rainbow_face_count))
print("stdev number of rainbow faces:", np.std(rainbow_face_count))
sem = (np.std(rainbow_face_count) * np.sqrt(args.n / (args.n - 1))) / np.sqrt(args.n)
print("std error of mean:", sem)
print("mean prop of rainbow faces:", np.mean(rainbow_face_proportions))
print("stdev prop of rainbow faces:", np.std(rainbow_face_proportions))

# expected proportion if every face had free choice of colors
limiting_p = np.prod(np.arange(dim) + 1) / ((dim + 1) ** dim)
limiting_mean = limiting_p * n_faces
limiting_variance = (limiting_p * (1 - limiting_p)) / n_faces
limiting_stdev = np.sqrt(limiting_variance)
print("limiting mean", limiting_p)
print("limiting stdev", limiting_stdev)

plt.hist(rainbow_face_count)
plt.title(f"distribution of # of rainbow faces in $\\Delta_{args.dim}$ with {sub_divisions} barycentric subdivisions")
ylim = plt.ylim()
(pop_mean,) = plt.plot([np.mean(rainbow_face_count)] * 2, ylim, linestyle="--", label="population mean $\\bar{\\mu}$")
plt.fill_betweenx(
    ylim,
    [np.mean(rainbow_face_count) - sem] * 2,
    [np.mean(rainbow_face_count) + sem] * 2,
    color=pop_mean._color,
    alpha=0.5,
    label="$\\bar{\\mu}\\pm$ Standard Error",
)
plt.plot([expected_mean] * 2, ylim, linestyle="--", label="expected mean")

plt.ylim(ylim)
plt.legend()
# max_cnt = np.max(rainbow_face_count)
# plt.bar(np.arange(max_cnt) + 1, [np.sum(rainbow_face_count == i) for i in range(1, 1 + max_cnt)])
if args.save:
    plt.savefig(args.save)
if not args.dont_show:
    plt.show()
