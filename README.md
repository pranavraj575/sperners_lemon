# [Sperner's Lemma](https://en.wikipedia.org/wiki/Sperner%27s_lemma)
**Statement:** In a Sperner coloring of a subdivided simplex ($\Delta_n$), there is an odd number of $n$-dimensional faces that are 'rainbow' (every vertex has a unique color).

### Proof Sketch
One way to prove this involves building a graph whose vertices consist of a vertex for each maximal dimensional face of the subdivided $\Delta_n$, and one specical vertex at 'infinity'.
Edges connect two vertices if their faces share a $(n-1)$-dimensional subface that is rainbow with colors $\\{1,2,...,n\\}$,
  and connect vertices to infinity if their faces have a $(n-1)$-dimensional rainbow subface with colors $\\{1,2,...,n\\}$ on the $\\{1,2,...,n\\}$ face of $\Delta_n$.
By derfintion of a Sperner coloring, the $\\{1,2,...,n\\}$ face of $\Delta_n$ is a valid Sperner coloring of $\Delta_{n-1}$, and by induction has an odd number of rainbow subfaces.
This implies that the vertex at infinity has odd degree, and thus that there an odd number of other vertices with odd degree.
It is simple to show that an $n$-dimensional face is rainbow if and only if its vertex has odd degree, and thus there are an odd number of rainbow $n$-dimensional faces.

### Visualizations
This package creates a subdivided simplex, chooses a Sperner coloring uniformly at random, then builds the associated graph.
Since this is a very visual proof for $\Delta_1$ and $\Delta_2$, in these cases the coloring, graph, and rainbow faces are plotted in matplot.
It is also relatively efficient to choose a Sperner coloring ($O(v)$ for $v$ the number of vertices), determine if each face is rainbow ($O(m)$ for $m$ the number of faces), and other similar statistics.
This lets us sample a large number of Sperner colorings and analyze these values statistically
