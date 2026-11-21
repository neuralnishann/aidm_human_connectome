# 🧠 Human Brain Connectome Analysis

Network-science analysis of the **human brain structural connectome** from the [Human Connectome Project (HCP)](https://www.humanconnectome.org/). This is Laboratory Work 1 of the AIDM course.

The brain is modelled as an undirected graph (nodes = brain regions / parcels, edges = structural connections, each node with a 3D position). The notebook answers four questions:

1. **Graph analysis:** is the brain a *complex network*? Clustering, path lengths, degree distribution, hubs.
2. **Labeling:** can community detection / clustering recover the anatomical labels?
3. **Assortativity:** how do hubs and regions mix, globally and inside detected communities?
4. **TMS simulation:** if one region is stimulated, how far does the signal spread and what does it reach?

![Python](https://img.shields.io/badge/python-3.9%2B-blue)
![Jupyter](https://img.shields.io/badge/notebook-Jupyter-orange)
![NetworkX](https://img.shields.io/badge/graphs-NetworkX-green)
![Plotly](https://img.shields.io/badge/3D%20viz-Plotly-purple)

## Contents

- [Dataset](#dataset)
- [Helper utilities](#helper-utilities)
- [Problem 0: 3D visualisation](#problem-0-3d-visualisation)
- [Problem 1: Graph analysis (25 pts)](#problem-1-graph-analysis-25-pts)
- [Problem 2: Connectome labeling (30 pts)](#problem-2-connectome-labeling-30-pts)
- [Problem 3: Assortativity (15 pts)](#problem-3-assortativity-15-pts)
- [Problem 4: TMS simulation (30 pts)](#problem-4-tms-simulation-30-pts)
- [Getting started](#getting-started)
- [Repository structure](#repository-structure)
- [Limitations and caveats](#limitations-and-caveats)
- [References](#references)

---

## Dataset

| Item | Details |
|---|---|
| Source | Human Connectome Project structural connectome ([paper](https://link.springer.com/article/10.1007/s11571-017-9445-1)) |
| File used | `brain_scale500.graphml` (highest resolution; lower-resolution versions are in the [course Drive folder](https://drive.google.com/drive/folders/1aLmWQfYGP5jOa0cojA5zgpPP3sMeso_Li7r9KZaJ6GVSpMupYsIjeBmKa3oZ69LXiuSaWyul?usp=sharing)) |
| Graph size | 1,015 nodes in total; the **largest connected component (LCC) has 748 nodes and 11,261 edges** |
| Node attributes | `dn_position_x/y/z`, `dn_region` (cortical / subcortical), `dn_name` (anatomical area), `dn_hemisphere` |

All analyses from Problem 1.1 onwards use the LCC (`G_lcc`).

**Anatomical ground-truth labels.** Labels are built from `dn_name` after normalisation (`canonize_key`: lower-cases, strips `lh`/`rh`, drops numeric suffixes). This gives **43 anatomical labels**:

- 29 cortical areas, merged across hemispheres (for example Superiorfrontal has 82 nodes, Precentral 67, Superiorparietal 58, Postcentral 53, Inferiorparietal 47, Precuneus 44).
- 14 subcortical structures with exactly one node each (left/right thalamus, caudate, putamen, pallidum, hippocampus, amygdala, right accumbens, brain stem).

> The data file is not in this repository. Download it from the Drive link and put it next to the notebook.

---

## Helper utilities

Defined at the top of the notebook and reused throughout.

| Function | Purpose |
|---|---|
| `canonize_key(key)` | Normalises names such as `lh_superiorfrontal_12` to `Superiorfrontal` (hemisphere tags removed, digits stripped, capitalised) |
| `brain_area_function(name, map=None)` | Maps ~49 area names to a short functional description (for example *Precentral* → motor cortex, *Insula* → body sensations / emotions / pain). Falls back to the name if unknown |
| `get_neuron_info(node_data)` | Returns `(region, name, hemisphere, function)` for a node |
| `get_node_colors(attrs)` | Assigns one RGB colour per unique attribute using Matplotlib's `turbo` colormap |
| `plot_connectome(G, node_colors, show_legend)` | Interactive Plotly 3D scatter of nodes (hover shows region, name, hemisphere, function) with semi-transparent edges and an optional legend grouped by colour |

---

## Problem 0: 3D visualisation

`brain_scale500.graphml` is loaded with `nx.read_graphml`, nodes are coloured by anatomical name, and the whole connectome is drawn in 3D with `plot_connectome`. The same function is reused later to highlight hubs (Problem 1.4), cluster assignments (Problem 2) and TMS activation (Problem 4).

---

## Problem 1: Graph analysis (25 pts)

### 1.1 Global characteristics and clustering

**Method**

- Convert to undirected, keep the **largest connected component**.
- One all-pairs BFS (`nx.all_pairs_shortest_path_length`) gives per-node eccentricity, hence **diameter**, **radius** and the **average shortest path** at once.
- Clustering: `nx.average_clustering` (mean of local coefficients) and `nx.transitivity` (global).
- **Small-world check against an Erdős–Rényi reference:** `C_rand = <k>/n`, `L_rand = ln(n)/ln(<k>)`, `σ = (C/C_rand) / (L/L_rand)`.
- Plots: histogram of local clustering coefficients (with mean and median lines) and clustering vs. degree (linear and log x-axis, with a log-binned mean trend and Pearson correlation).

**Results**

| Metric | Value |
|---|---|
| Nodes (LCC / total) | 748 / 1015 |
| Edges | 11,261 |
| Mean degree | 30.11 |
| Average clustering | 0.639 |
| Transitivity | 0.443 |
| Diameter / radius | 9 / 5 |
| Average shortest path | 3.097 |
| `C_rand` / `L_rand` | 0.0403 / 1.94 |
| Small-world σ | 9.96 |
| Degree vs. clustering (Pearson r) | -0.455 (p ≈ 1.9e-39) |

### 1.2 Path length

- Global average shortest path: `nx.average_shortest_path_length` = **3.0968** (matches the BFS-based value).
- Per-vertex average distance to all *other* vertices, plotted against degree. Degrees are grouped into **log-spaced bins** (bins with fewer than 3 vertices are skipped).
- Two views: a scatter of all vertices (binned mean + graph average) and a histogram-style bar chart of mean path length per degree range (whiskers = std).

### 1.3 Degree distribution

- Histogram of degrees (min 1, max 205, mean 30.11).
- **SciPy fits + Kolmogorov–Smirnov tests** for power-law, exponential and normal distributions:

| Distribution | KS D | p-value |
|---|---|---|
| Power-law | 0.3017 | 0.0000 |
| Exponential | 0.0594 | 0.0098 |
| Normal | 0.1388 | 0.0000 |

- **`powerlaw` package** (discrete fit with automatic `x_min`): log-log CCDF of the empirical data against power-law, exponential and lognormal fits.

| Quantity | Value |
|---|---|
| Power-law exponent α | 2.978 |
| `x_min` | 32 |
| KS D of the fit | 0.134 |
| Power-law vs. exponential | R = -20.8, p ≈ 0 |
| Power-law vs. lognormal | R = -24.3, p ≈ 0 |

> **How to read R.** In the `powerlaw` package, R > 0 means the *first* distribution (power law) is the better fit and R < 0 means the *second* one is. Here R is negative for both comparisons, so exponential and lognormal describe the tail **better** than a pure power law, even though the tail is heavy and the degree distribution is clearly far from Gaussian. The degrees are therefore best described as *heavy-tailed with hubs*, not as a strict scale-free power law.

### 1.4 Most important nodes (hubs)

Centrality measures computed on the LCC:

- degree centrality
- betweenness centrality (normalised)
- eigenvector centrality (`max_iter=500`)
- closeness centrality
- a **combined score** (sum of degree, betweenness and eigenvector) and a **harmonic mean** of betweenness, closeness and eigenvector (set to 0 if any component is 0)

For each measure the top-10 nodes are printed with region, name, hemisphere and function, then highlighted in red on the 3D connectome.

**Top 10 by combined score**

| Rank | Node | Area | Function |
|---|---|---|---|
| 1 | 1008 | Left thalamus | Relay of sensory signals to cortex |
| 2 | 502 | Right thalamus | Relay of sensory signals to cortex |
| 3 | 503 | Right caudate | Movement control, learning |
| 4 | 1009 | Left caudate | Movement control, learning |
| 5 | 504 | Right putamen | Movement, skill formation |
| 6 | 1010 | Left putamen | Movement, skill formation |
| 7 | 770 | Left superior parietal | Spatial perception, attention |
| 8 | 505 | Right pallidum | Regulation of motor functions |
| 9 | 487 | Right insula | Body sensations, emotion, pain |
| 10 | 285 | Right inferior parietal | Integration of hearing, vision, attention |

Node **1008** (left thalamus) ranks first by betweenness and by the combined / harmonic-mean scores. Node **502** (right thalamus) ranks first by eigenvector centrality. Degree centrality peaks at the right caudate (503). Subcortical structures dominate every ranking; cortical hubs include superior parietal, insula, inferior parietal, lingual, supramarginal and caudal middle frontal areas.

### 1.5 Is the brain a complex network?

Yes. The notebook argues from four properties:

1. **Small-world:** short paths (3.10) combined with high clustering (0.639).
2. **Hierarchical organisation:** negative degree-clustering correlation (-0.455), meaning hubs bridge local clusters.
3. **Heavy-tailed degrees:** a few nodes with very high degree (up to 205).
4. **Hub neurons:** thalamus, caudate and putamen carry most shortest paths and influence.

---

## Problem 2: Connectome labeling (30 pts)

**Goal:** partition the LCC into groups automatically, tune the parameters, and compare against the 43 anatomical labels and against each other.

### Algorithms and parameters

`run_clustering` runs 12 methods:

| Method | Implementation | Parameters |
|---|---|---|
| Louvain | `python-louvain` `best_partition` | resolution swept over `np.linspace(0.5, 2.0, 10)`, best chosen by ARI against the anatomical labels (selected: 1.8333) |
| Leiden | `cdlib.algorithms.leiden` | library defaults |
| DeepWalk (node2vec) + GMM | `node2vec` with `p=q=1`, then `GaussianMixture` | 128 dims, walk length 40, 20 walks per node, window 10, 43 spherical components |
| Node2Vec | `node2vec` then `SpectralClustering` | 128 dims, walk length 80, 10 walks, window 10, RBF affinity, 43 clusters |
| Hierarchical: ward / complete / average | `AgglomerativeClustering` on the dense adjacency matrix | 43 clusters each |
| Async-LPA | `networkx.asyn_lpa_communities` | defaults (stochastic, single run) |
| LPA | `networkx.label_propagation_communities` | defaults |
| MCL | `markov_clustering` | inflation 2.0 (range 1.8 to 2.2 tried: lower gave over 30 clusters, higher gave too few) |
| Edge betweenness | `networkx.girvan_newman` | stop at the first split with up to 43 communities; 60 s time limit |
| EgoSplitting | `cdlib.algorithms.ego_networks` | defaults; 30 s time limit |

Each run is wrapped in try/except, and a failed or timed-out method falls back to a single cluster. A 2D **t-SNE** projection of the DeepWalk embedding is used as the common layout for the cluster plots.

### Evaluation

`evaluate_clusters` returns:

- **ARI** and **NMI** against the anatomical labels
- **Modularity** and **coverage** (`networkx` community functions)
- **Silhouette**, computed on the node features `[degree, local clustering]` (not on the embeddings), set to -1 when only one cluster exists
- **Clusters**, and **VS Anatomy** = 43 − number of clusters

### Results (sorted by modularity)

| Algorithm | Param | ARI | NMI | Modularity | Coverage | Silhouette | Clusters |
|---|---|---|---|---|---|---|---|
| Leiden | - | 0.0770 | 0.2556 | **0.6235** | 0.8745 | -0.0467 | 4 |
| Louvain | 1.8333 | 0.1390 | 0.3909 | 0.6036 | 0.7451 | -0.3489 | 9 |
| Async-LPA | - | 0.0871 | 0.3186 | 0.6027 | 0.8194 | -0.3777 | 8 |
| LPA | - | 0.0502 | 0.2187 | 0.5850 | 0.8813 | -0.3196 | 6 |
| MCL | 2.0 | 0.0131 | 0.1602 | 0.5283 | 0.9504 | -0.3769 | 13 |
| Node2Vec | - | 0.1275 | 0.4938 | 0.4755 | 0.5809 | -0.4650 | 43 |
| DeepWalk + GMM | - | **0.1887** | **0.5701** | 0.3548 | 0.4058 | -0.5545 | 43 |
| Hierarchical-complete | - | 0.0566 | 0.4116 | 0.2113 | 0.2502 | -0.1960 | 43 |
| Hierarchical-ward | - | 0.0732 | 0.4611 | 0.1939 | 0.2207 | -0.2820 | 43 |
| Hierarchical-average | - | 0.0245 | 0.3245 | 0.1626 | 0.2747 | 0.0394 | 43 |
| EgoSplitting | - | 0.0606 | 0.6181 | 0.0795 | 0.0863 | -0.8696 | 264 |
| EdgeBetweenness | - | 0.0000 | 0.0000 | 0.0000 | 1.0000 | -1.0000 | 1 |

### Consensus clustering

A co-association matrix is built from the top 3 methods by modularity (Leiden, Louvain, Async-LPA) and clustered with `SpectralClustering(n_clusters=15, affinity='precomputed')`:

| ARI | NMI | Modularity | Coverage | Silhouette | Clusters |
|---|---|---|---|---|---|
| 0.1437 | 0.4175 | 0.5576 | 0.6809 | -0.4021 | 15 |

### Findings

- No method reproduces the anatomical labels well (ARI stays below 0.19).
- **Modularity vs. anatomy trade-off:** Leiden has the highest modularity but only 4 near-equal clusters; Louvain is the most stable structural method and the closest to anatomy; embedding methods match anatomy best but have low modularity and negative silhouette.
- Hierarchical methods give low modularity (≤ 0.22). MCL, LPA, EgoSplitting (264 micro-clusters) and edge betweenness are unreliable here.
- The edge-betweenness result (one cluster, silhouette -1) matches the notebook's fallback for a failed or timed-out run, so it most likely reflects the 60 s limit rather than a real partition.

---

## Problem 3: Assortativity (15 pts)

### Procedure

1. Attach the anatomical label as node attribute `regional_label`.
2. **Global:** `degree_assortativity_coefficient` and `attribute_assortativity_coefficient` on the whole LCC.
3. **Local:** for every community of every method, take the induced subgraph. Only communities with **more than 5 nodes and more than 5 edges** are scored for degree and regional assortativity.
4. Methods analysed: Louvain, Leiden, Hierarchical (ward / complete / average), Async-LPA, MCL, EdgeBetweenness, EgoSplitting.
5. Six plots (A to F): boxplot of degree assortativity by method, regional assortativity histogram, size vs. assortativity scatter, regional vs. degree scatter, violin plot, and global vs. local bars.

### Results

| Scope | Degree assortativity | Regional assortativity |
|---|---|---|
| **Whole connectome** | **-0.024** | **0.186** |
| Leiden (mean over communities) | -0.0380 | 0.1594 |
| Louvain | -0.1177 | 0.0969 |
| Async-LPA | -0.0874 | 0.1249 |
| Ward | -0.2155 | -0.0279 |
| Complete | -0.2077 | -0.0408 |
| Average | -0.1885 | -0.0344 |
| MCL | -0.2084 | 0.0944 |
| EdgeBetweenness (single community) | -0.0236 | 0.1856 |
| EgoSplitting | -0.4237 | -0.1861 |

### Interpretation

- Hubs tend to link to lower-degree nodes (slight disassortativity), supporting integration across modules.
- Nodes of the same anatomical area are more connected to each other than chance (positive regional assortativity).
- Large communities (Leiden, Louvain, Async-LPA) show positive regional assortativity of about 0.10 to 0.19; small hierarchical or ego clusters are noisy and often negative.
- Edge betweenness gives one 748-node community, so its local numbers just repeat the global ones.

---

## Problem 4: TMS simulation (30 pts)

[Transcranial magnetic stimulation](https://en.wikipedia.org/wiki/Transcranial_magnetic_stimulation) is simulated by injecting "activation" into a target region and spreading it over the weighted connectome.

### Edge weighting

For every edge, `weight = 1 / (Euclidean distance between node positions + 1e-5)`, then **min-max normalised** to [0, 1]. Short connections get large weights.

### Propagation models

| Model | Mechanism | Parameter grid |
|---|---|---|
| **Heat diffusion** | `exp(-L·t)` (sparse `expm`) applied to a one-hot vector on the target nodes, `L` = unnormalised weighted Laplacian, result normalised by max | `time` ∈ {0.3, 0.5, 1.0} |
| **Random walk with restart** | Walkers start at random target nodes and move proportionally to edge weights; restart at a target with probability `p`; visit counts are the activation | `steps` ∈ {500, 2000}, `walk_length` ∈ {5, 15}, `restart_prob` ∈ {0.1, 0.2} |
| **Message passing** (SIR-like) | Targets start infected; an infected node infects a susceptible neighbour with probability `infection_rate × weight` and recovers with `recovery_rate` | `infection_rate` ∈ {0.5, 0.9}, `recovery_rate` ∈ {0.1, 0.3}, `iterations` ∈ {3, 7} |
| **Betweenness propagation** | Targets set to 1; nodes whose weighted betweenness exceeds a threshold push `betweenness × weight` to their neighbours | `threshold` ∈ {0.005, 0.01, 0.02} |

All activations are normalised to [0, 1]. `plot_activation` renders them in 3D (red, with opacity proportional to activation above a threshold of 0.2) and `analyze_regional_impact` summarises the functions of the top-10 activated nodes.

### Stimulation targets and expected regions

| Condition | Target (hemisphere) | Expected regions for validation |
|---|---|---|
| Stroke (motor) | Precentral (L) | Precentral, Postcentral, Superiorparietal |
| Depression | Caudal middle frontal (L) | Caudalmiddlefrontal, Amygdala, Anteriorcingulate |
| Parkinson's | Precentral + Superior frontal (L) | Precentral, Superiorfrontal, Putamen |
| Tinnitus | Superior temporal (R) | Superiortemporal, Transversetemporal, Insula |
| Schizophrenia | Caudal middle frontal + Superior temporal (L) | Caudalmiddlefrontal, Superiortemporal, Anteriorcingulate |
| PTSD | Caudal middle frontal (L) | Caudalmiddlefrontal, Amygdala, Insula |
| Eating disorder | Lateral orbitofrontal (L) | Lateralorbitofrontal, Amygdala, Insula |
| Aphasia | Pars opercularis + Superior temporal (L) | Parsopercularis, Superiortemporal, Superiorfrontal |

Targets are motivated by the rTMS guidelines paper listed in the references.

### Validation and parameter search

For every condition x model x parameter combination, the **unique areas among the top-5 most activated nodes** are compared to the expected set:

```
score = max(0, matches / |expected| - 0.15 × (number of unexpected areas))
```

The best parameters per model are kept, then the best model per condition.

### Results

| Condition | Best model | Parameters | Score | Key regions |
|---|---|---|---|---|
| Stroke (motor) | Heat diffusion | t = 0.3 | 33.3% | Precentral |
| Depression | Heat diffusion | t = 0.3 | 33.3% | Caudalmiddlefrontal |
| Parkinson's | Heat diffusion | t = 0.3 | 66.7% | Superiorfrontal, Precentral |
| Tinnitus | Random walk | 500 steps, length 5, restart 0.1 | **100%** | Superiortemporal, Transversetemporal, Insula |
| Schizophrenia | Heat diffusion | t = 0.3 | 66.7% | Superiortemporal, Caudalmiddlefrontal |
| PTSD | Heat diffusion | t = 0.3 | 33.3% | Caudalmiddlefrontal |
| Eating disorder | Heat diffusion | t = 1.0 | 51.7% | Lateralorbitofrontal, Left-putamen, Insula |
| Aphasia | Heat diffusion | t = 0.3 | 66.7% | Superiortemporal, Parsopercularis |

Heat diffusion with a short diffusion time wins for 7 of 8 conditions and mostly stays local around the stimulated area, while message passing spreads the farthest (it often reaches precuneus, superior frontal, posterior cingulate and other distant areas).

---

## Getting started

### Install

```bash
git clone https://github.com/<your-username>/<your-repo>.git
cd <your-repo>

python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate

pip install networkx numpy pandas scipy scikit-learn matplotlib seaborn plotly \
            powerlaw cdlib leidenalg python-louvain markov-clustering node2vec jupyter
```

The notebook also contains `!pip install ...` cells for these packages (Problems 1.3 and 2).

### Data

Place `brain_scale500.graphml` in the project root. It is read from `./brain_scale500.graphml`.

### Run

```bash
jupyter notebook AIDM_lab_1_nishan.ipynb
```

Run the cells **top to bottom**; later sections reuse earlier objects:

| Object | Created in | Used in |
|---|---|---|
| `G_lcc` | Problem 1 | Problems 2, 3, 4 |
| `original_labels` | Problem 2 | Problem 2 evaluation, Problem 3 |
| `cluster_results`, `results_df` | Problem 2 | Problem 3 |
| `G_lcc` edge `weight` attribute | Problem 4 | Problem 4 models |

**Runtime and platform notes**

- Betweenness centrality, all-pairs BFS, Node2Vec and the matrix exponential are the slowest steps.
- The Girvan–Newman / EgoSplitting time limits use `signal.SIGALRM`, which works on **Linux and macOS only** (on Windows run the notebook in WSL or remove the time limits).
- Plotly 3D figures are interactive in Jupyter.

---

## Repository structure

```
.
├── AIDM_lab_1_nishan.ipynb    # code, plots and written analysis
├── brain_scale500.graphml     # data (download separately)
└── README.md
```

---

## Limitations and caveats

- **Nodes are regions, not neurons.** The notebook calls nodes "neurons"; at this resolution each is a small brain parcel.
- **Louvain tuning uses the ground truth.** The resolution was chosen by maximising ARI against the anatomical labels, so its ARI is optimistic.
- **Silhouette uses `[degree, clustering]` features**, not the embedding or connectivity, so it says little about community quality for graph methods.
- **TMS validation is lenient and partly unreachable.** The stimulated region itself counts as a match when it is in the expected set. Some expected names (`Amygdala`, `Putamen`, `Anteriorcingulate`) do not equal the label strings in the graph (`Left-amygdala`, `Left-putamen`, `Caudalanteriorcingulate`), so they can never match. Parameters are also tuned on the same score they are reported with. Treat the scores as qualitative.
- **Randomness.** Async-LPA, Node2Vec, GMM, t-SNE, random walk and message passing are stochastic and no seeds are set, so numbers can change slightly between runs.
- **Rendering.** The saved notebook is about 20 MB because of embedded Plotly figures; GitHub may not render it. Use [nbviewer](https://nbviewer.org/) or clear outputs before committing.

---

## References

- Human Connectome Project: <https://www.humanconnectome.org/>
- Connectome paper: <https://link.springer.com/article/10.1007/s11571-017-9445-1>
- Lefaucheur et al., *Evidence-based guidelines on the therapeutic use of repetitive transcranial magnetic stimulation (rTMS)*: <https://www.sciencedirect.com/science/article/pii/S1388245719312799>
- Brain region atlas: <https://atlas.brain-map.org/atlas?atlas=265297126>
- Alstott et al., `powerlaw` package, for the distribution comparison method.

## Author

Nishan, AIDM Laboratory Work 1.
