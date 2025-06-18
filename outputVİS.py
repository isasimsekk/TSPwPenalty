import matplotlib.pyplot as plt

# === File paths ===
input_file = "C:\\Users\\OMEN\\CLionProjects\\AlgorithmProject\\example-input-2.txt"
mst_file = "C:\\Users\\OMEN\\CLionProjects\\AlgorithmProject\\cmake-build-debug\\outputPATH.txt"
multi_file = "C:\\Users\\OMEN\\CLionProjects\\AlgorithmProject\\cmake-build-debug\\outputMULTIGRAPH.txt"
euler_file = "C:\\Users\\OMEN\\CLionProjects\\AlgorithmProject\\cmake-build-debug\\outputEULER.txt"
tsp_file = "C:\\Users\\OMEN\\CLionProjects\\AlgorithmProject\\cmake-build-debug\\outputTSP.txt"

# --- Load node coordinates ---
nodes = {}
with open(input_file, 'r') as f:
    lines = f.readlines()
    for line in lines[1:]:
        parts = line.strip().split()
        if len(parts) == 3:
            node_id, x, y = map(int, parts)
            nodes[node_id] = (x, y)

# --- Load edges ---
def read_edges(filename):
    edges = []
    with open(filename, 'r') as f:
        for line in f.readlines()[1:]:
            if '-' in line:
                u, rest = line.strip().split(' - ')
                v = rest.split()[0]
                edges.append((int(u), int(v)))
    return edges

# --- Load path ---
def read_path(filename):
    path = []
    with open(filename, 'r') as f:
        for line in f.readlines()[1:]:
            line = line.strip()
            if line.isdigit():
                path.append(int(line))
    return path

edges_mst = read_edges(mst_file)
edges_multi = read_edges(multi_file)
euler_path = read_path(euler_file)
tsp_path = read_path(tsp_file)

# === Plot 1: MST + Multigraph ===
plt.figure(figsize=(14, 10))
plt.title("MST and Multigraph with Arrows")

# Multigraph (faint arrows)
for u, v in edges_multi:
    x1, y1 = nodes[u]
    x2, y2 = nodes[v]
    dx, dy = x2 - x1, y2 - y1
    plt.arrow(x1, y1, dx, dy, head_width=2, head_length=3, fc='red', ec='red', alpha=0.2, length_includes_head=True)

# MST (solid arrows)
for u, v in edges_mst:
    x1, y1 = nodes[u]
    x2, y2 = nodes[v]
    dx, dy = x2 - x1, y2 - y1
    plt.arrow(x1, y1, dx, dy, head_width=2.5, head_length=4, fc='red', ec='red', linewidth=1.5, length_includes_head=True)

for node_id, (x, y) in nodes.items():
    plt.plot(x, y, 'ko', markersize=3)
    plt.text(x + 2, y + 2, str(node_id), fontsize=6)

plt.xlabel("X")
plt.ylabel("Y")
plt.grid(True)
plt.gca().set_aspect('equal')
plt.tight_layout()
plt.show()

# === Plot 2: Eulerian Tour ===
plt.figure(figsize=(14, 10))
plt.title("Eulerian Tour with Arrows")

for i in range(1, len(euler_path)):
    u = euler_path[i - 1]
    v = euler_path[i]
    x1, y1 = nodes[u]
    x2, y2 = nodes[v]
    dx, dy = x2 - x1, y2 - y1
    plt.arrow(x1, y1, dx, dy, head_width=2.5, head_length=4, fc='blue', ec='blue', linewidth=1, length_includes_head=True)

for node_id, (x, y) in nodes.items():
    plt.plot(x, y, 'ko', markersize=3)
    plt.text(x + 2, y + 2, str(node_id), fontsize=6)

plt.xlabel("X")
plt.ylabel("Y")
plt.grid(True)
plt.gca().set_aspect('equal')
plt.tight_layout()
plt.show()

# === Plot 3: TSP Approximation ===
plt.figure(figsize=(14, 10))
plt.title("TSP Approximate Tour with Arrows")

for i in range(1, len(tsp_path)):
    u = tsp_path[i - 1]
    v = tsp_path[i]
    x1, y1 = nodes[u]
    x2, y2 = nodes[v]
    dx, dy = x2 - x1, y2 - y1
    plt.arrow(x1, y1, dx, dy, head_width=2.5, head_length=4, fc='green', ec='green', linewidth=2, length_includes_head=True)

for node_id, (x, y) in nodes.items():
    plt.plot(x, y, 'ko', markersize=3)
    plt.text(x + 2, y + 2, str(node_id), fontsize=6)

plt.xlabel("X")
plt.ylabel("Y")
plt.grid(True)
plt.gca().set_aspect('equal')
plt.tight_layout()
plt.show()