import pandas as pd
import numpy as np
import networkx as nx

# Load your CSV - REPLACE 'your_file.csv' with actual filename
df = pd.read_csv(r'D:\Desktop\FYDP\FYDP 3\New folder (2)\DNE-v1.1.0\rui-yan-DNE-afa94dd\data\c_elegans\edge_list.csv')
# Check structure
print("Columns:", df.columns.tolist())
print(df.head())

# Create graph (adjust column names as needed)
# Try one of these:
try:
    graph = nx.from_pandas_edgelist(df, 'source', 'target')
except:
    try:
        graph = nx.from_pandas_edgelist(df, 'node1', 'node2')
    except:
        # If it's an adjacency matrix
        adj_matrix = df.to_numpy()
        graph = nx.from_numpy_array(adj_matrix)

# Get adjacency matrix
adj = nx.adjacency_matrix(graph).toarray()

# CN calculation
cn_counts = adj @ adj 
non_edge_mask = ~adj.astype(bool)
non_edge_cn = cn_counts[non_edge_mask]

print("\n=== CN Percentile Results ===")
print(f"95th percentile CN among non-edges: {np.percentile(non_edge_cn, 95)}")
print(f"99th percentile CN among non-edges: {np.percentile(non_edge_cn, 99)}")
print(f"Max CN among non-edges: {np.max(non_edge_cn)}")
print(f"Mean CN among non-edges: {np.mean(non_edge_cn):.2f}")
print(f"Total non-edges: {len(non_edge_cn)}")