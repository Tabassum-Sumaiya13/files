# code/diagnose_split.py
import sys
sys.path.append(".")
from collections import Counter
from sklearn import model_selection
from dataset import GraphDataset
from utils.edge_splitter import EdgeSplitter

def split_train_test_edges(graph, test_size, seed):
    edge_splitter_test = EdgeSplitter(graph)
    graph_test, X_test, Y_test = edge_splitter_test.train_test_split(
        p=0.1, keep_connected=False, seed=seed
    )
    edge_splitter_train = EdgeSplitter(graph_test)
    graph_train, X, Y = edge_splitter_train.train_test_split(
        p=0.1, keep_connected=False, seed=seed
    )
    X_train, X_valid, Y_train, Y_valid = model_selection.train_test_split(
        X, Y, test_size=test_size, random_state=seed
    )
    return graph_train, X_test, X_train, X_valid

graph_data = GraphDataset("../data")
graph_data.load_graph("c_elegans", add_feats=False)
graph = graph_data.graph

degs = [d for _, d in graph.degree()]
print(f"Original: {graph.number_of_nodes()} nodes, {graph.number_of_edges()} edges")
print(f"Degree-1 nodes: {sum(1 for d in degs if d==1)}, "
      f"degree-2 nodes: {sum(1 for d in degs if d==2)}\n")

for seed in range(5):
    graph_train, X_test, X_train, X_valid = split_train_test_edges(graph, test_size=0.6, seed=seed)
    isolated = {n for n in graph_train.nodes() if graph_train.degree(n) == 0}
    touching = sum(1 for u, v in X_test if u in isolated or v in isolated)
    print(f"seed={seed}: isolated_nodes={len(isolated)}  "
          f"test_edges_touching_isolated={touching}/{len(X_test)}  "
          f"train_examples={len(X_train)}  valid_examples={len(X_valid)}")

    # --- hub-node diagnostic ---
    node_freq_in_test = Counter()
    for u, v in X_test:
        node_freq_in_test[u] += 1
        node_freq_in_test[v] += 1
    isolated_hub_damage = sum(node_freq_in_test[n] for n in isolated if n in node_freq_in_test)
    max_isolated_freq = max((node_freq_in_test[n] for n in isolated if n in node_freq_in_test), default=0)
    print(f"         isolated nodes' total appearances in test edges = {isolated_hub_damage}  "
          f"(max single-node freq = {max_isolated_freq})")