from constants import LANGUAGES
import json
import os
import networkx as nx
from collections import Counter


os.makedirs("data/degree-sequences", exist_ok=True)
for lang in LANGUAGES:
    print(f'Processing {lang.name}...')
    src_file = f"data/dependency-graphs/{lang.name}.json"
    dst_file = f"data/degree-sequences/{lang.name}_degree_sequence.txt"

    with open(src_file, "r", encoding="utf-8") as f:
        dependency_graphs = json.load(f)

    G = nx.Graph() # Create an undirected graph
    for i, graph in enumerate(dependency_graphs):
        for child, parent in graph.items():
            G.add_edge(f'{i}-{child}', f'{i}-{parent}')
    
    with open(dst_file, "w", encoding="utf-8") as f:
        for node in G.nodes():
            f.write(f"{G.degree(node)}\n")

    #degree_counter = Counter(G.degree(node) for node in G.nodes())


