from constants import LANGUAGES
from constants import ConlluColumns as Cols
import json
import os
import re

os.makedirs("data/dependency-graphs", exist_ok=True)
for lang in LANGUAGES:
    print(f'Processing {lang.name}...')
    src_file = f"data/conllu/{lang.name}.conllu"
    dst_file = f"data/dependency-graphs/{lang.name}.json"
    with open(src_file, "r", encoding="utf-8") as f:
        conllu_data = f.readlines()
    # Process the CoNLL-U data and generate dependency graphs
    graphs = []
    for line in conllu_data:
        if not line.strip():
            continue
        if line.startswith("#"):
            if line.startswith("# text ="):
                graphs.append({})
            continue
        columns = line.split("\t")
        if len(columns) != Cols.TOTAL:
            print(f"WARNING: Skipping malformed line in {lang.name}: {line.strip()}")
            continue
        assert graphs, f"Data error in {lang.name}: token line before sentence text"
        id = columns[Cols.ID]
        head = columns[Cols.HEAD]
        if id.isnumeric() and head.isnumeric():
            graphs[-1][id] = head
    with open(dst_file, "w", encoding="utf-8") as f:
        json.dump(graphs, f, ensure_ascii=False, indent=4)