from constants import LANGUAGES
import os

os.makedirs("data/conllu", exist_ok=True)
for lang in LANGUAGES:
    src_file = f"raw_data/ud-treebanks-v2.16/UD_{lang.name}-PUD/{lang.id}_pud-ud-test.conllu"
    dst_file = f"data/conllu/{lang.name}.conllu"
    with open(src_file, "rb") as src, open(dst_file, "wb") as dst:
        dst.write(src.read())