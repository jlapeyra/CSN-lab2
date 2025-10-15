import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from languages import LANGUAGES, Language

ordered_language_names = [lang.name for lang in LANGUAGES]
print(*ordered_language_names, sep=", ")


# --- 1. Load table ---
# Read fixed-width formatted text, skipping the separator line
df = pd.read_fwf("plot/tables/log-likelihood.txt")

# --- 3. Set language as index ---
df = df.set_index('Language')

df = df.sort_index(key=lambda x: [ordered_language_names.index(name) for name in x])

# --- 4. Convert all numeric columns to floats (safe conversion) ---
df = df.apply(pd.to_numeric, errors='coerce')

# --- 5. Plot heatmap ---
plt.figure(figsize=(10, 8))

#df_log = np.log10(df.replace(0, np.nan))
sns.heatmap(df, annot=True, fmt=".2f", cbar_kws={'label': 'log-likelihood'})
plt.title("log-likelihood Heatmap")
plt.ylabel("Language")
plt.xlabel("Model")
plt.tight_layout()
plt.show() 
