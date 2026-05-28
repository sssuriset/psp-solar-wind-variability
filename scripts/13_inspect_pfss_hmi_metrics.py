from pathlib import Path

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd

base = Path(__file__).resolve().parents[1]

metrics_path = base / "data" / "processed" / "pfss" / "hmi" / "pfss_hmi_selected_crs_rss2.5_metrics_summary.csv"
figdir = base / "outputs" / "figures" / "pfss" / "hmi"
table_dir = base / "outputs" / "tables"

figdir.mkdir(parents=True, exist_ok=True)
table_dir.mkdir(parents=True, exist_ok=True)

if not metrics_path.exists():
    raise FileNotFoundError(f"Missing metrics summary: {metrics_path}")

df = pd.read_csv(metrics_path).sort_values("cr")

print("PFSS HMI metrics summary:")
print(df[[
    "cr",
    "rss",
    "nrho",
    "input_mean_abs_br",
    "source_surface_min_br",
    "source_surface_max_br",
    "source_surface_mean_abs_br",
    "source_surface_finite_fraction",
]].to_string(index=False))

bad = df[df["source_surface_finite_fraction"] < 0.99]

if len(bad) > 0:
    raise RuntimeError("At least one PFSS output has finite fraction below 0.99.")

summary_table = table_dir / "pfss_hmi_metrics_summary_compact.csv"
df[[
    "cr",
    "rss",
    "nrho",
    "input_mean_abs_br",
    "source_surface_min_br",
    "source_surface_max_br",
    "source_surface_mean_abs_br",
    "source_surface_finite_fraction",
]].to_csv(summary_table, index=False)

print()
print(f"Saved compact table: {summary_table}")

fig_path = figdir / "pfss_hmi_source_surface_mean_abs_br_by_cr.png"

plt.figure(figsize=(8, 4))
plt.plot(df["cr"], df["source_surface_mean_abs_br"], marker="o")
plt.xlabel("Carrington rotation")
plt.ylabel("Mean unsigned source-surface Br")
plt.title("PFSS source-surface mean unsigned Br by Carrington rotation")
plt.tight_layout()
plt.savefig(fig_path, dpi=200)
plt.close()

print(f"Saved figure: {fig_path}")

print()
print("Status: pass")
