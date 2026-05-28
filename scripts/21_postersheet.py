from pathlib import Path
import math
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

root = Path("outputs/poster_figures")
out = root / "poster_figure_contact_sheet.png"

items = [
    ("Main scatter: PFSS eq abs Br vs OMNI |B|", root / "fig_main_scatter_ballistic_eqabs_bmag.png"),
    ("Null test: circular phase shifts", root / "fig_null_test_eqabs_bmag.png"),
    ("Ballistic lag by CR", root / "fig_ballistic_lag_by_cr.png"),
    ("PFSS proxy by phase and CR", root / "fig_pfss_eqabs_phase_by_cr.png"),
    ("OMNI |B| by phase and CR", root / "fig_omni_bmag_phase_by_cr.png"),
]

pics = [(lab, path, mpimg.imread(path)) for lab, path in items if path.exists()]

cols = 2
rows = math.ceil(len(pics) / cols)

fig, axs = plt.subplots(rows, cols, figsize=(14, 4.8 * rows))
axs = axs.ravel() if hasattr(axs, "ravel") else [axs]

for ax, (lab, path, img) in zip(axs, pics):
    ax.imshow(img)
    ax.set_title(lab, fontsize=12)
    ax.axis("off")

for ax in axs[len(pics):]:
    ax.axis("off")

fig.tight_layout()
fig.savefig(out, dpi=200)
print(f"saved {out}")
