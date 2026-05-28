from pathlib import Path

base = Path(".")
scripts = base / "scripts"

heights = [2.0, 2.5, 3.0]

def ftag(r):
    return f"rss{r:.1f}"

def stag(r):
    return f"rss{r:.1f}".replace(".", "p")

def must_replace(text, old, new, label):
    if old not in text:
        raise RuntimeError(f"Could not find expected text for {label}:\n{old}")
    return text.replace(old, new)

for r in heights:
    sf = stag(r)
    rf = ftag(r)

    t12 = (scripts / "12_run_pfss_hmi_selected_crs.py").read_text()
    t12 = must_replace(t12, "rss = 2.5", f"rss = {r:.1f}", "script 12 rss")
    (scripts / f"12_run_pfss_hmi_selected_crs_{sf}.py").write_text(t12)

    t14 = (scripts / "14_build_pfss_hmi_longitude_proxies.py").read_text()
    t14 = must_replace(t14, "r = 2.5", f"r = {r:.1f}", "script 14 r")
    t14 = must_replace(
        t14,
        'cp = td / "pfss_hmi_longitude_proxy_summary_compact.csv"',
        f'cp = td / "pfss_hmi_longitude_proxy_summary_compact_{rf}.csv"',
        "script 14 compact summary name",
    )
    (scripts / f"14_build_pfss_hmi_longitude_proxies_{sf}.py").write_text(t14)

    t17 = (scripts / "17_build_pfss_omni_ballistic_match.py").read_text()
    t17 = must_replace(t17, "rss = 2.5", f"rss = {r:.1f}", "script 17 rss")
    t17 = must_replace(
        t17,
        "pfss_hmi_selected_crs_rss2.5_longitude_profiles_all.csv",
        f"pfss_hmi_selected_crs_{rf}_longitude_profiles_all.csv",
        "script 17 PFSS longitude input",
    )
    t17 = must_replace(
        t17,
        'fd = b / "outputs" / "figures" / "comparison"',
        f'fd = b / "outputs" / "figures" / "comparison" / "{rf}"',
        "script 17 comparison figure folder",
    )
    t17 = must_replace(
        t17,
        'td = b / "outputs" / "tables"',
        f'td = b / "outputs" / "tables" / "{rf}"',
        "script 17 table folder",
    )

    needle = f'td = b / "outputs" / "tables" / "{rf}"\n'
    if needle in t17:
        t17 = t17.replace(
            needle,
            needle + "fd.mkdir(parents=True, exist_ok=True)\ntd.mkdir(parents=True, exist_ok=True)\n",
            1,
        )

    (scripts / f"17_build_pfss_omni_ballistic_match_{sf}.py").write_text(t17)

print("Created height-safe scripts:")
for r in heights:
    sf = stag(r)
    print(f"  scripts/12_run_pfss_hmi_selected_crs_{sf}.py")
    print(f"  scripts/14_build_pfss_hmi_longitude_proxies_{sf}.py")
    print(f"  scripts/17_build_pfss_omni_ballistic_match_{sf}.py")
