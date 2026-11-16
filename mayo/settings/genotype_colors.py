"""
ONE place for genotype names, colors, plotting order and the reference genotype.

Used by the main notebook, cluster_analysis.py, cluster_coordination.py,
cluster_association_chi2.py and the plotting functions in mayo/base_functions.py.

To ADD A NEW GENOTYPE:
    1. Add its samples to mayo/settings/sample_genotypes.py   ("SAMPLE": "genotype name")
    2. Add one line below in GENOTYPE_COLORS                 ("genotype name": "color")
   That's it: it is then included in every "all genotypes" list and gets the same
   color in every plot. Position in GENOTYPE_COLORS = position in plots.

To CHANGE A COLOR: edit it below; every plot picks it up.
   Colors can be matplotlib names ("tab:blue", "gold") or hex codes ("#1f77b4").

To LEAVE A GENOTYPE OUT of the main analysis but keep its color, add it to
EXCLUDE_FROM_ANALYSIS.
"""

# genotype name  ->  color        (order here = order in plots)
GENOTYPE_COLORS: dict[str, str] = {
    "ung1∆":             "tab:blue",
    "ung1∆fun30∆":       "tab:olive",
    "ung1∆mph1∆":        "tab:brown",
    "ung1∆msh2∆":        "tab:pink",
    "ung1∆msh3∆":        "tab:cyan",
    "ung1∆msh6∆":        "tab:orange",
    "ung1∆mlh1∆":        "tab:gray",
    "ung1∆rnh201∆":      "tab:cyan",
    "ung1∆rnh1∆":        "tab:red",
    "ung1∆rnh1∆rnh201∆": "tab:green",   # was missing from cluster_analysis.py; change if you like
    "ung1∆ 37C_A3A":     "tab:orange",
    "ung1∆ 37C_EV":      "tab:orange",
    "ung1∆mlh2∆":        "tab:green",
    "ung1∆SNPless":      "tab:purple",
    "ung1∆SNPlessmsh2∆": "tab:green",
    "ung1∆pol3_01":      "tab:olive",
    "ung1∆pol3_5DV":     "tab:brown",
}

# Genotypes that have a color but are NOT part of the standard "all genotypes" analysis.
EXCLUDE_FROM_ANALYSIS: list[str] = [
    "ung1∆ 37C_A3A",
    "ung1∆ 37C_EV",
]

# Genotype every other genotype is compared against in statistical tests.
# Scripts can override this per run (see REFERENCE_GENOTYPE at the top of each script).
REFERENCE_GENOTYPE: str = "ung1∆"

# Color for anything not listed above (a warning is logged when it is used).
FALLBACK_COLOR: str = "lightgray"


# ---------------------------------------------------------------------------
# Helpers (no need to edit below)
# ---------------------------------------------------------------------------

# All genotypes used in the main analysis, in plotting order.
ANALYSIS_GENOTYPES: list[str] = [g for g in GENOTYPE_COLORS if g not in EXCLUDE_FROM_ANALYSIS]

# Suffixes used in some files (e.g. meiotic_frequencies_mayomate.csv: "ung1∆msh2∆ A3A")
_KNOWN_SUFFIXES = (" A3A", " EV")


def genotype_color(genotype: str) -> str:
    """Color for one genotype. Accepts names with ' A3A' / ' EV' suffixes and line breaks."""
    import logging
    name = str(genotype).replace("\n", "")
    if name in GENOTYPE_COLORS:
        return GENOTYPE_COLORS[name]
    for suffix in _KNOWN_SUFFIXES:
        if name.endswith(suffix) and name[: -len(suffix)] in GENOTYPE_COLORS:
            return GENOTYPE_COLORS[name[: -len(suffix)]]
    logging.warning(f"No color set for genotype '{genotype}' in mayo/settings/genotype_colors.py; using {FALLBACK_COLOR}")
    return FALLBACK_COLOR


def genotype_palette(genotypes=None) -> dict:
    """
    Palette dict for seaborn (palette=...). With no argument, returns all known genotypes.
    With a list of genotypes (or any iterable), returns colors for exactly those, so seaborn
    never fails with 'palette dictionary is missing keys'.
    """
    if genotypes is None:
        return dict(GENOTYPE_COLORS)
    return {g: genotype_color(g) for g in dict.fromkeys(genotypes)}


def genotype_order(genotypes) -> list:
    """Sort genotypes into the plotting order defined by GENOTYPE_COLORS (unknown ones go last)."""
    rank = {g: i for i, g in enumerate(GENOTYPE_COLORS)}
    return sorted(dict.fromkeys(genotypes), key=lambda g: rank.get(g, len(rank)))


def genotype_query(genotypes=None) -> str:
    """'ung1∆+ung1∆fun30∆+...' string for get_samples_with_genotype(); defaults to ANALYSIS_GENOTYPES."""
    return "+".join(ANALYSIS_GENOTYPES if genotypes is None else genotypes)
