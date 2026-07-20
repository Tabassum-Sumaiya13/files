"""
markers.py — Cross-cohort marker-name resolution.

WHY THIS EXISTS
---------------
Marker columns are named differently in every cohort, and the old validator
compared them with exact string equality (`"CD45" in marker_cols`). That silently
failed on the CRC cohort:

    panel actually contains : 'CD45 - hematopoietic cells:Cyc_4_ch_2'
    validator reported      : "0/10 canonical lineage markers present: []"

so the one expression-aware sanity check on the cell-type -> lineage map reported
"cannot verify" for a cohort that unquestionably has CD45, and the map went in
unchallenged. This module fixes that.

THE THREE PROBLEMS, AND WHICH ONE NEEDS A HUMAN
-----------------------------------------------
1. Decoration    'CD45 - hematopoietic cells:Cyc_4_ch_2'  -> 'CD45'
                 Mechanical: strip the CODEX channel tag after ':' and the
                 free-text description after ' - '. Generalises to any cohort.
2. Formatting    'FOXP3' / 'FoxP3',  'PD-L1' / 'PDL1'     -> 'foxp3' / 'pdl1'
                 Mechanical: casefold + drop non-alphanumerics. Generalises.
3. True aliases  'Cytokeratin' == 'PanCK',  'CD3' == 'CD3e'
                 NOT mechanical. No rule derives these; they are curated below
                 in SYNONYMS and must be extended by hand for a new panel.

Matching is on the WHOLE normalised token, never a substring, so 'CD4' can never
match 'CD45' and 'CD3' can never match 'CD31'.

WHEN A MARKER CANNOT BE RESOLVED
--------------------------------
`resolve_panel` returns it in `missing`, by canonical name. Callers report those
names so the gap is visible. A silent "0/10" is the failure mode this module was
written to remove.

VERSIONING
----------
Bump VOCABULARY_VERSION whenever SYNONYMS changes. Validation reports record it,
so a report can be traced to the vocabulary that produced it.
"""
from __future__ import annotations

import re
from typing import Dict, Iterable, List, Optional, Tuple

VOCABULARY_VERSION = "1.0.0"

# ---------------------------------------------------------------------------
# Curated cross-cohort aliases — canonical name -> other names for the SAME
# protein. This is the part that cannot be derived; every entry is a human
# judgement and carries its reason.
#
# To add a cohort whose panel uses a novel alias: add it here, bump
# VOCABULARY_VERSION. Do NOT add near-synonyms that are different proteins
# (e.g. CD3 vs CD3e are the same complex read the same way -> fine;
# CD4 vs CD40 are not -> never).
# ---------------------------------------------------------------------------
SYNONYMS: Dict[str, List[str]] = {
    # pan-cytokeratin: the epithelial/tumour anchor. UPMC calls it PanCK,
    # Schurch CRC calls it Cytokeratin.
    "PanCK": ["Cytokeratin", "panCK", "pan-cytokeratin", "PanCytokeratin", "CK"],
    # CD3 epsilon chain: UPMC names the chain, CRC names the complex.
    "CD3e": ["CD3", "CD3-epsilon", "CD3eps"],
    # checkpoint markers, hyphenation differs by vendor
    "PDL1": ["PD-L1", "PDL-1", "CD274"],
    "PD1": ["PD-1", "PDCD1", "CD279"],
    "HLA-DR": ["HLADR", "MHC-II", "MHCII"],
    "GranzymeB": ["Granzyme B", "GZMB", "GranzymeB"],
    "FoxP3": ["FOXP3", "Foxp3"],
    "CollagenIV": ["Collagen IV", "Col4", "COL4A1"],
    "aSMA": ["SMA", "alpha-SMA", "ACTA2", "SMActin"],
    "Podoplanin": ["PDPN", "D2-40"],
    "Vimentin": ["VIM"],
    "MUC1": ["MUC-1"],
    "CD45RO": ["CD45R0"],
    "Ki67": ["Ki-67", "MKI67"],
}


def normalise(name: str) -> str:
    """Reduce a raw column name to a comparable token.

    'CD45 - hematopoietic cells:Cyc_4_ch_2' -> 'cd45'
    'Granzyme B - cytotoxicity:Cyc_13_ch_2' -> 'granzymeb'
    'PD-L1 - checkpoint:Cyc_5_ch_3'         -> 'pdl1'
    'FoxP3'                                 -> 'foxp3'
    """
    s = str(name)
    s = s.split(":")[0]        # drop CODEX channel tag  ':Cyc_4_ch_2'
    s = s.split(" - ")[0]      # drop free-text description  ' - hematopoietic cells'
    return re.sub(r"[^a-z0-9]", "", s.strip().casefold())


def _alias_tokens(canonical: str) -> List[str]:
    """Every normalised token that should match `canonical`."""
    return [normalise(canonical)] + [normalise(a) for a in SYNONYMS.get(canonical, [])]


def resolve(canonical: str, available: Iterable[str]) -> Optional[str]:
    """Find the column in `available` that carries marker `canonical`.

    Returns the ACTUAL column name (so the caller can index the dataframe), or
    None. Exact normalised match first, then curated aliases — so a cohort that
    names the marker plainly never depends on the alias table.
    """
    cols = list(available)
    index: Dict[str, str] = {}
    for c in cols:
        index.setdefault(normalise(c), c)  # first occurrence wins, stable
    for token in _alias_tokens(canonical):
        if token in index:
            return index[token]
    return None


def resolve_panel(canonicals: Iterable[str],
                  available: Iterable[str]) -> Tuple[Dict[str, str], List[str]]:
    """Resolve a whole panel.

    Returns
    -------
    found   : {canonical_name: actual_column_name}
    missing : [canonical_name, ...]   — reported by name, never silently dropped
    """
    cols = list(available)
    found, missing = {}, []
    for m in canonicals:
        hit = resolve(m, cols)
        if hit is None:
            missing.append(m)
        else:
            found[m] = hit
    return found, missing
