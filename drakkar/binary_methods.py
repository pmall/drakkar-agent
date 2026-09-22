"""Binary detection-method definition shared by dataset scripts.

A method counts as **binary** when it reports a direct 1:1 physical contact
between the two proteins, as opposed to co-complex / affinity capture,
proximity labelling, cross-linking, enzymatic-reaction or nucleic-acid methods.

``BINARY_PSIMI_IDS`` holds stable PSI-MI ontology terms (``methods.psimi_id``,
also exposed as ``dataset.psimi_id``). Never match on ``methods.id`` -- those
internal ids are not stable across database versions. The allowlist was built
from the PSI-MI detection-method branch, restricted to methods actually
present in valid vh rows.

Deliberately EXCLUDED as non-binary: mass spectrometry of complexes
(MI:0069), all co-immunoprecipitation / pull-down / affinity chromatography /
TAP, proximity labelling & BioID, proximity ligation assay, cross-linking
studies, EM, co-sedimentation / co-migration / gel filtration / native PAGE,
light / X-ray / neutron scattering, DLS, EPR, EMSA / mobility shift, antibody
array, western blot, protein three-hybrid, and enzymatic-reaction methods
(phosphorylation, ubiquitination, cleavage, ...).

Borderline calls (protein array, far-western, ELISA, HDX-MS, disulfide bond,
CD, DSC) are included.

An interaction is considered binary when at least one supporting description
used one of these methods -- see :func:`is_binary_psimi` and
:func:`format_binary_methods` for the reusable API.
"""

from __future__ import annotations

# The allowlist as data: (group label, ((PSI-MI id, method label), ...)).
# Labels mirror methods.name; the PSI-MI id is the stable key.
BINARY_METHOD_GROUPS: tuple[tuple[str, tuple[tuple[str, str], ...]], ...] = (
    (
        "two-hybrid and reporter-recruitment assays",
        (
            ("MI:0018", "two hybrid"),
            ("MI:0397", "two hybrid array"),
            ("MI:0398", "two hybrid pooling approach"),
            ("MI:1112", "two hybrid prey pooling approach"),
            ("MI:1356", "validated two hybrid"),
            ("MI:2215", "barcode fusion genetics two hybrid"),
            ("MI:0728", "gal4 vp16 complementation"),
            ("MI:0231", "mammalian protein protein interaction trap"),
            ("MI:0655", "lambda repressor two hybrid"),
        ),
    ),
    (
        "protein-fragment complementation",
        (
            ("MI:0090", "protein complementation assay"),
            ("MI:0809", "bimolecular fluorescence complementation"),
            ("MI:0112", "ubiquitin reconstruction"),
            ("MI:1203", "split luciferase complementation"),
            ("MI:1037", "split renilla luciferase complementation"),
            ("MI:2170", "bimolecular luminiscence complementation"),
            ("MI:2219", "gaussia luciferase protein tag"),
            ("MI:0014", "adenylate cyclase complementation"),
            ("MI:0011", "beta lactamase complementation"),
            ("MI:0010", "beta galactosidase complementation"),
            ("MI:0729", "luminescence based mammalian interactome mapping"),
            ("MI:0012", "bioluminescence resonance energy transfer"),
        ),
    ),
    (
        "resonance energy transfer / homogeneous proximity",
        (
            ("MI:0055", "fluorescent resonance energy transfer"),
            ("MI:0510", "homogeneous time resolved fluorescence"),
            (
                "MI:0905",
                "amplified luminescent proximity homogeneous assay",
            ),
        ),
    ),
    (
        "direct biophysical binding of purified partners",
        (
            ("MI:0107", "surface plasmon resonance"),
            ("MI:0921", "surface plasmon resonance array"),
            ("MI:0065", "isothermal titration calorimetry"),
            ("MI:0969", "bio-layer interferometry"),
            ("MI:1247", "microscale thermophoresis"),
            ("MI:0053", "fluorescence polarization spectroscopy"),
            ("MI:0052", "fluorescence correlation spectroscopy"),
            ("MI:0017", "classical fluorescence spectroscopy"),
            ("MI:0051", "fluorescence technology"),
            ("MI:1235", "thermal shift binding"),
            ("MI:1311", "differential scanning calorimetry"),
            ("MI:2196", "quartz crystal microbalance"),
            ("MI:0968", "biosensor"),
            ("MI:0049", "filter binding"),
            ("MI:0405", "competition binding"),
            ("MI:0892", "solid phase assay"),
            ("MI:0099", "scintillation proximity assay"),
            ("MI:0016", "circular dichroism"),
        ),
    ),
    (
        "structural methods implying atomic contact",
        (
            ("MI:0114", "x-ray crystallography"),
            ("MI:0077", "nuclear magnetic resonance"),
            ("MI:1103", "solution state nmr"),
            ("MI:0824", "x-ray powder diffraction"),
            ("MI:0944", "mass spectrometry study of hydrogen/deuterium exchange"),
            ("MI:0408", "disulfide bond"),
        ),
    ),
    (
        "immobilised-bait binding / display selection",
        (
            ("MI:0089", "protein array"),
            ("MI:0081", "peptide array"),
            ("MI:0084", "phage display"),
            ("MI:0048", "filamentous phage display"),
            ("MI:0066", "lambda phage display"),
            ("MI:0108", "t7 phage display"),
            ("MI:0115", "yeast display"),
            ("MI:0047", "far western blotting"),
            ("MI:0411", "enzyme linked immunosorbent assay"),
            ("MI:0695", "sandwich immunoassay"),
        ),
    ),
)

# Flat allowlist, derived from the groups above -- the single source of truth.
BINARY_PSIMI_IDS: tuple[str, ...] = tuple(
    psimi for _, members in BINARY_METHOD_GROUPS for psimi, _ in members
)

_BINARY_SET: frozenset[str] = frozenset(BINARY_PSIMI_IDS)


def is_binary_psimi(psimi_id: str) -> bool:
    """Return True when a PSI-MI detection-method term counts as binary."""
    return psimi_id in _BINARY_SET


def format_binary_methods() -> str:
    """Render the grouped allowlist as text, for reports and papers."""
    lines: list[str] = []
    for group, members in BINARY_METHOD_GROUPS:
        lines.append(f"{group}:")
        lines += [f"  {psimi} -- {label}" for psimi, label in members]
    return "\n".join(lines) + "\n"
