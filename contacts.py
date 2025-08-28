"""Helpers for drawing contact lines and highlighting residues.

This module provides utilities used by a visualization script to render
dashed contact lines between atoms and highlight interacting residues in a
uniform gray color.
"""

from __future__ import annotations

import numpy as np


# Default styling parameters -------------------------------------------------

# Maximum number of contacts to draw
N_LABEL = 10

# Cylinder radius depending on contact type
RAD_POLAR = 0.2
RAD_OTHER = 0.2

# Dashed line appearance
DASH_LEN = 0.90
GAP_LEN = 0.55
OPACITY = 0.75


def add_dashed_cylinders(
    v,
    p: dict,
    q: dict,
    color: str,
    radius: float,
    opacity: float = OPACITY,
    dash_len: float = DASH_LEN,
    gap_len: float = GAP_LEN,
):
    """Draw dashed cylinders between two points.

    Parameters
    ----------
    v : object
        py3Dmol view object.
    p, q : dict
        Dictionaries with ``x``, ``y`` and ``z`` coordinates for the start and
        end points of the dashed line.
    color : str
        Hex color string for the cylinders.
    radius : float
        Cylinder radius.
    opacity : float, optional
        Cylinder opacity.
    dash_len, gap_len : float, optional
        Length of the visible and invisible segments.
    """

    P = np.array([p['x'], p['y'], p['z']], float)
    Q = np.array([q['x'], q['y'], q['z']], float)
    V = Q - P
    L = float(np.linalg.norm(V))
    if L < 1e-6:
        return

    u = V / L
    n = max(3, int(L / max(0.1, (dash_len + gap_len))))
    for i in range(n):
        a = P + (dash_len + gap_len) * i * u
        b = a + dash_len * u
        t0 = max(0.0, min(L, (a - P) @ u))
        t1 = max(0.0, min(L, (b - P) @ u))
        if t1 <= t0:
            continue
        A = P + t0 * u
        B = P + t1 * u
        v.addCylinder(
            {
                'start': {'x': float(A[0]), 'y': float(A[1]), 'z': float(A[2])},
                'end': {'x': float(B[0]), 'y': float(B[1]), 'z': float(B[2])},
                'radius': radius,
                'color': color,
                'opacity': opacity,
            }
        )


def add_contacts(v, cs, line_color, res_color=None, nlabel=N_LABEL):
    """Add dashed lines for contacts and highlight residues in gray.

    Parameters
    ----------
    v : object
        py3Dmol view object.
    cs : list
        Contact definitions as produced by ``find_contacts``.
    line_color : str
        Color for the dashed contact lines.
    res_color : str, optional
        Ignored; kept for backwards compatibility with older scripts.
    nlabel : int, optional
        Maximum number of contacts to render.
    """

    shown = set()
    for c in cs[:nlabel]:
        aP, aL = c['prot_atom'], c['lig_atom']
        p = {
            'x': float(aP.coord[0]),
            'y': float(aP.coord[1]),
            'z': float(aP.coord[2]),
        }
        q = {
            'x': float(aL.coord[0]),
            'y': float(aL.coord[1]),
            'z': float(aL.coord[2]),
        }
        rad = RAD_POLAR if c['kind'] == 'polar' else RAD_OTHER
        add_dashed_cylinders(v, p, q, color=line_color, radius=rad)

        chain_id = aP.get_parent().get_parent().id
        resi = int(aP.get_parent().id[1])
        key = (chain_id, resi)
        if key not in shown:
            v.addStyle(
                {'model': 0, 'chain': chain_id, 'resi': resi},
                {'stick': {'radius': 0.24, 'color': '0x808080'}},
            )
            shown.add(key)
