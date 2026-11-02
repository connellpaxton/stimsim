"""Validate an SWC tree without changing its geometry."""

import numpy as np


def load_swc(path, require_positive_radii=True):
    """Return SWC rows, parent row indices, and a parent-first traversal."""
    rows = np.loadtxt(path, comments='#', ndmin=2)
    if rows.shape[1] != 7 or not np.isfinite(rows).all():
        raise ValueError(f'{path}: expected seven finite columns')
    for col in (0, 1, 6):
        if not np.equal(rows[:, col], np.floor(rows[:, col])).all():
            raise ValueError(f'{path}: non-integer identifier or type')
    ids = rows[:, 0].astype(np.int64)
    lookup = {int(node): i for i, node in enumerate(ids)}
    if len(lookup) != len(ids):
        raise ValueError(f'{path}: duplicate node IDs')
    if require_positive_radii and (rows[:, 5] <= 0).any():
        raise ValueError(f'{path}: non-positive radii')
    parent = np.array([-1 if p == -1 else lookup.get(int(p), -2) for p in rows[:, 6]])
    roots = np.flatnonzero(parent == -1)
    if len(roots) != 1 or (parent == -2).any():
        raise ValueError(f'{path}: requires one root and valid parents')
    children = [[] for _ in ids]
    for i, p in enumerate(parent):
        if p >= 0:
            children[p].append(i)
    order = [int(roots[0])]
    for i in order:
        order.extend(children[i])
    if len(order) != len(ids):
        raise ValueError(f'{path}: cycle or disconnected tree')
    return rows, parent, np.asarray(order)
