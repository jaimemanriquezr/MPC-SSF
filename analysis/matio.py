"""Read MPC-SSF .mat results from Python, in either MATLAB format.

Two formats are on disk and both must keep working:

  -v7   (138 files, written before 2026-08-25) -- readable by scipy.io.loadmat,
        contains `rec` (derived metrics) and `snap` only.
  -v7.3 (HDF5, written after) -- NOT readable by scipy at all. Contains the full
        `results` object, a flattened `results_py`, plus `rec` and `snap`.

`results` is a MATLAB object: it serialises as MCOS references that neither
scipy nor h5py can follow, so it is MATLAB-only by construction. `results_py`
exists for this module.

Use load_rec() for the derived metrics (works on every file) and load_py() for
the full species profiles (v7.3 only).
"""
import numpy as np

__all__ = ["mat_version", "load_rec", "load_py", "has_full"]


def mat_version(path):
    """'v7.3' or 'v7', from the file header."""
    with open(path, "rb") as fh:
        return "v7.3" if b"MATLAB 7.3" in fh.read(128) else "v7"


def _h5_str(ds):
    """MATLAB char -> str. Stored as uint16 code points."""
    a = np.asarray(ds).ravel()
    return "".join(chr(int(c)) for c in a)


def _h5_val(g, key):
    """One field out of an HDF5 struct group, with MATLAB's ordering undone."""
    ds = g[key]
    a = np.array(ds)
    # MATLAB writes column-major; h5py hands back the reversed shape. Undoing it
    # with .T restores the MATLAB orientation for 2-D and higher.
    if a.ndim >= 2:
        a = a.T
    return a


def load_rec(path):
    """The `rec` metrics struct as a plain dict. Works on v7 and v7.3."""
    if mat_version(path) == "v7":
        import scipy.io as sio
        d = sio.loadmat(path, squeeze_me=True, struct_as_record=False)["rec"]
        return {k: getattr(d, k) for k in d._fieldnames}
    import h5py
    with h5py.File(path, "r") as f:
        g = f["rec"]
        out = {}
        for k in g.keys():
            try:
                v = _h5_val(g, k)
            except Exception:
                continue
            # char fields come back as uint16; MATLAB marks them, but the simple
            # test is a small integer array that decodes to printable text.
            if v.dtype == np.uint16 and v.size and v.max() < 0x3000:
                out[k] = _h5_str(v)
            else:
                out[k] = v.squeeze() if v.size > 1 else (v.ravel()[0] if v.size else v)
        return out


def has_full(path):
    """True if this file carries the flattened results_py block."""
    if mat_version(path) == "v7":
        return False
    import h5py
    with h5py.File(path, "r") as f:
        return "results_py" in f


def load_py(path):
    """The flattened `results_py` block: every species profile, plus provenance.

    Returns a dict with
        time (nt,), z (nz,), porosity (nz,), dz, n0, delta, flag
        particles (nz, nt, nSpecies, 3)  phases Matrix|Enclosed|Flowing
        liquids   (nz, nt, nSpecies, 2)  phases Enclosed|Flowing
        water     (nz, nt)
        particleNames, liquidNames, phaseNames : list[str]
        opt : dict of the SolverOptions provenance
    """
    if not has_full(path):
        raise ValueError(f"{path} is {mat_version(path)} with no results_py; "
                         "only load_rec() is available for it")
    import h5py
    with h5py.File(path, "r") as f:
        P = f["results_py"]
        out = {}
        for k in ("time", "z", "porosity", "water", "particles", "liquids",
                  "velocityBiofilm"):
            if k in P:
                out[k] = _h5_val(P, k).squeeze()
        for k in ("dz", "n0", "delta"):
            if k in P:
                out[k] = float(np.array(P[k]).ravel()[0])
        for k in ("flag", "particleNames", "liquidNames", "phaseNames"):
            if k in P:
                out[k] = _h5_str(P[k])
        for k in ("particleNames", "liquidNames", "phaseNames"):
            if k in out:
                out[k] = out[k].split("|")
        out["opt"] = {}
        for k in P.keys():
            if not k.startswith("opt_"):
                continue
            v = np.array(P[k])
            name = k[4:]
            if v.dtype == np.uint16:
                out["opt"][name] = _h5_str(v)
            else:
                out["opt"][name] = float(v.ravel()[0])
        return out


def species(full, name):
    """Pull one species out of a load_py() dict as {phase: (nz, nt)}."""
    if name in full["liquidNames"]:
        i = full["liquidNames"].index(name)
        a = full["liquids"][:, :, i, :]
        return {"Enclosed": a[..., 0], "Flowing": a[..., 1]}
    i = full["particleNames"].index(name)
    a = full["particles"][:, :, i, :]
    return {"Matrix": a[..., 0], "Enclosed": a[..., 1], "Flowing": a[..., 2]}
