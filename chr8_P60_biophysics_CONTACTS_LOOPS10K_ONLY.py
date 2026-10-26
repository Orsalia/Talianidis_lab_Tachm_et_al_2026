#!/usr/bin/env python3
"""
chr8_P60_CONTACTS_LOOPS10K_ONLY.py

Whole mouse chr8 coarse-grained chromatin model WITHOUT OpenMM.

Data roles
----------
--contacts    P60_chr8_matrix...       dense Hi-C / contact pairs
--loops10k    P60loops2mb...           10-kb loops

The chromosome is a continuous 20-kb polymer.  The script:
1) derives 100-kb A/B-like states from distance-normalized Hi-C O/E;
2) maps these states to 20-kb beads;
3) selects enriched long-range O/E contacts rather than treating every Hi-C
   observation as an equal spring;
4) maps 10-kb loops to the 20-kb polymer;
5) uses uniform polymer bending because TAD boundaries are not supplied;
6) uses excluded volume and weak same-compartment attraction;
7) adds weak generic short-range chromatin cohesion for finite-density packing;
8) uses Hi-C PC1-like states without external epigenomic orientation;
10) uses a soft anisotropic territory rather than spherical confinement;
11) relaxes the polymer using NumPy/SciPy forces.

This is a plausible data-constrained ensemble model, NOT a measured unique
single-cell chromosome conformation.

Dependencies: numpy, pandas, scipy
"""

import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree, ConvexHull, QhullError


def table(path):
    return pd.read_csv(path, sep=r"\s+", header=None, comment="#")


def pairs(path, chrom):
    d = table(path)
    if d.shape[1] < 6:
        raise ValueError(f"{path}: expected >=6 columns: chr s e chr s e [value]")
    return d[(d[0].astype(str) == chrom) & (d[3].astype(str) == chrom)].copy()


def intervals(path, chrom):
    d = table(path)
    if d.shape[1] < 3:
        raise ValueError(f"{path}: expected >=3 columns")
    return d[d[0].astype(str) == chrom].copy()


def bin_index(pos, start, res, n):
    return int(np.clip((int(pos) - start) // res, 0, n - 1))


def lamin_track(d, start, end, res, n):
    out = np.zeros(n, float)
    cov = np.zeros(n, float)
    for _, r in d.iterrows():
        a, b = max(start, int(r[1])), min(end, int(r[2]))
        if b <= a:
            continue
        score = 1.0
        if len(r) >= 4:
            z = pd.to_numeric(r[3], errors="coerce")
            score = float(z) if np.isfinite(z) else 0.0
        i0 = max(0, (a-start)//res)
        i1 = min(n-1, (b-1-start)//res)
        for i in range(i0, i1+1):
            bs, be = start+i*res, min(end, start+(i+1)*res)
            ov = max(0, min(b, be)-max(a, bs))
            out[i] += score*ov
            cov[i] += ov
    out = np.divide(out, cov, out=np.zeros_like(out), where=cov > 0)
    positive = out[out > 0]
    if len(positive):
        scale = np.percentile(positive, 95)
        if scale > 0:
            out = np.clip(out/scale, 0, 1)
    return out


def h3k27ac_track(d, start, end, res, n):
    """Map H3K27ac BED/narrowPeak signal to model beads and scale to 0..1.

    narrowPeak: column 7 (0-based index 6, signalValue) is preferred.
    BED/bedGraph: column 4 (0-based index 3) is used when numeric.
    If neither is available, every peak contributes weight 1.
    """
    out = np.zeros(n, float)
    cov = np.zeros(n, float)
    for _, r in d.iterrows():
        a, b = max(start, int(r[1])), min(end, int(r[2]))
        if b <= a:
            continue
        score = 1.0
        if len(r) >= 7:
            z = pd.to_numeric(r.iloc[6], errors="coerce")
            if np.isfinite(z) and float(z) >= 0:
                score = float(z)
        elif len(r) >= 4:
            z = pd.to_numeric(r.iloc[3], errors="coerce")
            if np.isfinite(z) and float(z) >= 0:
                score = float(z)
        i0 = max(0, (a-start)//res)
        i1 = min(n-1, (b-1-start)//res)
        for i in range(i0, i1+1):
            bs, be = start+i*res, min(end, start+(i+1)*res)
            ov = max(0, min(b, be)-max(a, bs))
            out[i] += score*ov
            cov[i] += ov
    out = np.divide(out, cov, out=np.zeros_like(out), where=cov > 0)
    positive = out[out > 0]
    if len(positive):
        scale = np.percentile(positive, 95)
        if scale > 0:
            out = np.clip(out/scale, 0, 1)
    return out


def structural_stats(X, excluded):
    """Return global structural statistics for the final polymer."""
    center = X.mean(axis=0)
    Y = X-center
    rg = float(np.sqrt(np.mean(np.sum(Y*Y, axis=1))))
    end_to_end = float(np.linalg.norm(X[-1]-X[0]))
    span = np.ptp(X, axis=0)
    try:
        hull_volume = float(ConvexHull(X).volume)
    except (QhullError, ValueError):
        hull_volume = float("nan")
    compactness = float(len(X)/hull_volume) if np.isfinite(hull_volume) and hull_volume > 0 else float("nan")
    bead_radius = excluded/2.0
    bead_volume = 4.0*np.pi*(bead_radius**3)/3.0
    packing_fraction = float(len(X)*bead_volume/hull_volume) if np.isfinite(hull_volume) and hull_volume > 0 else float("nan")
    bonds = np.linalg.norm(np.diff(X, axis=0), axis=1)
    return {
        "beads": len(X),
        "radius_of_gyration": rg,
        "end_to_end_distance": end_to_end,
        "convex_hull_volume": hull_volume,
        "compactness_beads_per_hull_volume": compactness,
        "approx_packing_fraction": packing_fraction,
        "span_x": float(span[0]), "span_y": float(span[1]), "span_z": float(span[2]),
        "median_bond": float(np.median(bonds)), "max_bond": float(bonds.max())
    }


def boundary_bins(d, start, end, res, n):
    ans = set()
    for _, r in d.iterrows():
        x = (int(r[1])+int(r[2]))//2
        if start <= x < end:
            ans.add(bin_index(x, start, res, n))
    return ans


def hic_oe(d, start, end, cres):
    """Aggregate contacts to compartment resolution and distance-normalize."""
    n = int(np.ceil((end-start)/cres))
    C = np.zeros((n, n), np.float64)

    x = pd.to_numeric(d[1], errors="coerce").to_numpy()
    y = pd.to_numeric(d[4], errors="coerce").to_numpy()
    good = np.isfinite(x) & np.isfinite(y)
    x, y = x[good], y[good]
    i = ((x-start)//cres).astype(int)
    j = ((y-start)//cres).astype(int)
    good = (i >= 0) & (i < n) & (j >= 0) & (j < n)
    i, j = i[good], j[good]

    # If a seventh column exists and is numeric, use it as contact weight.
    if d.shape[1] >= 7:
        raw = pd.to_numeric(d.iloc[:, 6], errors="coerce").fillna(0).to_numpy()
        raw = raw[np.where(np.isfinite(pd.to_numeric(d[1], errors="coerce").to_numpy()) &
                           np.isfinite(pd.to_numeric(d[4], errors="coerce").to_numpy()))[0]]
        w = np.maximum(raw[good], 0)
    else:
        w = np.ones(len(i), float)

    np.add.at(C, (i, j), w)
    off = i != j
    np.add.at(C, (j[off], i[off]), w[off])

    OE = np.zeros_like(C)
    for sep in range(n):
        v = np.diag(C, k=sep)
        nz = v[v > 0]
        if not len(nz):
            continue
        expected = np.mean(nz)
        ii = np.arange(n-sep)
        jj = ii+sep
        OE[ii, jj] = C[ii, jj]/expected
        OE[jj, ii] = OE[ii, jj]
    np.fill_diagonal(OE, 1.0)
    return C, OE


def compartment_pc1(OE, lamin100=None, h3_100=None):
    """PC1 from Hi-C O/E; orient positive as A-like using H3K27ac and/or Lamin."""
    X = np.log1p(OE)
    good = np.std(X, axis=1) > 1e-12
    R = np.eye(len(X))
    if good.sum() >= 3:
        rg = np.corrcoef(X[good])
        R[np.ix_(good, good)] = np.nan_to_num(rg)
    vals, vecs = np.linalg.eigh(R)
    pc = vecs[:, -1].copy()

    # PC sign is arbitrary. Positive is oriented toward active/H3K27ac-rich
    # chromatin and away from Lamin-rich chromatin.
    orient = 0.0
    pos = pc >= 0
    neg = ~pos
    if h3_100 is not None and np.any(h3_100 > 0) and pos.any() and neg.any():
        orient += float(np.mean(h3_100[pos]) - np.mean(h3_100[neg]))
    if lamin100 is not None and np.any(lamin100 > 0) and pos.any() and neg.any():
        orient += float(np.mean(lamin100[neg]) - np.mean(lamin100[pos]))
    if orient < 0:
        pc = -pc

    s = np.percentile(np.abs(pc), 95)
    if s > 0:
        pc = np.clip(pc/s, -1, 1)
    return pc

def map_pc(pc, start, res, cres, n):
    centers = start+(np.arange(n)+0.5)*res
    ci = np.clip(((centers-start)//cres).astype(int), 0, len(pc)-1)
    bpc = pc[ci]
    state = np.where(bpc >= 0, 1, -1)
    return bpc, state


def select_oe(OE, start, cres, res, n, threshold, max_pairs, seed):
    out = []
    for i in range(len(OE)):
        for j in np.where(OE[i, i+3:] >= threshold)[0]+i+3:
            bi = bin_index(start+(i+.5)*cres, start, res, n)
            bj = bin_index(start+(j+.5)*cres, start, res, n)
            if bi != bj:
                out.append((bi, bj, float(OE[i, j])))
    if len(out) > max_pairs:
        rng = np.random.default_rng(seed)
        weight = np.array([max(z-1, 1e-6) for _, _, z in out])
        weight /= weight.sum()
        idx = rng.choice(len(out), max_pairs, replace=False, p=weight)
        out = [out[k] for k in idx]
    return out


def mapped_loops(d, start, res, n):
    ans = set()
    for _, r in d.iterrows():
        i = bin_index(r[1], start, res, n)
        j = bin_index(r[4], start, res, n)
        if abs(i-j) < 2:
            continue
        if i > j:
            i, j = j, i
        ans.add((i, j))
    return list(ans)


def initial_chain(n, bond, axes, seed):
    rng = np.random.default_rng(seed)
    X = np.zeros((n, 3))
    v = rng.normal(size=3)
    v /= np.linalg.norm(v)
    for i in range(1, n):
        z = rng.normal(size=3)
        z /= np.linalg.norm(z)
        v = .72*v + .28*z
        v /= np.linalg.norm(v)
        X[i] = X[i-1] + bond*v
    X -= X.mean(0)
    ax, ay, az = axes
    q = np.sqrt((X[:,0]/ax)**2+(X[:,1]/ay)**2+(X[:,2]/az)**2)
    X *= .65/(np.percentile(q, 95)+1e-12)
    return X


def spring(F, X, i, j, k, target):
    dv = X[j]-X[i]
    r = np.linalg.norm(dv)+1e-12
    # one-sided attraction: do not compress already-close beads
    if r > target:
        f = k*(r-target)*dv/r
        F[i] += f
        F[j] -= f


def calculate_forces(X, state, lamin, h3k27ac, boundaries, oe, loops10, loops1, a):
    n = len(X)
    F = np.zeros_like(X)

    # 1. Polymer connectivity
    dv = X[1:]-X[:-1]
    rr = np.linalg.norm(dv, axis=1)+1e-12
    fb = a.k_bond*(rr-a.bond)[:,None]*dv/rr[:,None]
    F[:-1] += fb
    F[1:] -= fb

    # 2. Semiflexible chain; TAD boundaries locally reduce bending penalty.
    lap = 2*X[1:-1]-X[:-2]-X[2:]
    kb = np.full(n-2, a.k_bend)
    for b in boundaries:
        if 1 <= b < n-1:
            kb[b-1] *= a.boundary_flex
    F[1:-1] -= kb[:,None]*lap

    # 3. Hi-C O/E, 10-kb loops, 1-kb loops
    for i,j,w in oe:
        spring(F, X, i, j, a.k_oe*min(max(w,1),3), a.oe_target)
    for i,j in loops10:
        spring(F, X, i, j, a.k_loop10, a.loop10_target)
    for i,j in loops1:
        spring(F, X, i, j, a.k_loop1, a.loop1_target)

    # 4. Local nonbonded interactions: excluded volume + compartments.
    tree = cKDTree(X)
    pp = tree.query_pairs(a.comp_cutoff, output_type="ndarray")
    if len(pp) > a.max_nonbonded:
        stride = int(np.ceil(len(pp)/a.max_nonbonded))
        pp = pp[::stride]
    for i,j in pp:
        i,j = int(i),int(j)
        if abs(i-j) <= a.exclude_neighbors:
            continue
        dv = X[j]-X[i]
        r = np.linalg.norm(dv)+1e-12
        u = dv/r
        if r < a.excluded:
            f = a.k_excluded*(a.excluded-r)
            F[i] -= f*u
            F[j] += f*u
        else:
            # Weak generic chromatin attraction keeps the territory continuously
            # occupied instead of forming compact nodes joined by long strings.
            if r < a.chromatin_range:
                activity = 0.5*(h3k27ac[i]+h3k27ac[j])
                f0 = a.k_chromatin*(1.0-a.h3_decompact*activity)*(a.chromatin_range-r)
                F[i] += f0*u
                F[j] -= f0*u

            # A/B-like microphase organization is deliberately weaker than the
            # generic packing term: same-state attraction, very weak cross-state
            # repulsion. This creates density domains without isolated globules.
            if r < a.comp_range:
                if state[i] == state[j]:
                    f = a.k_comp*(a.comp_range-r)
                    F[i] += f*u
                    F[j] -= f*u
                else:
                    f = a.k_cross*(a.comp_range-r)
                    F[i] -= f*u
                    F[j] += f*u

    # 5. Very soft anisotropic territory.
    ax, ay, az = a.axes
    q = np.sqrt((X[:,0]/ax)**2+(X[:,1]/ay)**2+(X[:,2]/az)**2)+1e-12
    grad = np.column_stack((X[:,0]/ax**2, X[:,1]/ay**2, X[:,2]/az**2))/q[:,None]
    outside = q > 1
    F[outside] -= a.k_conf*(q[outside]-1)[:,None]*grad[outside]

    # 6. H3K27ac marks active chromatin. Give acetylated beads a weak interior
    # preference and reduce their generic cohesion (above), producing mild
    # local decompaction without breaking polymer continuity.
    F += (-a.k_h3k27ac*h3k27ac*q)[:,None]*grad

    mag = np.linalg.norm(F, axis=1)
    big = mag > a.force_cap
    F[big] *= (a.force_cap/mag[big])[:,None]
    return F


def write_pdb(path, X, lamin, state, chrom, start, end, res):
    with open(path, "w") as f:
        f.write(f"REMARK {chrom}:{start}-{end}; {res} bp per bead\n")
        f.write("REMARK occupancy: A-like=1 B-like=0; B-factor: Lamin x100\n")
        # PDB residue number is only 4 columns, so wrap residue IDs safely.
        for k, ((x,y,z), L, s) in enumerate(zip(X,lamin,state), 1):
            resid = ((k-1) % 9999)+1
            occ = 1.0 if s > 0 else 0.0
            f.write(f"ATOM  {k:5d}  CA  GLY A{resid:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{100*L:6.2f}           C\n")
        for k in range(1, len(X)):
            f.write(f"CONECT{k:5d}{k+1:5d}\n")
        f.write("END\n")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--contacts", required=True)
    p.add_argument("--loops10k", required=True)
    p.add_argument("--chrom", default="chr8")
    p.add_argument("--start", type=int, default=0)
    p.add_argument("--end", type=int, default=129401213)
    p.add_argument("--res", type=int, default=20000)
    p.add_argument("--comp-res", type=int, default=100000)
    p.add_argument("--out", default="chr8_P60_basic_biophysics")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--steps", type=int, default=12000)
    p.add_argument("--dt", type=float, default=.002)
    p.add_argument("--noise", type=float, default=.015)

    p.add_argument("--bond", type=float, default=1.0)
    p.add_argument("--k-bond", type=float, default=18)
    p.add_argument("--k-bend", type=float, default=.48)
    p.add_argument("--boundary-flex", type=float, default=.25)
    p.add_argument("--excluded", type=float, default=.72)
    p.add_argument("--k-excluded", type=float, default=5)
    p.add_argument("--k-chromatin", type=float, default=.018,
                   help="weak generic short-range chromatin attraction")
    p.add_argument("--chromatin-range", type=float, default=2.05)

    p.add_argument("--oe-threshold", type=float, default=2.0)
    p.add_argument("--max-oe", type=int, default=6000)
    p.add_argument("--k-oe", type=float, default=.08)
    p.add_argument("--oe-target", type=float, default=2.25)
    p.add_argument("--k-loop10", type=float, default=.40)
    p.add_argument("--loop10-target", type=float, default=1.55)
    p.add_argument("--k-loop1", type=float, default=.65)
    p.add_argument("--loop1-target", type=float, default=1.30)

    p.add_argument("--k-comp", type=float, default=.022)
    p.add_argument("--k-cross", type=float, default=.010)
    p.add_argument("--comp-range", type=float, default=2.35)
    p.add_argument("--comp-cutoff", type=float, default=2.60)
    p.add_argument("--exclude-neighbors", type=int, default=5)
    p.add_argument("--max-nonbonded", type=int, default=120000)

    p.add_argument("--axes", type=float, nargs=3, default=[30,24,20],
                   metavar=("X","Y","Z"))
    p.add_argument("--k-conf", type=float, default=.38)
    p.add_argument("--lamin-shell", type=float, default=.88)
    p.add_argument("--k-h3k27ac", type=float, default=.020, help="weak active-chromatin interior bias")
    p.add_argument("--h3-decompact", type=float, default=.45, help="fractional reduction of generic cohesion at H3K27ac-rich beads")
    p.add_argument("--force-cap", type=float, default=20)
    a = p.parse_args()

    n = int(np.ceil((a.end-a.start)/a.res))
    H = pairs(a.contacts, a.chrom)
    L10 = pairs(a.loops10k, a.chrom)
    # Only dense Hi-C contacts and 10-kb loops are used.
    lamin = np.zeros(n, float)
    h3k27ac = np.zeros(n, float)
    boundaries = set()

    print("Calculating distance-normalized Hi-C O/E...")
    C, OE = hic_oe(H, a.start, a.end, a.comp_res)
    lam100 = np.zeros(len(OE), float)
    h3_100 = np.zeros(len(OE), float)
    # Without H3K27ac/Lamin the PC1 sign is arbitrary; states are Hi-C PC1-like.
    pc = compartment_pc1(OE, None, None)
    bead_pc, state = map_pc(pc, a.start, a.res, a.comp_res, n)

    oe = select_oe(OE, a.start, a.comp_res, a.res, n,
                   a.oe_threshold, a.max_oe, a.seed)
    loops10 = mapped_loops(L10, a.start, a.res, n)
    loops1 = []

    print(f"Whole {a.chrom}: {n:,} beads at {a.res:,} bp")
    print(f"Hi-C rows: {len(H):,}; selected O/E restraints: {len(oe):,}")
    print(f"10-kb loops: {len(loops10):,}")
    print("1-kb loops: not used")
    print("TAD boundaries: not used")
    print("H3K27ac/Lamin: not used")
    print(f"A-like beads: {(state>0).sum():,}; B-like beads: {(state<0).sum():,}")

    X = initial_chain(n, a.bond, a.axes, a.seed)
    rng = np.random.default_rng(a.seed)

    for step in range(1, a.steps+1):
        F = calculate_forces(X, state, lamin, h3k27ac, boundaries, oe, loops10, loops1, a)
        X += a.dt*F + a.noise*np.sqrt(a.dt)*rng.normal(size=X.shape)
        X -= X.mean(axis=0)
        if step % 500 == 0 or step == 1:
            bonds = np.linalg.norm(np.diff(X,axis=0),axis=1)
            rg = np.sqrt(np.mean(np.sum(X*X,axis=1)))
            print(f"step {step:,}/{a.steps:,}: Rg={rg:.3f}; "
                  f"median bond={np.median(bonds):.3f}; max bond={bonds.max():.3f}")

    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    pdb = out.with_suffix(".pdb")
    write_pdb(pdb, X, lamin, state, a.chrom, a.start, a.end, a.res)

    starts = a.start+np.arange(n,dtype=np.int64)*a.res
    ends = np.minimum(starts+a.res, a.end)
    pd.DataFrame({
        "bin_index":np.arange(n),"chrom":a.chrom,"start":starts,"end":ends,
        "x":X[:,0],"y":X[:,1],"z":X[:,2],"PC1":bead_pc,
        "compartment":np.where(state>0,"A-like","B-like"),
        "lamin":lamin,"H3K27ac":h3k27ac,
        "TAD_boundary":[int(i in boundaries) for i in range(n)]
    }).to_csv(out.with_name(out.name+"_coords.csv"), index=False)

    cs = a.start+np.arange(len(pc),dtype=np.int64)*a.comp_res
    pd.DataFrame({
        "chrom":a.chrom,"start":cs,"end":np.minimum(cs+a.comp_res,a.end),
        "PC1_oriented":pc,
        "state":np.where(pc>=0,"A-like","B-like"),
        "lamin":lam100,"H3K27ac":h3_100
    }).to_csv(out.with_name(out.name+"_compartments.csv"), index=False)

    stats = structural_stats(X, a.excluded)
    stats.update({
        "OE_restraints":len(oe), "loops_10kb":len(loops10), "loops_1kb":len(loops1),
        "TAD_boundaries":len(boundaries),
        "H3K27ac_positive_beads":int((h3k27ac>0).sum())
    })
    pd.DataFrame([stats]).to_csv(out.with_name(out.name+"_QC.csv"),index=False)

    print("\nFinal structural statistics")
    print(f"  Radius of gyration: {stats['radius_of_gyration']:.4f}")
    print(f"  End-to-end distance: {stats['end_to_end_distance']:.4f}")
    print(f"  Convex hull volume: {stats['convex_hull_volume']:.4f}")
    print(f"  Compactness (beads/hull volume): {stats['compactness_beads_per_hull_volume']:.6g}")
    print(f"  Approx. packing fraction: {stats['approx_packing_fraction']:.6g}")

    print("\nFinished.")
    print("PDB:", pdb)
    print("Coordinates:", out.with_name(out.name+"_coords.csv"))
    print("Compartments:", out.with_name(out.name+"_compartments.csv"))
    print("QC:", out.with_name(out.name+"_QC.csv"))


if __name__ == "__main__":
    main()
