#!/usr/bin/env python3
"""
Generate an interactive 3D GenomeFlow HTML from two PDB files.

Panels:
1. HNF4A 3D molecule
2. P60 3D molecule
3. 3D alignment: HNF4A aligned onto P60 + displacement toggle
4. Comparison statistics

Requirements:
    pip install numpy pandas scipy

Usage:
    python make_interactive_genomeflow_3d.py \
      --hnf4a chr8_104860000_106860000_HNF4A_GenomeFlow_REGEN.pdb \
      --p60 chr8_104860000_106860000_P60_GenomeFlow_REGEN.pdb \
      --chrom chr8 \
      --start 104860000 \
      --end 106860000 \
      --out-prefix chr8_HNF4A_P60_interactive
"""

import argparse
import json
import math
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd


def read_coords(path: Path) -> np.ndarray:
    pts = []
    for line in path.read_text(errors="ignore").splitlines():
        if line.startswith(("ATOM", "HETATM")):
            parsed = False
            try:
                pts.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                parsed = True
            except Exception:
                pass

            if not parsed:
                parts = line.split()
                for idx in [(6, 7, 8), (5, 6, 7)]:
                    try:
                        pts.append([float(parts[idx[0]]), float(parts[idx[1]]), float(parts[idx[2]])])
                        break
                    except Exception:
                        continue
    return np.asarray(pts, float)


def pairwise(X: np.ndarray) -> np.ndarray:
    G = X @ X.T
    d = np.diag(G)
    D2 = d[:, None] + d[None, :] - 2 * G
    D2[D2 < 0] = 0
    return np.sqrt(D2)


def normalize(X: np.ndarray) -> np.ndarray:
    D = pairwise(X)
    vals = D[np.triu_indices(len(X), 1)]
    vals = vals[vals > 0]
    med = np.median(vals) if vals.size else 1.0
    return X / med if med > 0 else X


def contact_bfactor(X: np.ndarray) -> np.ndarray:
    center = X.mean(axis=0)
    radius = np.linalg.norm(X - center, axis=1)
    b = 100 * (1 - (radius - radius.min()) / (np.ptp(radius) + 1e-12))
    if len(b) > 7:
        b = np.convolve(b, np.ones(7) / 7, mode="same")
        b = 100 * (b - b.min()) / (np.ptp(b) + 1e-12)
    return b


def kabsch(P: np.ndarray, Q: np.ndarray) -> np.ndarray:
    """Rigidly align P onto Q."""
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    U, S, Vt = np.linalg.svd(Pc.T @ Qc)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    return Pc @ R + Q.mean(axis=0)


def write_pdb(path: Path, X: np.ndarray, B: np.ndarray | None = None, chain: str = "A") -> None:
    if B is None:
        B = np.zeros(len(X))
    with open(path, "w") as f:
        for i, ((x, y, z), b) in enumerate(zip(X, B), start=1):
            # Standard PDB-style formatting for 3Dmol.js/PyMOL compatibility
            f.write(
                f"ATOM  {i:5d}  CA  GLY {chain}{i:4d}    "
                f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{b:6.2f}           C\n"
            )
        f.write("END\n")


def radius_of_gyration(X: np.ndarray) -> float:
    c = X.mean(axis=0)
    return float(np.sqrt(((X - c) ** 2).sum(axis=1).mean()))


def convex_hull_volume(X: np.ndarray) -> float:
    try:
        from scipy.spatial import ConvexHull
        return float(ConvexHull(X).volume)
    except Exception:
        return float("nan")


def make_html(
    out_html: Path,
    hnf4a_pdb_text: str,
    p60_pdb_text: str,
    hnf4a_aligned_pdb_text: str,
    p60_aligned_pdb_text: str,
    displacement_pdb_text: str,
    stats_rows: list,
    chrom: str,
    start: int,
    end: int,
    n_atoms: int,
) -> None:
    region = f"{chrom}:{start/1e6:.2f}–{end/1e6:.2f} Mb"

    def stat_rows_html() -> str:
        s = ""
        for metric, h, p, d, interp in stats_rows:
            htxt = f"{h:.3f}" if not math.isnan(h) else "NA"
            ptxt = f"{p:.3f}" if not math.isnan(p) else "NA"
            dtxt = f"{d:+.3f}" if not math.isnan(d) else "NA"
            if metric == "Convex hull volume":
                htxt += " (norm.)"
            if metric == "Mean displacement":
                s += (
                    f"<tr><td>{metric}</td><td class='blue'>{htxt}</td>"
                    f"<td class='red'>max {ptxt}</td><td class='purple'>mean {htxt}</td><td>{interp}</td></tr>"
                )
            else:
                s += (
                    f"<tr><td>{metric}</td><td class='blue'>{htxt}</td>"
                    f"<td class='red'>{ptxt}</td><td class='purple'>{dtxt}</td><td>{interp}</td></tr>"
                )
        return s

    html = f"""<!doctype html>
<html>
<head>
<meta charset="utf-8">
<title>{chrom} HNF4A vs P60 interactive 3D molecules</title>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<style>
body{{margin:0;font-family:Arial,Helvetica,sans-serif;background:#fff;color:#111}}
.header{{font-size:30px;font-weight:900;padding:15px 18px 5px 18px}}
.meta{{position:absolute;right:22px;top:20px;font-weight:700;font-size:14px}}
.grid{{display:grid;grid-template-columns:1fr 1fr;gap:12px;padding:12px;border-top:1px solid #bbb}}
.panel{{border:1px solid #bbb;border-radius:12px;padding:12px;background:#fbfbfb}}
.badge{{display:inline-flex;align-items:center;justify-content:center;width:42px;height:42px;border-radius:50%;color:white;font-size:27px;font-weight:900;margin-right:10px}}
.bluebg{{background:#0615c9}}.redbg{{background:red}}.purplebg{{background:#6200cc}}.greenbg{{background:green}}
.title{{font-size:22px;font-weight:900;vertical-align:middle}}
.blue{{color:blue;font-weight:900}}.red{{color:red;font-weight:900}}.purple{{color:#6a00d4;font-weight:900}}
.viewer{{height:560px;border:1px solid #ddd;border-radius:10px;background:white;overflow:hidden;margin-top:10px}}
.controls{{text-align:center;margin-top:8px}}
button{{margin:3px;padding:6px 10px;border:1px solid #999;border-radius:6px;background:white;font-weight:700;cursor:pointer}}
button:hover{{background:#eee}}
.colorbarwrap{{display:flex;align-items:center;justify-content:center;gap:16px;font-weight:900;margin-top:10px}}
.colorbar{{width:280px;height:22px;border:1px solid #999;background:linear-gradient(to right,blue,cyan,lime,yellow,orange,red)}}
.alignment{{display:grid;grid-template-columns:1fr 1fr;gap:12px;padding:0 12px 12px 12px}}
.viewer-small{{height:500px;border:1px solid #ddd;border-radius:10px;background:white;overflow:hidden;margin-top:10px}}
table{{border-collapse:collapse;width:100%;font-size:16px;margin-top:15px}}
th,td{{border:1px solid #aaa;padding:12px 8px;text-align:center}}
th{{background:#f5f5f5;font-weight:900}}td:first-child{{font-weight:900}}
.footer{{padding:8px 16px 14px 16px;font-size:14px}}
</style>
</head>
<body>
<div class="header">{region} | HNF4A vs P60 | Interactive 3D GenomeFlow Models</div>
<div class="meta">Atoms/bins: {n_atoms} | Contact B-factor coloring | Alignment: Kabsch</div>

<div class="grid">
  <div class="panel">
    <div><span class="badge bluebg">1</span><span class="title blue">HNF4A 3D molecule</span></div>
    <div id="h_viewer" class="viewer"></div>
    <div class="controls">
      <button onclick="front(vH)">Front</button><button onclick="back(vH)">Back 180°</button><button onclick="spin(vH)">Spin</button>
      <button onclick="contactStyle(vH)">Contact color</button><button onclick="solidStyle(vH,'blue')">Solid blue</button>
    </div>
    <div class="colorbarwrap"><span>Low Contact<br>Frequency</span><span class="colorbar"></span><span>High Contact<br>Frequency</span></div>
  </div>

  <div class="panel">
    <div><span class="badge redbg">2</span><span class="title red">P60 3D molecule</span></div>
    <div id="p_viewer" class="viewer"></div>
    <div class="controls">
      <button onclick="front(vP)">Front</button><button onclick="back(vP)">Back 180°</button><button onclick="spin(vP)">Spin</button>
      <button onclick="contactStyle(vP)">Contact color</button><button onclick="solidStyle(vP,'orange')">Solid orange</button>
    </div>
    <div class="colorbarwrap"><span>Low Contact<br>Frequency</span><span class="colorbar"></span><span>High Contact<br>Frequency</span></div>
  </div>
</div>

<div class="alignment">
  <div class="panel">
    <div><span class="badge purplebg">3</span><span class="title purple">3D alignment: HNF4A aligned onto P60</span></div>
    <div id="align_viewer" class="viewer-small"></div>
    <div class="controls">
      <button onclick="front(vA)">Front</button><button onclick="back(vA)">Back 180°</button><button onclick="spin(vA)">Spin</button>
      <button onclick="alignmentMode()">Blue/orange alignment</button><button onclick="displacementMode()">White/red displacement</button>
    </div>
  </div>

  <div class="panel">
    <div><span class="badge greenbg">4</span><span class="title" style="color:green">Comparison statistics</span></div>
    <table>
      <thead><tr><th>Metric</th><th>HNF4A</th><th>P60 / max</th><th>Δ / mean</th><th>Interpretation</th></tr></thead>
      <tbody>{stat_rows_html()}</tbody>
    </table>
  </div>
</div>

<div class="footer"><b>Notes:</b> P60 uses a standard GLY-chain PDB format and includes a solid-orange fallback if browser gradient coloring fails.</div>

<script>
const pdbH = {json.dumps(hnf4a_pdb_text)};
const pdbP = {json.dumps(p60_pdb_text)};
const pdbHA = {json.dumps(hnf4a_aligned_pdb_text)};
const pdbPA = {json.dumps(p60_aligned_pdb_text)};
const pdbD = {json.dumps(displacement_pdb_text)};

function baseStyle(v, scheme) {{
  v.setStyle({{}}, {{
    sphere: {{radius:0.48, colorscheme:scheme}},
    stick: {{radius:0.18, colorscheme:scheme}}
  }});
  v.render();
}}

function makeSingle(id, pdb) {{
  let v = $3Dmol.createViewer(id, {{backgroundColor:"white"}});
  v.addModel(pdb, "pdb");
  contactStyle(v);
  front(v);
  v.__spinning=false;
  return v;
}}

function contactStyle(v) {{
  baseStyle(v, {{prop:"b", gradient:"roygb", min:0, max:100}});
}}

function solidStyle(v, colorname) {{
  v.setStyle({{}}, {{sphere:{{radius:0.48,color:colorname}}, stick:{{radius:0.18,color:colorname}}}});
  v.render();
}}

function makeAlignment(id) {{
  let v = $3Dmol.createViewer(id, {{backgroundColor:"white"}});
  v.addModel(pdbHA, "pdb");
  v.addModel(pdbPA, "pdb");
  alignmentStyle(v);
  front(v);
  v.__spinning=false;
  return v;
}}

function alignmentStyle(v) {{
  v.setStyle({{model:0}}, {{sphere:{{radius:0.33,color:"blue"}}, stick:{{radius:0.12,color:"blue"}}}});
  v.setStyle({{model:1}}, {{sphere:{{radius:0.33,color:"orange"}}, stick:{{radius:0.12,color:"orange"}}}});
  v.render();
}}

function front(v) {{
  v.stopAnimate();
  v.__spinning=false;
  v.zoomTo();
  v.rotate(18, {{x:1,y:0,z:0}});
  v.rotate(-58, {{x:0,y:1,z:0}});
  v.render();
}}

function back(v) {{
  front(v);
  v.rotate(180, {{x:0,y:1,z:0}});
  v.render();
}}

function spin(v) {{
  if(v.__spinning) {{ v.stopAnimate(); v.__spinning=false; }}
  else {{ v.spin("y",1); v.__spinning=true; }}
}}

let vH = makeSingle("h_viewer", pdbH);
let vP = makeSingle("p_viewer", pdbP);
let vA = makeAlignment("align_viewer");

function alignmentMode() {{
  vA.removeAllModels();
  vA.addModel(pdbHA, "pdb");
  vA.addModel(pdbPA, "pdb");
  alignmentStyle(vA);
  front(vA);
}}

function displacementMode() {{
  vA.removeAllModels();
  vA.addModel(pdbD, "pdb");
  vA.setStyle({{}}, {{
    sphere:{{radius:0.48,colorscheme:{{prop:"b",gradient:"rwb",min:0,max:100}}}},
    stick:{{radius:0.18,colorscheme:{{prop:"b",gradient:"rwb",min:0,max:100}}}}
  }});
  front(vA);
}}
</script>
</body>
</html>
"""
    out_html.write_text(html)


def build(args):
    outdir = Path(args.out_prefix + "_files")
    outdir.mkdir(exist_ok=True)

    Hraw = read_coords(Path(args.hnf4a))
    Praw = read_coords(Path(args.p60))
    if len(Hraw) == 0 or len(Praw) == 0:
        raise RuntimeError("Could not parse one or both PDB files.")

    n = min(len(Hraw), len(Praw))
    H = normalize(Hraw[:n])
    P = normalize(Praw[:n])
    H_al = kabsch(H, P)
    disp = np.linalg.norm(H_al - P, axis=1)
    disp_b = 100 * disp / (disp.max() if disp.max() > 0 else 1.0)

    HB = contact_bfactor(H)
    PB = contact_bfactor(P)

    h_pdb = outdir / "HNF4A_contact.pdb"
    p_pdb = outdir / "P60_contact.pdb"
    ha_pdb = outdir / "HNF4A_aligned.pdb"
    pa_pdb = outdir / "P60_aligned.pdb"
    pd_pdb = outdir / "P60_displacement.pdb"

    write_pdb(h_pdb, H, HB, "A")
    write_pdb(p_pdb, P, PB, "A")
    write_pdb(ha_pdb, H_al, np.zeros(n), "A")
    write_pdb(pa_pdb, P, np.zeros(n), "B")
    write_pdb(pd_pdb, P, disp_b, "B")

    DH, DP = pairwise(H), pairwise(P)
    vh, vp = DH[np.triu_indices(n, 1)], DP[np.triu_indices(n, 1)]
    Hvol, Pvol = convex_hull_volume(H), convex_hull_volume(P)
    H_e2e = float(np.linalg.norm(H[0] - H[-1]))
    P_e2e = float(np.linalg.norm(P[0] - P[-1]))

    stats_rows = [
        ["Radius of gyration (Rg)", radius_of_gyration(H), radius_of_gyration(P), radius_of_gyration(P)-radius_of_gyration(H), "P60 more compact" if radius_of_gyration(P) < radius_of_gyration(H) else "HNF4A more compact"],
        ["End-to-end distance", H_e2e, P_e2e, P_e2e-H_e2e, "P60 more compact" if P_e2e < H_e2e else "HNF4A more compact"],
        ["Convex hull volume", 1.0, Pvol/Hvol if Hvol and not math.isnan(Hvol) else float("nan"), (Pvol/Hvol-1.0) if Hvol and not math.isnan(Hvol) else float("nan"), "P60 more compact" if Hvol and not math.isnan(Hvol) and Pvol/Hvol < 1 else "HNF4A more compact"],
        ["Contact fraction", float(np.mean(vh < 1.0)), float(np.mean(vp < 1.0)), float(np.mean(vp < 1.0)-np.mean(vh < 1.0)), "P60 more contacts" if np.mean(vp < 1.0) > np.mean(vh < 1.0) else "HNF4A more contacts"],
        ["Mean displacement", float(np.mean(disp)), float(np.max(disp)), float(np.mean(disp)), "3D alignment shift"],
    ]

    stats_csv = outdir / "HNF4A_vs_P60_statistics.csv"
    pd.DataFrame(stats_rows, columns=["Metric","HNF4A","P60_or_max","Delta_or_mean","Interpretation"]).to_csv(stats_csv, index=False)

    html_path = Path(args.out_prefix + ".html")
    make_html(
        html_path,
        h_pdb.read_text(),
        p_pdb.read_text(),
        ha_pdb.read_text(),
        pa_pdb.read_text(),
        pd_pdb.read_text(),
        stats_rows,
        args.chrom,
        args.start,
        args.end,
        n,
    )

    zip_path = Path(args.out_prefix + ".zip")
    with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as z:
        for f in [html_path, h_pdb, p_pdb, ha_pdb, pa_pdb, pd_pdb, stats_csv]:
            z.write(f, arcname=f.name)

    print("Saved:")
    print(html_path)
    print(zip_path)
    print(stats_csv)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hnf4a", required=True)
    parser.add_argument("--p60", required=True)
    parser.add_argument("--chrom", default="chr8")
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--end", type=int, required=True)
    parser.add_argument("--out-prefix", default="HNF4A_P60_interactive")
    args = parser.parse_args()
    build(args)


if __name__ == "__main__":
    main()
