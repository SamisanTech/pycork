"""Side-by-side report of compare_out/orig.json vs compare_out/fast.json."""
import json, sys
from pathlib import Path

d = Path(r"E:\github.com\pycork\compare_out")
A = json.loads((d / "orig.json").read_text())
B = json.loads((d / "fast.json").read_text())

scalar = ["seconds", "V", "F", "unique_positions", "volume", "area", "shells", "unique_edges",
          "open_edges", "nonmanifold_edges", "duplicate_faces", "degenerate_faces", "euler_number",
          "watertight", "winding_consistent", "pycork_isSolid"]
vec = ["bbox_min", "bbox_max", "centroid", "center_mass"]

def fmt(v):
    if isinstance(v, bool): return str(v)
    if isinstance(v, int): return f"{v:,}"
    if isinstance(v, float): return f"{v:.6f}"
    return str(v)

def rng(runs, k):
    vals = [r[k] for r in runs]
    if all(isinstance(v, bool) for v in vals):
        return "/".join(str(v) for v in vals)
    lo, hi = min(vals), max(vals)
    mean = sum(vals) / len(vals)
    if isinstance(lo, int):
        return f"{lo:,} .. {hi:,}  (mean {mean:,.0f})"
    return f"{lo:.6f} .. {hi:.6f}  (mean {mean:.6f})"

lines = []
lines.append(f"INPUT 21.stl: V={A['input']['V']:,} F={A['input']['F']:,} open={A['input']['open_edges']} "
             f"nm={A['input']['nonmanifold_edges']} volume={A['input']['volume']:.6f} area={A['input']['area']:.6f} "
             f"shells={A['input']['shells']} watertight={A['input']['watertight']}")
lines.append("")
lines.append(f"{'metric':20s} | {'ORIGINAL pycork (3 runs)':44s} | {'OPTIMIZED pycork (3 runs)':44s} | rel. diff of means")
lines.append("-" * 140)
for k in scalar:
    a, b = rng(A["runs"], k), rng(B["runs"], k)
    va = [r[k] for r in A["runs"]]; vb = [r[k] for r in B["runs"]]
    if all(isinstance(v, bool) for v in va):
        diff = "same" if set(va) == set(vb) else "DIFF"
    else:
        ma, mb = sum(va)/len(va), sum(vb)/len(vb)
        diff = f"{(mb-ma)/ma*100:+.4f}%" if ma else f"{mb-ma:+g}"
    lines.append(f"{k:20s} | {a:44s} | {b:44s} | {diff}")
for k in vec:
    a = A["runs"][0][k]; b = B["runs"][0][k]
    lines.append(f"{k:20s} | {str([round(x,6) for x in a]):44s} | {str([round(x,6) for x in b]):44s} | "
                 f"max |d|={max(abs(x-y) for x,y in zip(a,b)):.2e}")
ha = A["runs"][0]["edge_hist"]; hb = B["runs"][0]["edge_hist"]
lines.append(f"{'edge_hist run0':20s} | {str(ha):44s} | {str(hb):44s} |")
lines.append(f"{'shell sizes top5':20s} | {str(A['runs'][0]['shell_sizes_top10'][:5]):44s} | {str(B['runs'][0]['shell_sizes_top10'][:5]):44s} |")
# run-to-run spread of the original itself, for reference
va = [r["V"] for r in A["runs"]]
lines.append("")
lines.append(f"Reference: original's own run-to-run spread in V is {max(va)-min(va):,} vertices "
             f"({(max(va)-min(va))/min(va)*100:.3f}%) because cork randomly perturbs positions each call.")
txt = "\n".join(lines)
print(txt)
(d / "REPORT.txt").write_text(txt)
