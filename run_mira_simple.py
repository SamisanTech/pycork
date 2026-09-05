"""Default pycork.repair on Mira Simple repair/*.stl → pycork_repair/."""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

PY = r"C:\Users\Karan\AppData\Local\Programs\Python\Python39\python.exe"
CHILD = Path(r"E:\github.com\pycork\run_nm_batch.py")
FOLDER = Path(
    r"E:\github.com\SamisanTech\manifold\test\Repair Test files"
    r"\Non Manifold\Mira Simple repair"
)
OUT_DIR = FOLDER / "pycork_repair"
OUT_JSON = Path(r"E:\github.com\pycork\mira_simple_results.json")
OUT_TXT = Path(r"E:\github.com\pycork\mira_simple_results.txt")
TIMEOUT = 180


def main():
    files = sorted(
        [p for p in FOLDER.glob("*.stl") if not p.name.endswith("_repair.stl")],
        key=lambda p: (p.stat().st_size, p.name),
    )
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    print(f"{len(files)} Mira Simple files, timeout {TIMEOUT}s, out {OUT_DIR}", flush=True)
    results = []
    for i, path in enumerate(files, 1):
        mb = path.stat().st_size / 1e6
        print(f"\n=== [{i}/{len(files)}] {path.name}  {mb:.1f} MB ===", flush=True)
        rec = {"file": path.name, "ok": False, "in_bytes": path.stat().st_size}
        try:
            proc = subprocess.run(
                [PY, "-u", str(CHILD), str(path), str(OUT_DIR)],
                capture_output=True,
                text=True,
                timeout=TIMEOUT,
            )
        except subprocess.TimeoutExpired:
            rec["error"] = f"timeout {TIMEOUT}s"
            print(f"TIMEOUT {path.name}", flush=True)
            results.append(rec)
            continue
        rec["returncode"] = proc.returncode
        rec["stderr_tail"] = (proc.stderr or "")[-300:]
        line = ""
        for ln in (proc.stdout or "").splitlines()[::-1]:
            if ln.startswith("{"):
                line = ln
                break
        if line:
            try:
                rec = json.loads(line)
            except json.JSONDecodeError:
                rec["error"] = "bad json"
        elif proc.returncode:
            rec["error"] = f"crash rc={proc.returncode}"
        else:
            rec["error"] = "no json"
        if rec.get("ok"):
            print(
                f"  {rec['repair_s']:.3f}s  "
                f"F {rec['in_F']:,}->{rec['out_F']:,}  "
                f"open {rec['in_open']}->{rec['out_open']}  "
                f"nm {rec['in_nm']}->{rec['out_nm']}  "
                f"shells={rec['shells']}  wt={rec['watertight']}  "
                f"vol={rec.get('volume')}",
                flush=True,
            )
        else:
            print(f"  FAIL {rec.get('error')}", flush=True)
        results.append(rec)
        OUT_JSON.write_text(json.dumps(results, indent=2), encoding="utf-8")

    lines = [
        "Mira Simple repair — default pycork.repair",
        f"STLs: {OUT_DIR}\\{{stem}}_repair.stl",
        "",
        f"{'file':<22} {'s':>7} {'F in':>10} {'F out':>10} {'open':>12} {'nm':>12} "
        f"{'sh':>3} {'wt':>3} {'vol':>12}",
    ]
    for r in results:
        if not r.get("ok"):
            lines.append(f"{r['file']:<22}  {r.get('error', 'FAIL')}")
            continue
        lines.append(
            f"{r['file']:<22} {r['repair_s']:7.3f} {r['in_F']:10,} {r['out_F']:10,} "
            f"{r['in_open']:5}->{r['out_open']:<5} {r['in_nm']:5}->{r['out_nm']:<5} "
            f"{r['shells']:3} {'YES' if r['watertight'] else 'no ':3} "
            f"{(r.get('volume') if r.get('volume') is not None else 0):12.4f}"
        )
    text = "\n".join(lines) + "\n"
    OUT_TXT.write_text(text, encoding="utf-8")
    print("\n" + text, flush=True)
    print("wrote", OUT_TXT, flush=True)


if __name__ == "__main__":
    main()
