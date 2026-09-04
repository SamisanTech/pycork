"""Isolated repair batch. One _repair.stl per input, in <folder>/pycork_repair."""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

PY = r"C:\Users\Karan\AppData\Local\Programs\Python\Python39\python.exe"
CHILD = Path(r"E:\github.com\pycork\run_nm_batch.py")
OUT_JSON = Path(r"E:\github.com\pycork\working_results.json")
TIMEOUT = 180
# Hang-only skips. 27 / 28_danger finished last NM-all pass.
SKIP = {
    "Non Manifold": {"20.stl", "26.stl"},
    "Stitching": set(),
    "Impossible": {
        "MFUE2460(2) J2.stl",
        "G40601=P.stl",
        "29.stl",  # 11 MB, resolve hang >180s
        "19.stl",  # 22 MB, resolve hang >180s
        "1.stl",   # 27 MB, resolve hang >180s
    },
}
FOLDERS = [
    Path(r"E:\github.com\SamisanTech\manifold\test\Repair Test files\Non Manifold"),
    Path(r"E:\github.com\SamisanTech\manifold\test\Repair Test files\Stitching"),
    Path(r"E:\github.com\SamisanTech\manifold\test\Repair Test files\Impossible"),
]


def main():
    files = []
    for folder in FOLDERS:
        deny = SKIP.get(folder.name, set())
        for p in folder.glob("*.stl"):
            if p.name.endswith("_repair.stl") or p.name in deny:
                continue
            files.append(p)
    files.sort(key=lambda p: (p.stat().st_size, p.name))
    skipped = {k: sorted(v) for k, v in SKIP.items() if v}
    print(f"{len(files)} files, skip {skipped}, timeout {TIMEOUT}s", flush=True)
    results = []
    for i, path in enumerate(files, 1):
        out_dir = path.parent / "pycork_repair"
        out_dir.mkdir(parents=True, exist_ok=True)
        mb = path.stat().st_size / 1e6
        print(f"\n=== [{i}/{len(files)}] {path.parent.name}/{path.name}  {mb:.1f} MB ===", flush=True)
        rec = {"file": path.name, "folder": path.parent.name, "ok": False, "in_bytes": path.stat().st_size}
        try:
            proc = subprocess.run(
                [PY, "-u", str(CHILD), str(path), str(out_dir)],
                capture_output=True,
                text=True,
                timeout=TIMEOUT,
            )
        except subprocess.TimeoutExpired:
            rec["error"] = f"timeout {TIMEOUT}s"
            print(f"TIMEOUT {path.name}", flush=True)
            results.append(rec)
            OUT_JSON.write_text(json.dumps(results, indent=2), encoding="utf-8")
            continue
        rec["returncode"] = proc.returncode
        rec["stderr_tail"] = (proc.stderr or "")[-400:]
        line = ""
        for ln in (proc.stdout or "").splitlines()[::-1]:
            if ln.startswith("{"):
                line = ln
                break
        if line:
            try:
                rec = json.loads(line)
                rec["folder"] = path.parent.name
            except json.JSONDecodeError:
                rec["error"] = "bad json"
        elif proc.returncode:
            rec["error"] = f"crash rc={proc.returncode}"
        else:
            rec["error"] = "no json"
        if rec.get("ok"):
            print(
                f"OK  {rec['repair_s']}s  F {rec['in_F']:,}->{rec['out_F']:,}  "
                f"nm {rec['in_nm']}->{rec['out_nm']}  open {rec['in_open']}->{rec['out_open']}  "
                f"shells={rec['shells']} wt={rec['watertight']}",
                flush=True,
            )
        else:
            print(f"FAIL {rec.get('error')}", flush=True)
        results.append(rec)
        OUT_JSON.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(f"\nWrote {OUT_JSON}", flush=True)


if __name__ == "__main__":
    main()
