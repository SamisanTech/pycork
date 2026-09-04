"""Repair every Non Manifold STL. Isolated child per file. No skips."""
from __future__ import annotations

import json
import subprocess
from pathlib import Path

PY = r"C:\Users\Karan\AppData\Local\Programs\Python\Python39\python.exe"
CHILD = Path(r"E:\github.com\pycork\run_nm_batch.py")
FOLDER = Path(r"E:\github.com\SamisanTech\manifold\test\Repair Test files\Non Manifold")
OUT_JSON = Path(r"E:\github.com\pycork\nm_all_results.json")
TIMEOUT = 300


def main():
    files = sorted(
        (p for p in FOLDER.glob("*.stl") if not p.name.endswith("_repair.stl")),
        key=lambda p: (p.stat().st_size, p.name),
    )
    out_dir = FOLDER / "pycork_repair"
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"{len(files)} Non Manifold files, timeout {TIMEOUT}s, out {out_dir}", flush=True)
    results = []
    for i, path in enumerate(files, 1):
        mb = path.stat().st_size / 1e6
        print(f"\n=== [{i}/{len(files)}] {path.name}  {mb:.1f} MB ===", flush=True)
        rec = {
            "file": path.name,
            "folder": "Non Manifold",
            "ok": False,
            "in_bytes": path.stat().st_size,
        }
        try:
            proc = subprocess.run(
                [PY, "-u", str(CHILD), str(path), str(out_dir)],
                capture_output=True,
                text=True,
                timeout=TIMEOUT,
            )
        except subprocess.TimeoutExpired as e:
            rec["error"] = f"timeout {TIMEOUT}s"
            rec["ok"] = False
            print(f"TIMEOUT {path.name}", flush=True)
            results.append(rec)
            OUT_JSON.write_text(json.dumps(results, indent=2), encoding="utf-8")
            if getattr(e, "process", None) is not None:
                try:
                    e.process.kill()
                except OSError:
                    pass
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
                rec["folder"] = "Non Manifold"
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
                f"shells={rec['shells']} wt={rec['watertight']} vol={rec.get('volume')}",
                flush=True,
            )
        else:
            print(f"FAIL {rec.get('error')}", flush=True)
        results.append(rec)
        OUT_JSON.write_text(json.dumps(results, indent=2), encoding="utf-8")
    print(f"\nWrote {OUT_JSON}", flush=True)


if __name__ == "__main__":
    main()
