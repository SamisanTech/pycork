"""Drive run_nm_batch.py per file; isolate native crashes."""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

DIR = Path(r"E:\github.com\SamisanTech\manifold\test\Repair Test files\Non Manifold")
OUT_DIR = DIR / "pycork_repair"
PY = r"C:\Users\Karan\AppData\Local\Programs\Python\Python39\python.exe"
CHILD = Path(r"E:\github.com\pycork\run_nm_batch.py")
TIMEOUT = 600
OUT_JSON = Path(r"E:\github.com\pycork\nm_batch_results.json")
SKIP = {"20.stl"}


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    files = sorted(
        [
            p
            for p in DIR.glob("*.stl")
            if not p.name.endswith("_repair.stl") and p.name not in SKIP
        ],
        key=lambda p: p.stat().st_size,
    )
    print(f"{len(files)} files (skip {sorted(SKIP)}), out {OUT_DIR}, timeout {TIMEOUT}s", flush=True)
    results = []
    for i, path in enumerate(files, 1):
        mb = path.stat().st_size / 1e6
        print(f"\n=== [{i}/{len(files)}] {path.name}  {mb:.1f} MB ===", flush=True)
        rec = {
            "file": path.name,
            "in_bytes": path.stat().st_size,
            "ok": False,
        }
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
        rec["stderr_tail"] = (proc.stderr or "")[-800:]
        line = ""
        for ln in (proc.stdout or "").splitlines()[::-1]:
            if ln.startswith("{"):
                line = ln
                break
        if line:
            try:
                rec = json.loads(line)
                rec.setdefault("returncode", proc.returncode)
            except json.JSONDecodeError:
                rec["error"] = "bad json"
                rec["stdout_tail"] = (proc.stdout or "")[-400:]
        elif proc.returncode != 0:
            rec["error"] = f"crash rc={proc.returncode:#x}" if proc.returncode < 0 or proc.returncode > 255 else f"crash rc={proc.returncode}"
            if proc.returncode == 3221225477:
                rec["error"] = "native crash ACCESS_VIOLATION"
        else:
            rec["error"] = "no json output"
            rec["stdout_tail"] = (proc.stdout or "")[-400:]
        status = "OK" if rec.get("ok") else f"FAIL {rec.get('error')}"
        extra = ""
        if rec.get("ok"):
            extra = (
                f"  in V={rec['in_V']:,} F={rec['in_F']:,} open={rec['in_open']} nm={rec['in_nm']}"
                f"  -> V={rec['out_V']:,} F={rec['out_F']:,} open={rec['out_open']} nm={rec['out_nm']}"
                f" shells={rec['shells']} wt={rec['watertight']} {rec['repair_s']}s"
            )
        print(f"{status}{extra}", flush=True)
        results.append(rec)
        OUT_JSON.write_text(json.dumps(results, indent=2), encoding="utf-8")

    print(f"\nWrote {OUT_JSON}", flush=True)


if __name__ == "__main__":
    main()
