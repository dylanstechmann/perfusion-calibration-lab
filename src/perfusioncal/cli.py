"""Analyze a CSV or generate synthetic calibration traces."""

import argparse
import csv
import json
from pathlib import Path
import platform
import shutil
from tempfile import TemporaryDirectory

import numpy as np

from perfusioncal import __version__
from perfusioncal.analysis import analyze, markdown


def write_demo(path, seed=0):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    with path.open("x", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(["run_id", "time_s", "mass_mg", "target_flow_ul_min"])
        for target in [25, 50, 100]:
            for repeat in range(5):
                flow = target * (0.93 + rng.normal(0, 0.01))
                for t in range(0, 301, 10):
                    # Synthetic fluid density is exactly 1 mg/µL in this fixture.
                    writer.writerow([f"synthetic-{target}-{repeat}", t,
                                     7.0 + flow * t / 60 + rng.normal(0, 0.1), target])


def main(argv=None):
    parser = argparse.ArgumentParser(description="Offline gravimetric calibration for research syringe pumps")
    sub = parser.add_subparsers(dest="command", required=True)
    demo = sub.add_parser("demo")
    demo.add_argument("--out", required=True)
    demo.add_argument("--seed", type=int, default=0)
    run = sub.add_parser("analyze")
    run.add_argument("csv")
    run.add_argument("--density-mg-ul", type=float, required=True, help="density for your fluid and temperature")
    run.add_argument("--discard-seconds", type=float, default=0)
    run.add_argument("--seed", type=int, default=0)
    run.add_argument("--out", required=True, help="new output directory")
    args = parser.parse_args(argv)
    try:
        if args.command == "demo":
            write_demo(args.out, args.seed)
        else:
            report = analyze(args.csv, density_mg_ul=args.density_mg_ul,
                             discard_seconds=args.discard_seconds, seed=args.seed)
            report["environment"] = {"python": platform.python_version(), "numpy": np.__version__,
                                      "perfusioncal": __version__}
            output = Path(args.out)
            json_text = json.dumps(report, indent=2, allow_nan=False) + "\n"
            markdown_text = markdown(report)
            if output.exists() or output.is_symlink():
                raise FileExistsError(f"output already exists: {output}")
            output.parent.mkdir(parents=True, exist_ok=True)
            with TemporaryDirectory(prefix=f".{output.name}-", dir=output.parent) as temporary:
                staged = Path(temporary)
                (staged / "report.json").write_text(json_text, encoding="utf-8")
                (staged / "REPORT.md").write_text(markdown_text, encoding="utf-8")
                output.mkdir()
                try:
                    for artifact in (staged / "report.json", staged / "REPORT.md"):
                        artifact.rename(output / artifact.name)
                except BaseException:
                    shutil.rmtree(output)
                    raise
    except (ValueError, OSError) as exc:
        parser.error(str(exc))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
