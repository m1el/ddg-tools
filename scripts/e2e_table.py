#!/usr/bin/env python3
"""Markdown tables from the .jsonl files of scripts/e2e-developable.sh.

    python3 scripts/e2e_table.py runs_tmp/e2e-developable [--all-settings]
"""
import glob
import json
import os
import sys


def load(directory):
    runs = {}
    for path in sorted(glob.glob(os.path.join(directory, "*.jsonl"))):
        recs = [json.loads(line) for line in open(path)]
        budget = next((r["iterations"] for r in recs if r["stage"] != "flow"), None)
        for r in recs:
            if r["stage"] == "flow" and "iterations_done" not in r:  # older files: the field was overwritten
                r["iterations_done"] = r["iterations"]
                r["iterations"] = budget
            runs.setdefault((r["case"], r["subdivisions"], r["iterations"]), []).append(r)
    return runs


def stop_text(flow):
    s = flow.get("stopped")
    if s:
        return "stopped before round {}: {}".format(s["round"] + 1, s["reason"])
    statuses = [x["status"] for x in flow["rounds"]]
    short = {"gradient_tolerance": "grad", "energy_stalled": "stall", "max_iterations": "max",
             "line_search_failed": "LS-FAIL", "stationary": "stat", "cancelled": "cancel", "non_descent": "NON-DESC"}
    return "/".join(short.get(x, x) for x in statuses)


def main():
    directory = sys.argv[1]
    every = "--all-settings" in sys.argv
    runs = load(directory)
    print("| mesh | sub | tris | flow iters | rounds end | E before | E after | bad | flow s | "
          "setting | seams | terminals | charts | cones | overlap | splits | mirrored | UV s |")
    print("|" + "---|" * 18)
    for (case, sub, its), records in sorted(runs.items(), key=lambda kv: (kv[0][0], -kv[0][1], kv[0][2])):
        flow = next((r for r in records if r["stage"] == "flow"), None)
        for r in records:
            if r["stage"] != "seams":
                continue
            if not every and r["setting"] != "default":
                continue
            if "error" in r:
                print("| {} | {} | | | | | | | | {} | ERROR: {} |".format(
                    case, sub, r["setting"], r["error"][:80]))
                continue
            if flow and "error" not in flow:
                f = "{} | {} | {:.3g} | {:.3g} | {} | {} | {:.1f}".format(
                    flow["triangles_out"], flow["iterations_done"], flow["energy_initial"],
                    flow["energy_final"], stop_text(flow), flow["bad_faces"], flow["seconds"])
                f = "{} | {} | {} | {} | {} | {} | {} | {}".format(
                    case, sub, flow["triangles_out"], flow["iterations_done"], stop_text(flow),
                    "{:.3g}".format(flow["energy_initial"]), "{:.3g}".format(flow["energy_final"]),
                    flow["bad_faces"]) + " | {:.1f}".format(flow["seconds"])
            elif flow:
                f = "{} | {} | ERROR: {} | | | | | |".format(case, sub, flow["error"][:60])
            else:
                f = "{} | -1 | no flow | | | | | | ".format(case)
            ov = r["overlap"]
            print("| {} | {} | {} | {} | {} | {} | {} | {} | {} | {:.2f} |".format(
                f, r["setting"], r["seam_edges"], r["terminals"], r["charts"], r["cones"],
                "n/a" if ov is None else "{:.1f} %".format(100 * ov),
                "-" if r.get("splits") is None else r["splits"], r["mirrored"], r["seconds"]))


main()
