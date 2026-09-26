from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/prism-matplotlib")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path("analysis/california")
OUTPUT = ROOT / "output"
TABLES = OUTPUT / "tables"
DIAGNOSTICS = OUTPUT / "diagnostics"
FIGURES = OUTPUT / "figures"
SOLVERS = OUTPUT / "solver-logs"
RESULTS_ARCHIVE = Path("prism-uploads/california_trsc_results.zip")
RESULTS_ARCHIVE_SHA256 = "c86dc8fba654cc96571cb5f0656c382390439ed98220c705925ac228c5999923"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


fibers = pd.read_csv(TABLES / "evaluation_fibers.csv")
modulus = pd.read_csv(TABLES / "deployment_modulus.csv")
primary_modulus_zero = modulus[
    (modulus["geometry_role"] == "primary_calibration")
    & (modulus["epsilon"] == 0)
]

summary_rows = []
main_rows = []
for q in (1, 2):
    selection = pd.read_csv(TABLES / f"q{q}_selection.csv").iloc[0]
    radius = pd.read_csv(TABLES / f"q{q}_rmin_scale.csv").iloc[0]
    dynamic = pd.read_csv(TABLES / f"q{q}_dynamic_bounds.csv")
    base = dynamic[dynamic["c"] == dynamic["c"].min()]
    final = dynamic[dynamic["c"] == dynamic["c"].max()]
    fiber = fibers[
        (fibers["q"] == q)
        & (fibers["geometry_role"] == "primary_calibration")
    ]
    solver = pd.read_csv(SOLVERS / f"q{q}_dynamic_support_solver.csv")
    counts = solver["status"].value_counts()
    methods = solver["certificate_method"].value_counts()
    solver_status = "; ".join(f"{name}={count}" for name, count in counts.items())
    solver_status += "; certificates: " + ", ".join(
        f"{name}={count}" for name, count in methods.items()
    )

    summary_rows.append(
        {
            "q": q,
            "selected_kappa": selection["kappa"],
            "selected_h": int(selection["h"]),
            "validation_score": selection["validation_score"],
            "R_min_lower": radius["lower"],
            "R_min_upper": radius["upper"],
            "median_dynamic_width": base["width"].median(),
            "max_dynamic_width": base["width"].max(),
            "median_evaluation_fiber_width": fiber["width"].median(),
            "median_deployment_modulus_eps0": primary_modulus_zero["omega"].median(),
            "inference_tier": "Tier 1 plug-in identification only",
            "solver_status": solver_status,
            "dynamic_width_basis": "c=0; R=certified R_min upper",
            "median_dynamic_hull_fraction": base["hull_fraction"].median(),
            "max_dynamic_hull_fraction": base["hull_fraction"].max(),
            "largest_grid_c": final["c"].iloc[0],
            "largest_grid_R": final["R"].iloc[0],
            "largest_grid_median_width": final["width"].median(),
            "largest_grid_max_width": final["width"].max(),
            "largest_grid_median_hull_fraction": final["hull_fraction"].median(),
            "largest_grid_max_hull_fraction": final["hull_fraction"].max(),
        }
    )

    for c_value in (0, 1, 16, 1024):
        subset = dynamic[dynamic["c"] == c_value]
        if subset.empty:
            continue
        main_rows.append(
            {
                "q": q,
                "c": c_value,
                "R": subset["R"].iloc[0],
                "median_width": subset["width"].median(),
                "max_width": subset["width"].max(),
                "median_hull_fraction": subset["hull_fraction"].median(),
                "horizon_1_width": subset.loc[subset["horizon"] == 1, "width"].iloc[0],
                "horizon_3_width": subset.loc[subset["horizon"] == 3, "width"].iloc[0],
                "horizon_6_width": subset.loc[subset["horizon"] == 6, "width"].iloc[0],
                "horizon_12_width": subset.loc[subset["horizon"] == 12, "width"].iloc[0],
            }
        )

pd.DataFrame(summary_rows).to_csv(OUTPUT / "california_summary.csv", index=False)
pd.DataFrame(main_rows).to_csv(TABLES / "table_california_main_audited.csv", index=False)

solver_summary_rows = []
solver_flag_rows = []
solver_files = {
    "q1_dynamic": SOLVERS / "q1_dynamic_support_solver.csv",
    "q2_dynamic": SOLVERS / "q2_dynamic_support_solver.csv",
    "evaluation_fibers": SOLVERS / "evaluation_fibers_solver.csv",
    "deployment_modulus": SOLVERS / "deployment_modulus_solver.csv",
}
for component, path in solver_files.items():
    values = pd.read_csv(path)
    status = values["status"].value_counts()
    methods = values["certificate_method"].value_counts()
    solver_summary_rows.append(
        {
            "component": component,
            "endpoints": len(values),
            "optimal": status.get("Optimal solution found", 0),
            "close_to_optimal": status.get("Close to optimal solution found", 0),
            "numerical_problems": status.get("Ran into numerical problems", 0),
            "multiplier_cone": status.get("PROBLEM: Multipliers leaving the cone", 0),
            "max_iterations": status.get("Maximum number of iterations reached", 0),
            "structural_infeasible": status.get("structural_infeasible", 0),
            "rational_dual_residual": methods.get("rational_dual_residual", 0),
            "simplex_hull": methods.get("simplex_hull", 0),
            "max_primal_violation": values["primal_violation"].max(),
            "max_stationarity_residual": values["stationarity_residual"].max(),
        }
    )
    flagged = values[
        values["status"].isin(
            [
                "Ran into numerical problems",
                "PROBLEM: Multipliers leaving the cone",
                "Maximum number of iterations reached",
                "structural_infeasible",
            ]
        )
        | (values["certificate_method"] == "simplex_hull")
    ].copy()
    if len(flagged):
        flagged.insert(0, "component", component)
        solver_flag_rows.append(flagged)

pd.DataFrame(solver_summary_rows).to_csv(
    SOLVERS / "solver_audit_summary.csv", index=False
)
if solver_flag_rows:
    pd.concat(solver_flag_rows, ignore_index=True).to_csv(
        SOLVERS / "solver_flagged_endpoints.csv", index=False
    )

reconciliation = pd.DataFrame(
    [
        ("run_manifest.json", "stale", "Records the earlier Python-only run; completed R stages postdate it."),
        ("solver-logs/all_endpoint_status.csv", "stale", "Records TRSC as unexecuted despite completed endpoint logs."),
        ("diagnostics/fairness_audit.csv", "stale", "Pre-R audit language says TRSC stages were unexecuted."),
        ("diagnostics/figure_status.csv", "stale", "Predates the generated TRSC figures."),
        ("tables/table2_pre_treatment_selection.csv", "stale_placeholder", "Use q1_selection.csv and q2_selection.csv."),
        ("tables/table4_q_R_identification.csv", "stale_placeholder", "Use q1_dynamic_bounds.csv and q2_dynamic_bounds.csv."),
        ("tables/table5_evaluation_deployment.csv", "stale_placeholder", "Use evaluation_fibers.csv and deployment_modulus.csv."),
        ("tables/table6_inference_decomposition.csv", "stale_placeholder", "No Tier-2 or Tier-3 inference was run."),
        ("executed package tree", "checksum_unavailable", "The macOS run manually patched the labelled-state reader and q=2 boundary helper after the stale manifest was written."),
    ],
    columns=["file", "audit_status", "resolution"],
)
reconciliation.to_csv(OUTPUT / "audit_reconciliation.csv", index=False)

legacy_manifest_path = OUTPUT / "run_manifest.json"
audited_manifest_path = OUTPUT / "run_manifest_completed_audit.json"
if legacy_manifest_path.exists():
    old_manifest = json.loads(legacy_manifest_path.read_text())
elif audited_manifest_path.exists():
    prior_audit = json.loads(audited_manifest_path.read_text())
    old_manifest = {
        "package_base_git_sha": prior_audit.get("package_base_git_sha"),
        "package_working_tree_sha256": prior_audit.get(
            "pre_run_manifest_package_working_tree_sha256"
        ),
        "manuscript_sha256": prior_audit.get("manuscript_sha256"),
        "dataset_sha256": prior_audit.get("dataset_sha256"),
        "config_sha256": prior_audit.get("config_sha256"),
    }
else:
    raise FileNotFoundError(
        "Neither the archived run manifest nor the audited manifest is available."
    )
execution_provenance = pd.DataFrame(
    [
        ("execution_host", "macOS arm64", "From the supplied terminal history"),
        ("R_runtime", "R 4.5 arm64 framework; patch version not recorded", "Do not infer an exact patch release"),
        ("TRSC_version", "0.1.0", "Reported by R CMD INSTALL"),
        ("ECOSolveR", "0.6.1", "Binary package downloaded before execution"),
        ("gmp", "0.7-5.1", "Binary package downloaded before execution"),
        ("data_reader_patch", "haven::as_factor(raw$state)", "Required to decode labelled Stata state names"),
        ("boundary_helper_patch", "manually added and exported", "Saved h=1 output agrees with adjacent 1985/1984 fitted endpoints"),
        ("executed_package_tree_sha256", "unavailable", "The package source directory was not archived after the manual patches"),
    ],
    columns=["item", "value", "evidence"],
)
execution_provenance.to_csv(OUTPUT / "execution_provenance.csv", index=False)

completed_manifest = {
    "audit_timestamp_utc": datetime.now(timezone.utc).isoformat(),
    "source_results_archive": "prism-uploads/california_trsc_results.zip",
    "source_results_archive_sha256": (
        sha256(RESULTS_ARCHIVE) if RESULTS_ARCHIVE.exists()
        else RESULTS_ARCHIVE_SHA256
    ),
    "original_manifest_preserved": True,
    "original_manifest_status": "stale_python_only_state",
    "package_base_git_sha": old_manifest.get("package_base_git_sha"),
    "executed_package_working_tree_sha256": None,
    "pre_run_manifest_package_working_tree_sha256": old_manifest.get("package_working_tree_sha256"),
    "execution_environment": {
        "host": "macOS arm64",
        "R": "4.5 arm64 framework; exact patch release not recorded",
        "TRSC": "0.1.0",
        "ECOSolveR": "0.6.1",
        "gmp": "0.7-5.1",
    },
    "manual_pre_run_fixes": [
        "decoded the haven-labelled Stata state variable with haven::as_factor",
        "added and exported trsc_boundary_readouts after the initial installed package lacked it",
    ],
    "manuscript_sha256": old_manifest.get("manuscript_sha256"),
    "dataset_sha256": old_manifest.get("dataset_sha256"),
    "config_sha256": old_manifest.get("config_sha256"),
    "completed_stage_status_files": {
        path.name: {
            key: (None if pd.isna(value) else value)
            for key, value in pd.read_csv(path).iloc[0].to_dict().items()
        }
        for path in sorted(DIAGNOSTICS.glob("[0-9][0-9]_*_status.csv"))
    },
    "selected_candidates": {
        f"q{q}": pd.read_csv(TABLES / f"q{q}_selection.csv").iloc[0][
            ["selected_candidate", "kappa", "h", "validation_score"]
        ].to_dict()
        for q in (1, 2)
    },
    "inference_status": pd.read_csv(
        OUTPUT / "certificates/inference_tier_status.csv"
    ).to_dict(orient="records"),
    "audit_reconciliation_file": "output/audit_reconciliation.csv",
    "execution_provenance_file": "output/execution_provenance.csv",
}
(OUTPUT / "run_manifest_completed_audit.json").write_text(
    json.dumps(completed_manifest, indent=2, default=lambda x: x.item())
)

plt.style.use("seaborn-v0_8-whitegrid")
wide = pd.read_csv(ROOT / "data/processed/cigsale_wide.csv").set_index("year")
role = pd.read_csv(DIAGNOSTICS / "static_role_matched_trajectory.csv").set_index("year")
q1 = pd.read_csv(TABLES / "q1_dynamic_bounds.csv")
q1 = q1[q1["c"] == q1["c"].min()].set_index("year")

fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(wide.index, wide["California"], color="black", linewidth=2, label="California observed")
ax.plot(role.index, role["synthetic"], color="#2C7FB8", linewidth=1.8, linestyle="--", label="Role-matched static SCM")
ax.fill_between(q1.index, q1["lower"], q1["upper"], color="#31A354", alpha=0.22, label="q=1 plug-in set at R=0")
ax.plot(q1.index, q1["lower"], color="#238B45", linewidth=1)
ax.plot(q1.index, q1["upper"], color="#238B45", linewidth=1)
ax.axvline(1989, color="grey", linestyle=":", linewidth=1.2)
ax.set(xlabel="Year", ylabel="Cigarette sales (packs per capita)", title="California, static SCM, and the q=1 plug-in counterfactual set")
ax.legend(frameon=False, loc="upper right")
fig.tight_layout()
fig.savefig(FIGURES / "figure_main_q1_counterfactual.png", dpi=240)
plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(12, 5.2), sharey=False)
colors = {1: "#1B9E77", 3: "#D95F02", 6: "#7570B3", 12: "#1F78B4"}
for ax, q in zip(axes, (1, 2)):
    dynamic = pd.read_csv(TABLES / f"q{q}_dynamic_bounds.csv")
    for horizon in (1, 3, 6, 12):
        z = dynamic[dynamic["horizon"] == horizon].sort_values("c")
        monotone_lower = np.minimum.accumulate(z["lower"].to_numpy())
        monotone_upper = np.maximum.accumulate(z["upper"].to_numpy())
        width = monotone_upper - monotone_lower
        ax.plot(z["c"], width, marker="o", markersize=3, linewidth=1.7,
                color=colors[horizon], label=f"h={horizon}")
    ax.set_xscale("symlog", linthresh=0.25, base=2)
    ax.set_xlabel("Sensitivity multiplier c")
    ax.set_ylabel("Monotone certified outer width")
    ax.set_title(f"q={q}")
    ax.legend(frameon=False)
fig.suptitle("Dynamic sensitivity (monotone closure of certified endpoint enclosures)")
fig.tight_layout()
fig.savefig(FIGURES / "figure_main_dynamic_sensitivity.png", dpi=240)
plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(12, 5.2))
for q, color in ((1, "#1B9E77"), (2, "#7570B3")):
    primary = fibers[(fibers.q == q) & (fibers.geometry_role == "primary_calibration")]
    diagnostic = fibers[(fibers.q == q) & (fibers.geometry_role == "full_pre_diagnostic")]
    axes[0].plot(primary.year, primary.width, color=color, linewidth=2, label=f"q={q}, primary")
    axes[0].plot(diagnostic.year, diagnostic.width, color=color, linewidth=1.5, linestyle="--", label=f"q={q}, full-pre diagnostic")
axes[0].set(xlabel="Deployment year", ylabel="Evaluation-fiber width", title="Evaluation ambiguity")
axes[0].legend(frameon=False, fontsize=9)

for epsilon, color in ((0, "#1B9E77"), (1, "#D95F02"), (5, "#7570B3")):
    z = modulus[(modulus.geometry_role == "primary_calibration") & (modulus.epsilon == epsilon)]
    axes[1].plot(z.year, z.omega, linewidth=2, color=color, label=fr"primary $\epsilon={epsilon}$")
z = modulus[(modulus.geometry_role == "full_pre_diagnostic") & (modulus.epsilon == 0)]
axes[1].plot(z.year, z.omega, color="black", linewidth=1.5, linestyle="--", label=r"full-pre diagnostic $\epsilon=0$")
axes[1].set(xlabel="Deployment year", ylabel="Deployment modulus", title="Deployment amplification")
axes[1].legend(frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig(FIGURES / "figure_main_evaluation_deployment.png", dpi=240)
plt.close(fig)

print("Created audited summaries, reconciliation manifest, and three figures.")
