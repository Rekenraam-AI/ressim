//! #116 barrier-jump probe: where does the residual jump live at the line-search stop?
//!
//! The line-search experiment in `spe5_n1_tests.rs` stops at a residual barrier without a trial
//! phase-label change. This probe takes the final failed iterate and direction and scans `tau`
//! down to roughly an ulp, recording the row-scaled residual change, its accumulation part and the
//! remaining face/well part. It also records the focus cells' flash path, `beta`, `sum ln^2 K`,
//! well-control set and incident-face upwind/gravity state, then reruns the scan under tighter
//! flash tolerances.
//!
//! Diagnostic only: nothing is accepted and the production path is unchanged.
use super::api::CompositionalCase;
use super::assembly::AssemblyResult;
use super::flux::{FaceFlux, face_flux};
use super::layout::CompositionalLayout;
use super::newton::{
    FailedLineSearch, NewtonProblem, NewtonReport, assemble_for_test, line_search_trial_for_test,
    solve_newton, take_failed_line_search_for_test, with_line_search_for_test,
};
use super::spe5_m0_branch_tests::restore;
use super::state::CompositionalState;
use crate::fluid::flash::{
    FlashState, PhaseState, flash, flash_for_cell, has_carried_k_for_test,
    with_difference_refinement_for_test, with_equilibrium_tolerance_for_test,
    with_seed_recovery_for_test, with_split_phase_decision_for_test,
};
use crate::fluid::units::bar_to_pa;
use serde_json::{Value, json};
use std::path::PathBuf;

/// The capture and the step used by the line-search comparison.
const DT_DAYS: f64 = 0.04402126141911289;

/// Largest row-scaled residual change and its accumulation / face-plus-well split.
#[derive(Clone, Debug)]
struct RowChange {
    cell: usize,
    row: usize,
    total: f64,
    accumulation: f64,
    faces_wells: f64,
}

/// Per-row scaled residual changes of `trial` relative to `current`.
fn row_changes(
    current: &AssemblyResult,
    trial: &AssemblyResult,
    layout: &CompositionalLayout,
) -> Vec<RowChange> {
    let n = layout.cell_block_size();
    let mut out = Vec::with_capacity(layout.cell_count() * n);
    for cell in 0..layout.cell_count() {
        let difference: Vec<f64> = (0..n)
            .map(|row| trial.residual[cell * n + row] - current.residual[cell * n + row])
            .collect();
        let scaled = current.scaling[cell].scale_residual(&difference);
        let inventory: Vec<f64> = trial.inventories[cell]
            .component_moles
            .iter()
            .zip(&current.inventories[cell].component_moles)
            .map(|(a, b)| a - b)
            .collect();
        let accumulation = current.scaling[cell].scale_residual(&inventory);
        for row in 0..n {
            out.push(RowChange {
                cell,
                row,
                total: scaled[row],
                accumulation: accumulation[row],
                faces_wells: scaled[row] - accumulation[row],
            });
        }
    }
    out
}

/// The largest `limit` rows by absolute total change.
fn top_rows(mut changes: Vec<RowChange>, limit: usize) -> Vec<RowChange> {
    changes.sort_by(|a, b| b.total.abs().total_cmp(&a.total.abs()));
    changes.truncate(limit);
    changes
}

fn rows_json(rows: &[RowChange]) -> Value {
    json!(
        rows.iter()
            .map(|row| json!({
                "cell": row.cell,
                "row": row.row,
                "total": row.total,
                "accumulation": row.accumulation,
                "faces_wells": row.faces_wells,
            }))
            .collect::<Vec<_>>()
    )
}

fn phase_changes(current: &AssemblyResult, trial: &AssemblyResult) -> Vec<usize> {
    trial
        .inventories
        .iter()
        .enumerate()
        .filter(|(cell, inventory)| inventory.phase_state != current.inventories[*cell].phase_state)
        .map(|(cell, _)| cell)
        .collect()
}

fn flash_at(problem: &NewtonProblem<'_>, cell: usize, state: &CompositionalState) -> FlashState {
    let c = state.cell(cell);
    flash_for_cell(
        cell,
        problem.spec,
        bar_to_pa(c.pressure_bar),
        problem.spec.reservoir_temperature_k(),
        &c.overall_composition(),
        None,
    )
    .expect("cell flash")
}

fn face_present(state: &FlashState, liquid: bool) -> bool {
    matches!(
        (state.phase_state, liquid),
        (PhaseState::TwoPhase, _)
            | (PhaseState::SingleLiquid, true)
            | (PhaseState::SingleVapour, false)
    )
}

fn sum_ln_k_sq(state: &FlashState) -> f64 {
    state.k.iter().map(|k| crate::math::ln(*k).powi(2)).sum()
}

/// The flash path that produced a cell's state, and its equilibrium diagnostics.
fn flash_json(
    problem: &NewtonProblem<'_>,
    cell: usize,
    state: &CompositionalState,
    recovery: bool,
) -> Value {
    let c = state.cell(cell);
    let p = bar_to_pa(c.pressure_bar);
    let t = problem.spec.reservoir_temperature_k();
    let z = c.overall_composition();
    let used = flash_for_cell(cell, problem.spec, p, t, &z, None).expect("used flash");
    let standard = with_split_phase_decision_for_test(false, || {
        flash(problem.spec, p, t, &z, None).expect("standard flash")
    });
    let has_hint = has_carried_k_for_test(cell);
    let path = if recovery && has_hint {
        "carried K hint"
    } else if standard.phase_state != used.phase_state {
        "N1"
    } else {
        "standard"
    };
    json!({
        "cell": cell,
        "p_bar": c.pressure_bar,
        "z": z,
        "has_carried_k_hint": has_hint,
        "path": path,
        "phase": format!("{:?}", used.phase_state),
        "beta": used.beta,
        "sum_ln_k_sq": sum_ln_k_sq(&used),
        "equilibrium_residual": used.equilibrium_residual,
        "iterations": used.iterations,
        "standard_phase": format!("{:?}", standard.phase_state),
        "standard_phase_differs": standard.phase_state != used.phase_state,
    })
}

fn well_json(problem: &NewtonProblem<'_>, assembly: &AssemblyResult) -> Value {
    json!(
        problem
            .wells
            .iter()
            .enumerate()
            .map(|(index, well)| {
                let elimination = assembly
                    .well_eliminations
                    .iter()
                    .find(|elimination| elimination.well == index);
                json!({
                    "id": well.id,
                    "control": format!("{:?}", well.control),
                    "bhp": assembly.well_bhp.get(index),
                    "on_bhp_limit": assembly.well_on_bhp_limit.get(index),
                    "active_rate_row": elimination.is_some(),
                    "rate_residual": elimination.map(|e| e.residual),
                    "rate_scale": elimination.map(|e| e.scale),
                    "rate_components": assembly.well_rates.get(index),
                })
            })
            .collect::<Vec<_>>()
    )
}

fn face_json(problem: &NewtonProblem<'_>, state: &CompositionalState, cell: usize) -> Value {
    let mut faces = Vec::new();
    for (index, face) in problem.faces.iter().enumerate() {
        if face.cell_i != cell && face.cell_j != cell {
            continue;
        }
        let i = state.cell(face.cell_i);
        let j = state.cell(face.cell_j);
        let flux: FaceFlux = face_flux(
            problem.spec,
            problem.relperm,
            face.geom_t,
            face.gravity,
            (face.cell_i, i),
            (face.cell_j, j),
        )
        .expect("face flux");
        let fi = flash_at(problem, face.cell_i, state);
        let fj = flash_at(problem, face.cell_j, state);
        let liquid_present = [face_present(&fi, true), face_present(&fj, true)];
        let vapour_present = [face_present(&fi, false), face_present(&fj, false)];
        faces.push(json!({
            "face": index,
            "i": face.cell_i,
            "j": face.cell_j,
            "gravity_enabled": face.gravity.enabled,
            "component_flux": flux.component_moles_per_day,
            "phase_rates": flux.phase_rates_m3_per_day,
            "upstream": flux.upstream.map(|side| format!("{side:?}")),
            "phase_potential_bar": flux.phase_potential_bar,
            "dp_bar": flux.potential_difference_bar,
            "liquid_present": liquid_present,
            "vapour_present": vapour_present,
            "head_one_sided": [
                liquid_present[0] != liquid_present[1],
                vapour_present[0] != vapour_present[1],
            ],
        }));
    }
    json!(faces)
}

fn trial(
    problem: &NewtonProblem<'_>,
    failure: &FailedLineSearch,
    tau: f64,
) -> Option<(CompositionalState, AssemblyResult)> {
    line_search_trial_for_test(problem, failure, tau)
        .map(|(state, _bhp, assembly)| (state, assembly))
}

/// Scan `tau` down to roughly an ulp and summarize where the residual change lives.
fn scan(problem: &NewtonProblem<'_>, failure: &FailedLineSearch) -> (Value, Vec<usize>, f64, f64) {
    let mut rows = Vec::new();
    let mut top_cells = Vec::new();
    let mut min_norm = f64::INFINITY;
    let mut min_tau = f64::NAN;
    for halving in 0..=60 {
        let tau = failure.tau_bound * 0.5f64.powi(halving);
        let Some((_state, trial_assembly)) = trial(problem, failure, tau) else {
            rows.push(json!({"tau": tau, "assembly_failed": true}));
            continue;
        };
        let norm = trial_assembly.scaled_residual_norm(problem.layout);
        let changes = phase_changes(&failure.current_assembly, &trial_assembly);
        let top = top_rows(
            row_changes(&failure.current_assembly, &trial_assembly, problem.layout),
            3,
        );
        if let Some(first) = top.first()
            && !top_cells.contains(&first.cell)
        {
            top_cells.push(first.cell);
        }
        if norm < min_norm {
            min_norm = norm;
            min_tau = tau;
        }
        rows.push(json!({
            "tau": tau,
            "norm": norm,
            "phase_changes": changes,
            "top_rows": rows_json(&top),
        }));
    }
    (json!(rows), top_cells, min_norm, min_tau)
}

fn details(
    problem: &NewtonProblem<'_>,
    failure: &FailedLineSearch,
    taus: &[f64],
    focus_cells: &[usize],
    recovery: bool,
) -> Value {
    let mut out = Vec::new();
    for &tau in taus {
        let Some((state, assembly)) = trial(problem, failure, tau) else {
            out.push(json!({"tau": tau, "assembly_failed": true}));
            continue;
        };
        out.push(json!({
            "tau": tau,
            "norm": assembly.scaled_residual_norm(problem.layout),
            "phase_changes": phase_changes(&failure.current_assembly, &assembly),
            "top_rows": rows_json(&top_rows(
                row_changes(&failure.current_assembly, &assembly, problem.layout),
                8,
            )),
            "wells": well_json(problem, &assembly),
            "cells": focus_cells
                .iter()
                .map(|&cell| json!({
                    "current": flash_json(problem, cell, &failure.current, recovery),
                    "trial": flash_json(problem, cell, &state, recovery),
                    "faces": face_json(problem, &state, cell),
                }))
                .collect::<Vec<_>>(),
        }));
    }
    json!(out)
}

fn tolerance_sweep(problem: &NewtonProblem<'_>, failure: &FailedLineSearch) -> Value {
    let taus: Vec<f64> = (0..=20)
        .map(|halving| failure.tau_bound * 0.5f64.powi(halving))
        .collect();
    let mut rows = Vec::new();
    for tolerance in [1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14] {
        let (current_norm, norms) = with_equilibrium_tolerance_for_test(tolerance, || {
            let current = assemble_for_test(
                problem,
                &failure.current,
                &failure.current_assembly.well_bhp,
            )
            .expect("current assembly under tolerance");
            let current_norm = current.scaled_residual_norm(problem.layout);
            let norms: Vec<f64> = taus
                .iter()
                .filter_map(|&tau| {
                    trial(problem, failure, tau)
                        .map(|(_, assembly)| assembly.scaled_residual_norm(problem.layout))
                })
                .collect();
            (current_norm, norms)
        });
        let min = norms.iter().copied().fold(f64::INFINITY, f64::min);
        let max = norms.iter().copied().fold(0.0f64, f64::max);
        rows.push(json!({
            "tolerance": tolerance,
            "current_norm": current_norm,
            "min": min,
            "max": max,
            "range": max - min,
            "norms": norms,
        }));
    }
    json!(rows)
}

fn probe_mode(case: &CompositionalCase, recovery: bool) -> Value {
    with_difference_refinement_for_test(true, || {
        with_split_phase_decision_for_test(true, || {
            with_seed_recovery_for_test(recovery, || {
                with_line_search_for_test(true, || {
                    case.with_diagnostic_problem(DT_DAYS, |problem, start, options| {
                        let (outcome, report) = solve_newton(problem, start, options);
                        assert!(outcome.is_err(), "the probe needs the failing solve");
                        let failure = take_failed_line_search_for_test()
                            .expect("a failed line search must leave its iterate and direction");
                        probe_failure(problem, &failure, &report, recovery)
                    })
                })
            })
        })
    })
}

fn probe_failure(
    problem: &NewtonProblem<'_>,
    failure: &FailedLineSearch,
    report: &NewtonReport,
    recovery: bool,
) -> Value {
    let current_norm = failure
        .current_assembly
        .scaled_residual_norm(problem.layout);
    let last = report.residual_history.last().copied().unwrap_or(f64::NAN);
    assert_eq!(
        current_norm, last,
        "the captured iterate must be the report's last"
    );
    let (scan_rows, top_cells, min_norm, min_tau) = scan(problem, failure);
    let mut focus_cells = vec![33usize, 39];
    for cell in top_cells {
        if !focus_cells.contains(&cell) && focus_cells.len() < 8 {
            focus_cells.push(cell);
        }
    }
    let mut taus = vec![failure.tau_bound, 0.5 * failure.tau_bound, 1e-3, 1e-6];
    if min_tau.is_finite() && !taus.iter().any(|tau| *tau == min_tau) {
        taus.push(min_tau);
    }
    json!({
        "recovery": recovery,
        "iteration": failure.iteration,
        "tau_bound": failure.tau_bound,
        "current_norm": current_norm,
        "report_last": last,
        "scan_min_norm": min_norm,
        "scan_min_tau": min_tau,
        "scan": scan_rows,
        "focus_cells": focus_cells,
        "details": details(problem, failure, &taus, &focus_cells, recovery),
        "tolerance_sweep": tolerance_sweep(problem, failure),
    })
}

#[test]
#[ignore = "#116 barrier-jump probe at the line-search stop; minutes; --release --ignored --nocapture"]
fn comp_spe5_barrier_jump_probe() {
    let case = restore("stall_n1_6612_836.json");
    let on = probe_mode(&case, true);
    let off = probe_mode(&case, false);
    let out = json!({
        "seed": "stall_n1_6612_836.json",
        "dt_days": DT_DAYS,
        "recovery_on": on,
        "recovery_off": off,
    });

    if let Ok(dir) = std::env::var("SPE5_TRACE_OUT") {
        let path = PathBuf::from(dir).join("barrier_jump_probe.json");
        std::fs::write(
            &path,
            serde_json::to_string_pretty(&out).expect("serialize"),
        )
        .expect("write the probe");
        println!("written {}", path.display());
    }

    for (label, mode) in [("on", &out["recovery_on"]), ("off", &out["recovery_off"])] {
        let current = mode["current_norm"].as_f64().expect("current norm");
        let min = mode["scan_min_norm"].as_f64().expect("scan min");
        let min_tau = mode["scan_min_tau"].as_f64().expect("scan min tau");
        println!(
            "recovery {label}: iteration {}, current {current:.6e}, scan min {min:.6e} \
             at tau {min_tau:.3e}",
            mode["iteration"]
        );
        let scan = mode["scan"].as_array().expect("scan");
        assert!(
            scan.iter()
                .all(|row| row["phase_changes"].as_array().is_none_or(Vec::is_empty)),
            "the final scan must not change a phase label"
        );
        assert!(
            min >= current * 0.99,
            "the final scan contains a material decrease ({min:e} < {current:e})"
        );
        let details = mode["details"].as_array().expect("details");
        let tau_bound = mode["tau_bound"].as_f64().expect("tau bound");
        let full = details
            .iter()
            .find(|detail| detail["tau"].as_f64() == Some(tau_bound))
            .expect("the full-step detail");
        let top = &full["top_rows"][0];
        assert!(
            top["faces_wells"].as_f64().expect("faces").abs()
                > 100.0 * top["accumulation"].as_f64().expect("accumulation").abs(),
            "the jump must be outside the accumulation term: {top}"
        );
        let expected_path = if label == "on" {
            "carried K hint"
        } else {
            "N1"
        };
        for (cell_index, expected_cell) in [(0usize, 33u64), (1, 39)] {
            assert_eq!(
                full["cells"][cell_index]["trial"]["cell"].as_u64(),
                Some(expected_cell)
            );
            let path = full["cells"][cell_index]["trial"]["path"]
                .as_str()
                .expect("flash path");
            assert_eq!(path, expected_path, "cell {expected_cell} flash path");
        }
        let sweep = mode["tolerance_sweep"].as_array().expect("sweep");
        for row in sweep {
            println!(
                "  flash {:.0e}: current {:.6e}, trials {:.6e}..{:.6e} range {:.6e}",
                row["tolerance"].as_f64().unwrap_or(f64::NAN),
                row["current_norm"].as_f64().unwrap_or(f64::NAN),
                row["min"].as_f64().unwrap_or(f64::NAN),
                row["max"].as_f64().unwrap_or(f64::NAN),
                row["range"].as_f64().unwrap_or(f64::NAN),
            );
            assert!(
                row["min"].as_f64().expect("min") > 1e-8,
                "a flash tolerance reached the acceptance test: {row}"
            );
            assert!(
                row["range"].as_f64().expect("range") > 1e-5,
                "the trial roughness disappeared with flash tolerance: {row}"
            );
        }
    }
    println!("{}", serde_json::to_string_pretty(&out).expect("serialize"));
}
