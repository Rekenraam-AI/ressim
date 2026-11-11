# SPE5 barrier-jump probe: where the line-search stop lives

**Status: diagnostic, test-only, no production behavior change.** Measured on `4085cab8`
(branch `diag/116-barrier-jump`, stacked on #135), release build.

#135's final line-search stop has no trial phase-label change. This probe takes the final failed
iterate and Newton direction and asks where the residual jump sits: inventory (accumulation) or
faces/wells, which cells, which flash path, and whether a tighter flash tolerance changes it.

## Method

- Capture the final failed line search through the test-only `FailedLineSearch` snapshot
  (`compositional/newton.rs`): the iterate, its assembly, the direction, the BHP corrections and
  the attempted `tau`.
- Scan `tau = tau_bound * 2^-k` for `k = 0..60` and assemble each trial. For every cell row record
  the row-scaled residual change, its accumulation part (`dn` from the trial inventory) and the
  remaining faces-plus-wells part.
- At the full step and selected small `tau`, record each focus cell's flash path (carried `K` or
  N1), phase, `beta`, `sum ln^2 K`, the active well-control set, and every incident face's
  upwind side, phase presence and phase potentials.
- Rerun a 21-point `tau` ladder under flash equilibrium tolerances `1e-9..1e-14` through the
  existing failure-B override.

The carried-`K` hint is the attempt-owned one from the failed solve; it is held fixed while the
probe changes `tau` and the flash tolerance.

## Result

At the final iterate:

| Recovery | Iteration | Current norm | Scan minimum | Minimum at `tau` | Phase changes in scan |
|---|---:|---:|---:|---:|---|
| on | 18 | `1.594464e-6` | `1.587411e-6` | `2.980e-8` | none |
| off | 14 | `2.403354e-6` | `2.403354e-6` | `1.863e-9` | none |

The recovery-on minimum is 0.44% below the current iterate at `tau ~ 3e-8`; the recovery-off
minimum is the current iterate at the smallest representable step. Neither is a line-search
decrease at the halvings actually tried, and no trial changes a phase label anywhere in the scan.

### The jump is in the face terms

At the full step, the largest scaled row changes are cell **39** rows 1–5 (recovery off: cell 33
at some `tau`). The accumulation part is tiny and the face/well part carries essentially all of
the change:

| Recovery | Worst full-step row | Total | Accumulation | Faces + wells |
|---|---|---:|---:|---:|
| on | cell 39, row 5 | `1.070329e-5` | `-7.428e-9` | `1.071072e-5` |
| off | cell 39, row 5 | `-8.467983e-6` | `1.596e-8` | `-8.483943e-6` |

The active well set is not the source: `PROD` is on its BHP limit with no rate row, `INJW` has one
active rate row with residual `-4.5e-13` (recovery on) / `1.4e-12` (recovery off), and `INJG` is
shut in. The remaining change is therefore the face terms of cells 33/39.

### Flash path and state at those cells

| Recovery | Cell | Path | Phase | `beta` range over details | `sum ln^2 K` range |
|---|---|---|---:|---:|---:|
| on | 33 | carried `K` hint | TwoPhase | `0.47457..0.47616` | `3.81e-5..3.86e-5` |
| on | 39 | carried `K` hint | TwoPhase | `0.47295..0.47616` | `3.69e-5..3.81e-5` |
| off | 33 | N1 | TwoPhase | `0.47308..0.47729` | `4.16e-5..4.25e-5` |
| off | 39 | N1 | TwoPhase | `0.47427..0.47873` | `3.81e-5..4.01e-5` |

Cell 40 is standard-flash TwoPhase; cell 82 is standard-flash SingleLiquid. On the one-sided
vapour faces (76 for cell 33, 109 for cell 39) the upstream is the first cell. Gravity is on, and
the faces carrying both phases show different liquid/vapour potentials (for example face 97,
`[0.03073, -0.05045] bar`), so the gravity head is active and phase-specific.

### Flash-tolerance sweep

The residual is tolerance-sensitive, but the roughness does not shrink monotonically and no
tolerance reaches `1e-8` on this ladder:

| Tolerance | Recovery on: current / trial min..max / range | Recovery off: current / trial min..max / range |
|---:|---|---|
| `1e-9` | `7.914e-6` / `1.390e-6..2.073e-5` / `1.934e-5` | `3.925e-6` / `3.026e-6..1.797e-5` / `1.494e-5` |
| `1e-10` | `5.967e-6` / **`1.864e-8..2.711e-5`** / `2.709e-5` | `2.403e-6` / `3.461e-6..2.014e-5` / `1.668e-5` |
| `1e-11` | `1.594e-6` / `4.155e-6..2.719e-5` / `2.303e-5` | `2.403e-6` / `4.147e-6..2.014e-5` / `1.600e-5` |
| `1e-12` | `6.845e-6` / `1.486e-6..2.086e-5` / `1.937e-5` | `4.148e-6` / `4.888e-7..1.513e-5` / `1.464e-5` |
| `1e-13` | `5.897e-6` / `1.580e-6..2.535e-5` / `2.377e-5` | `4.148e-6` / `2.241e-6..1.708e-5` / `1.484e-5` |
| `1e-14` | `6.537e-6` / `1.484e-6..1.153e-5` / `1.005e-5` | `3.684e-6` / `2.107e-6..2.063e-5` / `1.853e-5` |

The recovery-on `1e-10` trial at `1.864e-8` is the only near miss; the next tolerance (`1e-12`)
returns `1.486e-6`, so this is not a monotone tolerance fix. The current iterate's own norm also
moves with the tolerance (`1.59e-6` to `7.91e-6`), so the residual evaluation itself is
flash-sensitive at these cells.

## Verdict

**DIAGNOSTIC — the barrier is in the face terms of cells 33/39 on the current branch, with no
phase-label change in any trial. The well rows are not the source. Tightening the flash tolerance
changes the residual but does not remove the roughness; a flash-only fix is not established.**

The frozen-label root-existence check is still the next step: solve with cells 33/39 held on the
current two-phase branch, then with the single-phase branch, and check whether the flash agrees at
the solution. The face-term result says what to inspect if a root exists; if neither branch agrees,
the model discontinuity is the blocker and a continuous formulation or an explicit acceptance rule
is the choice.

## Reproduce

```bash
SPE5_TRACE_OUT=DIR cargo test --release --manifest-path src/lib/ressim/Cargo.toml --lib \
    comp_spe5_barrier_jump_probe -- --ignored --nocapture
```

The probe writes `DIR/barrier_jump_probe.json`. It is test-only: the failure snapshot and flash
path accessor are `cfg(test)`, and production results are unchanged.
