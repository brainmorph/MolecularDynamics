# Molecular Dynamics

Classical molecular dynamics in C#. The first phase is a Lennard-Jones (argon-like) fluid in a
periodic box, integrated with velocity Verlet in the NVE ensemble. Design decisions are in
[PLAN.md](PLAN.md).

Requires the .NET SDK (projects target `net8.0`). No third-party packages.

## Layout

| Path | Contents |
|---|---|
| `src/MD.Core` | Simulation library: LJ potential, periodic image-sum forces, velocity Verlet, FCC setup, units, trajectory I/O. No UI. |
| `src/MD.Cli` | `md` command-line program: runs a simulation from a JSON settings file, dumps trajectories. |
| `tests/MD.Tests` | Validation tests with a small in-repo test runner. |
| `docs/trajectory-format.md` | Binary trajectory format spec (the contract with the viewer). |
| `examples/` | Example settings files. |
| `viewer/` | Godot 4.7 (C#) trajectory viewer. Separate project; reads only the trajectory file. |

## Usage

```bash
dotnet build -c Release

# Run the tests (optional filter argument matches test names)
dotnet run -c Release --project tests/MD.Tests
dotnet run -c Release --project tests/MD.Tests -- Energy

# Run a simulation: writes out/argon-8.traj and out/argon-8.log
dotnet run -c Release --project src/MD.Cli -- run examples/argon-8.json         # cold row of three 16-atom clusters in a 168 x 84 x 84 box
dotnet run -c Release --project src/MD.Cli -- run examples/argon-8-liquid.json  # dense liquid, box 3.4 x 1.7 x 1.7

# Inspect a trajectory
dotnet run -c Release --project src/MD.Cli -- dump out/argon-8.traj
```

Paths in a settings file are relative to that file. Unknown settings keys are rejected.

## Viewer

```bash
GODOT=/home/purple-eden/Desktop/Godot_v4.7-stable_mono_linux_x86_64/Godot_v4.7-stable_mono_linux.x86_64
(cd viewer && dotnet build)
$GODOT --path viewer                          # plays out/argon-8.traj
$GODOT --path viewer -- /path/to/other.traj   # any trajectory

scripts/view-latest.sh            # builds the viewer and opens the newest .traj under out/
scripts/view-latest.sh some/dir   # newest .traj under some/dir (set GODOT to override the binary path)
```

| Key | Action |
|---|---|
| W / A / S / D | move forward / left / back / right (relative to view) |
| Q / E | move down / up |
| Shift | move faster |
| Mouse | look around |
| Esc | open/close menu (Restart, Quit); playback and camera pause while it is open |
| F1 | show/hide diagnostics: FPS, frame time, atoms on screen, KE/PE/E from the file, drift E − E0 |

Display scale is 1 Godot unit = 1σ, atoms drawn with radius 0.3σ. Camera field of view is 60° vertical. Playback runs at 30 trajectory
frames per second and stops on the last frame. Positions between frames are interpolated for
display only; energies always show the current frame's stored values.

## Settings

| Key | Default | Meaning (reduced LJ units) |
|---|---|---|
| `cells` | `[2, 1, 1]` | FCC cells per axis; atoms = 4 × product |
| `density` | `0.8442` | number density ρσ³ |
| `temperature` | `0.722` | starting temperature (velocities drawn, net momentum removed, rescaled exactly) |
| `timestep` | `0.005` | τ |
| `steps` | `20000` | number of steps |
| `seed` | `12345` | random seed for velocities |
| `cutoff` | `2.5` | LJ cutoff, σ (energy-shifted) |
| `trajectoryEvery` | `10` | steps between trajectory frames |
| `logEvery` | `100` | steps between log lines |
| `elements` | `["Ar"]` | labels assigned to atoms in turn (labels only; all atoms are argon-like) |
| `output` | `out/run.traj` | trajectory path |
| `log` | output with `.log` | log path |

## What is validated

- Two-particle energy, force and virial against the analytic LJ expressions.
- Forces equal the negative energy gradient (finite differences), in small and large boxes.
- Perfect-FCC energy per atom is independent of supercell and matches a direct lattice sum.
- Image summation equals minimum image when the box exceeds 2 × cutoff.
- Velocity Verlet energy error scales as dt² (bound dimer).
- NVE energy conservation (8 and 32 atoms), momentum conservation, time reversibility, determinism.
- Unit conversions against known argon values; trajectory round trip and byte layout.
- The viewer's own reader (compiled into the tests) reads core output bit-exactly.

## Comparison with published data

`validation/LjReference` is a separate program (not part of the core solution). It runs the core
engine and compares the results with published data. LAMMPS and GROMACS are not installed, so
those comparisons are still open.

```bash
dotnet run -c Release --project validation/LjReference            # Johnson et al. (1993) LJ table, 3 state points
dotnet run -c Release --project validation/LjReference -- argon   # LJ saturated liquid density vs real argon (NIST)
```

- **Johnson, Zollweg & Gubbins (1993), Table 2.** At ρ* = 0.80 and 0.90 and T* = 0.70 and 0.75 (500 atoms, r_c = 4σ),
  the raw pressure and energy agree with the paper's cut-shifted values within statistical error. The table is
  interpolated to each run's measured temperature. Energy drift over 100,000 NVE steps is below 0.001%.
- **Real argon (NIST WebBook, saturation at 86.485 K).** The LJ pressure equals argon's vapour pressure at
  ρ* ≈ 0.8334, which matches argon's saturated liquid density (1.4004 g/cm³). The statistical uncertainty from
  the pressure error is about ±0.0002 in ρ*. This is a model comparison: LJ is an approximation for argon, so a
  match does not show the LJ fluid is exact.

Limits: the paper used 864 atoms and this check uses 500. The tail corrections assume g(r) = 1 beyond r_c.
The argon run holds temperature during production, which is an NVT average. Plain NVE gave a shifted answer,
because the measured temperature drifted at low density.
