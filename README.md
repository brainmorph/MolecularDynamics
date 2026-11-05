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

## Usage

```bash
dotnet build -c Release

# Run the tests (optional filter argument matches test names)
dotnet run -c Release --project tests/MD.Tests
dotnet run -c Release --project tests/MD.Tests -- Energy

# Run a simulation: writes out/argon-8.traj and out/argon-8.log
dotnet run -c Release --project src/MD.Cli -- run examples/argon-8.json

# Inspect a trajectory
dotnet run -c Release --project src/MD.Cli -- dump out/argon-8.traj
```

Paths in a settings file are relative to that file. Unknown settings keys are rejected.

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

Not yet done: comparison against an established code (LAMMPS/GROMACS).
