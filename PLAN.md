# Molecular Dynamics Simulator: Design Plan

This file tracks design decisions for the project. Each decision has a status:

- **Decided**: settled, do not revisit without a reason
- **Proposed**: a default we lean toward, open to change
- **Open**: needs an answer before work on the affected area starts

Add new decisions at the bottom of the relevant section. When a decision changes, update it in place and note the date and reason.

---

## Working practices

**Status: Decided (2026-11-05)**

- **Platform: Linux.** Development and testing happen on Linux (bash shell). Paths, scripts, and build steps must work there. Windows-specific or macOS-specific tooling is not assumed.

- **No large downloads.** Do not download large files, datasets, or toolchains without asking first. Small, specific references (for example, a single parameter file) are acceptable when needed.
- **No open-ended task spawning.** Do not start complex or multi-part tasks that can branch into further sub-tasks without a clear boundary. Keep work scoped to the request, and report back before expanding it.

---

## 1. Goals and scope

**Status: Proposed**

- Simulate classical molecular dynamics (MD) of atomic and molecular systems.
- **Priorities, in order: correctness, then realism.** Performance comes after both are verified. Where a faster method trades away correctness or physical realism, we choose correctness.
- Out of scope for now: quantum mechanical methods, reactive force fields, GPU acceleration.

**Open:** What systems matter first (e.g., argon, water, proteins, polymers)? This drives the force field and the required features.

---

## 2. Implementation language and tooling

**Status: Decided (2026-11-05)**

- The core simulation is written in **C#**.
- Use .NET with `double` precision throughout the force and integration code. Single precision is not acceptable for the reference path, since energy drift from rounding would obscure physics errors.
- Keep the core as a plain library with no UI or I/O dependencies, so tests can drive it directly.
- **No third-party libraries.** The project uses only native C# (the .NET base class library) and native Godot. No NuGet packages beyond what the Godot C# template provides, and no external packages in the Godot project.
  - Consequences: no xUnit/NUnit/MSTest, so tests use a small in-repo test runner. Serialization uses `System.Text.Json` (part of the base library), not Newtonsoft. Math, linear algebra, and I/O are hand-written or use the base library.

**Run control (Decided 2026-11-05):** A command-line program reads a JSON settings file (atom count, temperature, timestep, output path, and so on) and writes a binary trajectory. It is a separate project from the core library, and it has no UI dependencies.

**Godot version (Decided 2026-11-05):** Godot 4.7 stable, Mono build (needed for C#). The binary is already installed locally at `/home/purple-eden/Desktop/Godot_v4.7-stable_mono_linux_x86_64`, so nothing needs to be downloaded.

---

## 3. Core simulation model

**Status: Proposed**

- **System representation:** particles with position, velocity, mass, type/charge, and a bond topology list.
- **Periodic boundary conditions (Decided 2026-11-05):** Atoms that leave one side re-enter from the opposite side. Orthorhombic box first; triclinic later if needed. Free-space and walled setups are not planned.
- **Atom count for first run (Decided 2026-11-05):** 8 atoms (two FCC cells). Chosen over 5 because 5 does not fit the grid layout. Atom count is a parameter, so larger runs can follow once the 8-atom checks pass. Eight atoms is too few for meaningful averages or clusters, so results at this size test correctness only.
- **Starting layout (Decided 2026-11-05):** Atoms start on a grid (FCC lattice), then receive random velocities drawn for a chosen starting temperature. Total momentum is removed so the box does not drift. Random placement is not used, since overlapping atoms make the first steps unstable.
- **Units (Decided 2026-11-05):** Reduced LJ units. Length in σ, energy in ε, mass in m, with k_B = 1. Used everywhere in the core. Conversion to real argon units (Å, fs, K) happens only at input/output, and the conversion is a documented, tested function.
- **Units for molecular runs (Decided 2026-11-06, WATER_PLAN.md D1):** Water uses real units (`RealUnits`): Å, kJ/mol, g/mol, e, with a core time unit of 100 fs so the integrator needs no factors; times are given in fs at input/output. LJ runs keep reduced units. The two systems are never mixed in one simulation.
- **Ensembles (Decided 2026-11-05):** NVE only (fixed energy, no thermostat) for the first phase, so the drift check is a clean correctness test. NVT and NPT thermostats/barostats come after NVE is verified.

---

## 4. Integrator

**Status: Decided for Milestone 1 (velocity Verlet). Timestep and constraints still Open.**

- Velocity Verlet as the default. It is time-reversible and symplectic, which gives good energy conservation.
- **Timestep (Decided 2026-11-05):** Default 0.005 τ in reduced LJ units. It stays a user input, so it can be changed for tests (e.g., 0.001 τ as a reference).
- Real-unit timestep guidance for all-atom systems is deferred until bonded systems are in scope.

**Decided (2026-11-05):** No constraint algorithms (SHAKE/LINCS) in the first phase. They are only relevant once bonded molecules are in scope.

**Decided (2026-11-06):** Water is flexible (WATER_PLAN.md D3), so constraints are still not needed. Water's timestep is about 0.5 fs, set by the energy-drift test in W3.

---

## 5. Force field and interactions

**Status: Open** (Lennard-Jones pair potential is Decided for Milestone 1)

**Physics caveat for Milestone 1:** A Lennard-Jones potential is a van der Waals (dispersion and repulsion) term. It has no bonds, no charges, and no angular preference. Free H and O atoms interacting only through LJ will form loose, transient clusters. They will not form covalent O–H bonds or a stable H₂O geometry. Bond-capable models may be added later as separate, explicitly labeled phases. See Milestone 1 below.

Initial plan:

1. Lennard-Jones 12-6 pair potential, **cut off at 2.5 σ with energy shift (Decided 2026-11-05)**, so the potential reaches zero at the cutoff and total energy stays conserved.
2. Coulomb electrostatics. Options: plain cutoff (simple, inaccurate), reaction field, or Ewald/PME (accurate, more work).
3. Bonded terms (bonds, angles, dihedrals) for molecular systems.

**Open (not blocking the first phase):**
- Force field file format. Not needed while parameters are hard-coded argon-like values.
- Electrostatics: not in the first phase, since LJ-only has no charges.

**Changed (2026-11-06):** Electrostatics is now in scope for the water phase. Method: Ewald summation (WATER_PLAN.md D2), PME later. The LJ-only argon path is unchanged and has no charges.

---

## 6. Neighbor lists

**Status: Proposed**

- Use a cell list to bin particles, then a Verlet neighbor list with a skin distance (e.g., 2 Å) rebuilt when any particle has moved more than half the skin.
- Build the neighbor list only after the O(N²) reference implementation exists, and use the O(N²) version to test the optimized one.

---

## 7. Input and output

**Status: Open**

- **Input (Decided 2026-11-05):** A JSON settings file (see Run control in Section 11, or the Section 7 run-control note). Starting coordinates are generated (FCC lattice), so no structure file is needed in the first phase.
- **Output:**
  - Trajectory: binary format (Decided 2026-11-05), with the spec documented in the repo.
  - Log: energies, temperature, pressure, and timing every N steps.
- **Checkpointing (Decided 2026-11-05):** Not in the first phase. Added when runs get long enough to need resuming.

---

## 8. Validation and testing

**Status: Proposed**

- **Energy conservation test:** NVE run should conserve total energy to within a stated tolerance.
- **Two-particle test:** compare force and energy against an analytic result.
- **Reference comparison (partly done, 2026-11-06):** LAMMPS and GROMACS are not installed, so the check uses published data instead. `validation/LjReference` compares against Johnson, Zollweg & Gubbins (1993) Table 2 at three state points, and against NIST argon saturation data. Both agree with the engine within statistical error. Real argon is not exactly LJ, so the argon comparison is a model comparison, not a pass/fail test. A direct check against LAMMPS or GROMACS is still open.
- **Time-reversibility test:** run forward, reverse velocities, run back, and check that the initial state is recovered.
- Every new force term gets a finite-difference check that forces equal the negative gradient of energy.

---

## 9. Parallelism

**Status: Open**

- **Single-threaded for the first phase (Decided 2026-11-05).**
- Likely next step: shared-memory parallelism over particles (threads or OpenMP/Rayon).
- Distributed memory (MPI) only if a target system size requires it.

---

## 10. Analysis and visualization

**Status: Open**

- **Initial analysis (Decided 2026-11-05):** energy and temperature time series only. g(r) and mean squared displacement are deferred.
- Visualization: see Section 11. The primary viewer is a Godot project. External tools (VMD, OVITO) remain useful for checking trajectories.

---

## 11. Visualization (Godot)

**Status: Decided for direction (Godot, 3D, separate from core). Format and tooling Open.**

**Priority:** Second. Visualization is built only after the core is verified. A visualization of incorrect physics is worse than none, so the viewer must never be used to validate the simulation.

**Separation rules (Decided):**

- The Godot project lives in its own directory and is not part of the C# core solution.
- The Godot project does not reference the core assemblies, does not call the core API, and does not run simulation logic.
- The only interface is a **documented trajectory file format** written by the core and read by Godot. Both sides depend on the format spec, not on each other's code.
- The format spec lives in the repo (e.g., `docs/trajectory-format.md`) and is versioned. Breaking changes require a version bump.
- Godot may compute display-only things (camera, colors, interpolation between frames). It must not compute or alter physical quantities such as energies, forces, or bond status.

**Trajectory format requirements (Open, but must include):**

- Per frame: time, positions for every atom, and atom identities (element, index).
- Bond/topology information per frame, so bonds forming or breaking can be shown. Without this, a reaction cannot be displayed faithfully.
- Per frame: kinetic energy, potential energy, and total energy, all computed by the core at the same time step as the positions. Godot must not recompute these from positions.
- Header: the energy at the start of the run (or the core writes the initial total energy in the header so drift can be measured without relying on frame 0 being present).
- Box dimensions when periodic boundaries are used, so atoms can be wrapped or unwrapped consistently.
- Units stated explicitly in the file header.

**Camera controls (Decided):**

- Free-fly camera: **W/A/S/D** to move forward/left/back/right, **mouse** to look around.
- **Esc** opens a menu, and pressing Esc again closes it.
  - The menu has **Restart** (restarts the simulation playback from the first frame) and **Quit** (closes the application).
  - While the menu is open, the cursor is released and the camera does not move.
  - Closing the menu recaptures the cursor.
- **Q/E** move down/up, and **Shift** moves faster (Decided 2026-11-05).

**Diagnostics overlay (F1) (Decided):**

- **F1** toggles a diagnostics panel. It is hidden by default. Pressing F1 again closes it.
- Shows:
  - Framerate (frames per second) and frame time in milliseconds, from Godot.
  - Number of atoms on screen (inside the camera's view frustum), from Godot's rendering data. This is display information, not physics.
  - Current kinetic, potential, and total energy, read from the trajectory frame.
  - Energy drift: total energy minus the starting total energy, shown as an absolute value and as a percentage of the starting value.
  - Simulation time and frame index.
- **Ownership of physics numbers:** the core computes all energies. Godot only displays them. Computing drift (`E(t) - E(0)`) is display arithmetic on core values, and the panel must label it as such.
- **Why this matters:** in a correct NVE run, total energy should stay nearly constant. A steady drift means the integrator, timestep, or force code is wrong. The overlay is a quick visual check, but the core's automated tests remain the authority on correctness.
- The panel is a diagnostic, so it should not change camera, playback, or simulation state.
- Movement is relative to the camera's view direction, so flying "forward" follows where you look.
- Mouse look is captured (cursor hidden and locked) while flying, so the view can rotate without the cursor leaving the window.

**Open:**

- **File format (Decided 2026-11-05): binary.** Compact and fast to read. Downside: files can't be checked by eye, so the format spec and a small dump tool (in the core) are needed for debugging. Text export remains an option to add later.
- **Viewer language (Decided 2026-11-05): C#.** The viewer must not reference the core assemblies. Reading the binary trajectory uses only the format spec, not core code. Godot version (4.x) still open.
- Rendering approach for atoms and bonds (instanced meshes are the usual choice for many atoms).

---

## Milestones

### Phase order

**Status: Proposed (2026-11-05).** Order follows the priorities: verify correctness first, then add realism, then speed.

| Phase | Content | Status |
|---|---|---|
| 1 | Milestone 1a: LJ-only engine validation (NVE, velocity Verlet, periodic image-sum forces, trajectory format v1, CLI, tests) | Done (2026-11-05) |
| 2 | Godot viewer: reads trajectory format v1, free-fly camera, Esc menu, F1 diagnostics (Section 11) | Built (2026-11-05); needs hands-on check of controls |
| 3 | Cell list + Verlet neighbour list, tested against the O(N²) reference; larger LJ runs | Not started |
| 4 | Analysis: g(r) and mean squared displacement, checked against published LJ data | Not started |
| 5 | NVT thermostat (then NPT), after NVE is verified | Not started |
| 6 | Water (Milestone 1b), planned in WATER_PLAN.md: flexible SPC/Fw, Ewald electrostatics, real units, Nosé–Hoover. Replaces the earlier rigid SPC/E plan (changed 2026-11-06) | W0 and W1 done (2026-11-06) |

Open items that can change the order: the cutoff treatment (energy-shifted, as decided, or force-shifted), and the LAMMPS/GROMACS reference comparison (Section 8), which needs software that is not installed.

### Milestone 1: Free atoms, Lennard-Jones, velocity Verlet

**Stated goal:** Simulate H and O atoms floating freely, interacting through Lennard-Jones forces, integrated with velocity Verlet, and observe a water molecule forming.

**Status: Goal needs revision before implementation.** As stated, the physics cannot produce the result:

- LJ has no covalent bonding. O–H bonds form only through electronic structure, which a classical fixed-potential model does not capture.
- LJ has no angular term, so it cannot produce the 104.5° H–O–H geometry of water.
- Forming H₂O from free atoms also requires a bond-breaking/forming model (reactive force field such as ReaxFF) or a quantum calculation. Neither is in scope.

**Options to decide between:**

1. **Milestone 1a (recommended first step):** Validate the engine with LJ-only systems of free atoms (e.g., argon or a generic LJ fluid). Checks: energy conservation, time reversibility, two-body force/energy against analytic values. This is needed under any option.
2. **Milestone 1b (realistic water):** Simulate water as rigid or flexible H₂O molecules with fixed O–H bonds and H–O–H angle, plus LJ on O and Coulomb charges on all sites (e.g., the SPC/E model). This is the standard way to get realistic water behavior. Requires the Coulomb method decision (Section 5).
3. **Milestone 1c (reactive):** Start from free H and O atoms and watch bonds form. Requires a reactive force field or quantum method, which is a much larger project.

**Requested by user:** Free H and O atoms with LJ and Verlet, producing water. Not yet confirmed which of the options above matches the intent.

**Revised (2026-11-05):** Starting point is LJ-only. No reactive force field, bond terms, or Coulomb charges in the first phase. The goal is to see how far a pure Lennard-Jones model can go.

**What LJ-only can and cannot produce:**

- Can: gas, liquid, and solid-like phases; droplet and cluster formation; realistic-looking noble-gas behavior (argon-like); correct energy conservation and diffusion.
- Cannot: covalent O–H bonds, bond breaking, the 104.5° water angle, or any molecule with a fixed shape.

**Decided (2026-11-05):** All atoms use argon-like LJ parameters. H and O are treated as argon-like balls for now, so results are checkable against known argon behavior. Any claim that LJ-only forms water would be incorrect and will not be made in the plan or the viewer.

---

## Decision log

| Date | Decision | Reason |
|---|---|---|
| 2026-11-05 | Created this plan document | Start tracking design decisions |
| 2026-11-05 | Core simulation written in C# (.NET, double precision) | User requirement |
| 2026-11-05 | Priorities: correctness first, then realism; performance later | User requirement |
| 2026-11-05 | Velocity Verlet integrator for Milestone 1 | Time-reversible and symplectic; suits energy-conservation tests |
| 2026-11-05 | Milestone 1 "free H/O atoms with LJ forming water" flagged as not achievable with LJ alone | LJ lacks bonds, charges and angular terms; see Milestone 1 |
| 2026-11-05 | Visualization in Godot (3D), second priority after correct simulation | User requirement |
| 2026-11-05 | Godot viewer is fully separate from C# core; communicates only via a documented, versioned trajectory file | User requirement; keeps physics code free of rendering concerns |
| 2026-11-05 | Trajectory format must include per-frame bond topology | Needed to show bond formation/breaking faithfully |
| 2026-11-05 | Viewer may not compute or alter physical quantities | Keeps the visualization from masking or creating physics errors |
| 2026-11-05 | Godot camera is free-fly: WASD + mouse look | User requirement |
| 2026-11-05 | Esc toggles a menu with Restart and Quit; Esc again closes it | User requirement |
| 2026-11-05 | F1 toggles a diagnostics panel: FPS, on-screen atom count, KE/PE/total energy, drift from start | User requirement, to track energy drift |
| 2026-11-05 | F1 pressed again closes the diagnostics panel | User requirement |
| 2026-11-05 | No third-party libraries: native C# (.NET base library) and native Godot only | User requirement; implies in-repo test runner and System.Text.Json |
| 2026-11-05 | No large downloads; no open-ended task spawning; scope work to the request | User requirement; see Working practices |
| 2026-11-05 | Platform is Linux; development and testing target Linux | User requirement |
| 2026-11-05 | Starting point is LJ-only; no reactive force field or bond terms in the first phase | User requirement: see how far LJ can go. Bond-forming physics deferred |
| 2026-11-05 | All atoms use argon-like LJ parameters | Well-known reference behavior to check against |
| 2026-11-05 | Reduced LJ units (σ, ε, m, k_B = 1) inside the core; real-unit conversion only at I/O | Standard for LJ testing; easy to check by hand |
| 2026-11-05 | Periodic boundary conditions (wrap-around); no free-space or wall setups | Bulk-like behavior with no artificial edge forces |
| 2026-11-05 | LJ cutoff at 2.5 σ with energy shift | Standard LJ choice; keeps total energy conserved |
| 2026-11-05 | NVE only (no thermostat) in the first phase | Clean energy-drift test; thermostats come after NVE is verified |
| 2026-11-05 | Default timestep 0.005 τ, user-adjustable | Standard for LJ fluids; accurate enough for energy conservation |
| 2026-11-05 | Start on FCC lattice with random velocities at chosen temperature; remove net momentum | Avoids overlapping atoms; prevents whole-box drift |
| 2026-11-05 | First run uses 8 atoms (two FCC cells); atom count is a parameter | User chose 8 after 5 did not fit the grid layout; correctness checks first |
| 2026-11-05 | Trajectory file format is binary | User choice; a dump tool and format spec are needed because files can't be read by eye |
| 2026-11-05 | Godot viewer written in C#, with no reference to core assemblies | User choice; the separation rule still applies, so the viewer reads only the format spec |
| 2026-11-05 | Camera: Q/E for down/up, Shift for faster movement | User choice |
| 2026-11-05 | No checkpointing in the first phase | 8-atom runs are short; add when needed |
| 2026-11-05 | Single-threaded in the first phase | Correctness before speed |
| 2026-11-05 | First-phase analysis: energy and temperature time series only | Needed for drift check; g(r) and MSD deferred |
| 2026-11-05 | Run control: command-line program reading a JSON settings file | Repeatable runs; settings recorded with each result |
| 2026-11-05 | Godot 4.7 stable, Mono build, installed locally | User choice; no download needed |
| 2026-11-05 | Core writes KE, PE, and total energy per frame; Godot displays them and never recomputes physics | Keeps energy correctness in the core, consistent with the separation rules |
| 2026-11-05 | Reference force code sums over all periodic images within the cutoff (not minimum image) | The 8-atom box (3.36 × 1.68 × 1.68 σ at ρ = 0.8442) is smaller than 2 × 2.5 σ, where minimum image gives wrong forces. Tested against a direct FCC lattice sum and against minimum image on a large box |
| 2026-11-05 | Trajectory format version 1 specified in `docs/trajectory-format.md` | Little-endian, f64 positions, per-frame energies, image counts, bond list (empty for now), settings JSON in header |
| 2026-11-05 | Projects target .NET 8 (`net8.0`) | LTS and installed locally; .NET 10 SDK builds it too |
| 2026-11-05 | Temperature uses 3N − 3 degrees of freedom; pressure from the virial, no tail correction | Net momentum is removed and conserved; potential is truncated and shifted |
| 2026-11-05 | Observed: energy-shifted cutoff causes small random energy kicks (force jumps at r_c). 8 atoms, dt 0.005, 100 τ: max \|E − E0\| ≈ 3 × 10⁻⁴ of \|E0\|, no systematic drift. Integrator itself verified second order on a bound dimer | A force-shifted potential would remove the kicks if tighter conservation is needed; not changed, since the energy shift is a decided item |
| 2026-11-05 | Viewer renders atoms as one instanced MultiMesh of spheres (radius 0.3σ, 1 unit = 1σ) plus a box outline | Resolves the Section 11 rendering question for now; one draw call for all atoms |
| 2026-11-05 | Viewer plays 30 trajectory frames/s, stops on the last frame, pauses while the Esc menu is open; positions interpolated (unwrapped) for display, energies shown from the current frame | Display-only choices; not specified by the user, open to change |
| 2026-11-05 | Core tests compile the viewer's trajectory reader source to check it against the core writer | Catches format disagreements; the viewer still never references core assemblies |
| 2026-11-06 | Reality check uses published data (Johnson et al. 1993 Table 2; NIST argon saturation) in `validation/LjReference`, outside the core solution | LAMMPS/GROMACS are not installed, and PLAN forbids large downloads. The check tests the core, not the viewer |
| 2026-11-06 | Published-data comparisons hold temperature during production (NVT average); the Johnson check stays NVE and interpolates the table to the measured T | Under NVE the measured T drifted by up to ±0.02 at low density, which shifts pressure enough to hide the answer |
| 2026-11-05 | Default 8-atom example is a dilute gas: density 0.0067536 (1/125 of 0.8442), box sides 5× longer (16.8 × 8.4 × 8.4 σ). Dense-liquid setup kept as `examples/argon-8-liquid.json` | User request: the liquid box was too small to watch. The gas box exceeds 2 × cutoff. Energy check: max \|E − E0\| ≈ 0.011 ε (0.15 % of the small, mostly kinetic E0) at dt 0.005, shrinking with dt |
| 2026-11-06 | Default example density reduced 10× in each dimension's box length (density 0.0000067536, box 168 × 84 × 84 σ, i.e. 1000× volume) | User request: atoms were crossing the walls too often to follow. Side effect: no pair comes within the cutoff, so PE = 0 and the atoms fly straight. Not a useful test of the force code; the liquid example and the tests are the physics checks |
| 2026-11-06 | Default example starts as a cold, compact 8-atom FCC cluster at the LJ minimum spacing (density 1.0, nearest neighbours 2^(1/6) σ, T = 0.05), centred in a 168 × 84 × 84 σ box via a new optional `box` setting | User request: start close to equilibrium so low-energy behaviour is visible, and keep the box large so the cluster never meets its images. Supersedes the dilute-gas example above. `density` is only the lattice spacing when `box` is set |
| 2026-11-06 | Water phase adopted per WATER_PLAN.md D1–D5: real units, Ewald, flexible molecules, Nosé–Hoover, SPC/Fw | User decisions. Electrostatics enters the plan for the first time; the rigid SPC/E plan in Phase 6 is replaced |
| 2026-11-06 | Real-unit core time unit is 100 fs (Å, kJ/mol, g/mol consistent set); fs only at input/output | Keeps velocity Verlet unchanged, as WATER_PLAN.md Section 4 requires. Tested by running the same argon liquid in reduced and real units (same path to 1e-7 σ over 200 steps) |
| 2026-11-06 | Water code lives in `MD.Core` as new classes (`RealUnits`, `Topology`); the LJ-only path is unchanged | Resolves WATER_PLAN.md W0 proposal |

---

## Open questions (summary)

1. What is the project's primary goal (learning, research, production)?
2. Which language to use?
3. Which systems and force field to support first?
4. Electrostatics method for version 1? *(Answered 2026-11-06: Ewald for water; see Section 5.)*
5. Unit system? *(Answered: reduced LJ units for argon, real units for water; see Section 3.)*
6. Input and output file formats?
7. Do we need constraints (SHAKE/LINCS) in version 1? *(Answered 2026-11-06: no; water is flexible.)*
8. Is parallelism needed in the first version?
