# Water Simulation Plan

This file plans the move from the LJ-only argon engine to a realistic liquid water simulation.
It follows the conventions in [PLAN.md](PLAN.md): each decision has a status (**Decided**,
**Proposed**, or **Open**), and changes are dated.

**Status: Decisions D1–D5 answered (2026-11-06); plan revised to match.** Nothing here has been
built yet.

---

## 1. Goal

A classical, fixed-charge simulation of liquid water that reproduces basic experimental behaviour
at room temperature and 1 bar:

- Density close to 1.0 g/cm³ at 298 K.
- Oxygen–oxygen radial distribution function with a first peak near 2.8 Å.
- Self-diffusion coefficient of the right order (about 2 × 10⁻⁵ cm²/s).
- Liquid-state energy close to the experimental value (the cohesive energy of liquid water is about 44 kJ/mol).

These targets come from memory and are **to be verified against published sources** before any
of them is used to judge the code (see Section 7).

Out of scope for this plan: chemical reactions (H and O atoms forming bonds from scratch),
polarizable force fields, and GPU acceleration. The reactive case is covered by Milestone 1c in
PLAN.md and is not planned here.

---

## 2. Decisions

| # | Question | Decision | Consequence |
|---|---|---|---|
| D1 | Units for water | **Real units** (Å, fs, kJ/mol, e) in a new unit system. The LJ core stays in reduced units. | Conversion functions between the two systems are needed and must be tested (W1). |
| D2 | Electrostatics | **Ewald** (full summation). PME later. | Ewald is the first electrostatics code to write and verify (W2). |
| D3 | Molecular rigidity | **Flexible** (harmonic bond and angle terms). | No constraint algorithm is needed. The timestep must be about 0.5 fs, not 2 fs. Replaces the SHAKE/RATTLE milestone. |
| D4 | Thermostat and barostat | **Nosé–Hoover.** | Berendsen can be used for equilibration only, if needed. |
| D5 | Water model | **SPC/Fw** (flexible SPC family). Chosen after D3 and D5 conflicted with SPC/E. | Parameters must come from the SPC/Fw paper, not from SPC/E. SPC/E values are not used. |

**Note on the D3/D5 conflict:** SPC/E was parameterised for rigid molecules. The first answer,
flexible molecules with SPC/E, would have been non-standard, so it was changed to SPC/Fw.

**Relation to PLAN.md:** Flexible molecules are consistent with PLAN.md's rule against constraints.
PLAN.md's rule that electrostatics is not in phase one is not, and it needs a dated change in PLAN.md
before W2 starts.

---

## 3. Milestones

Each milestone ends with tests that pass. Nothing moves on until the tests for the milestone pass.

### W0: Prerequisites (no physics changes)
- D1–D5 are answered (Section 2). Record the electrostatics change in PLAN.md with a date.
- Decide where the water code lives. Proposal: new classes in `MD.Core`, with the LJ-only path unchanged.
- Make sure the existing 37 tests still pass before starting.

### W1: Units and topology
- Add a unit system for real units, with conversion functions tested against the existing reduced-unit code.
- Add per-atom charge and type, and per-molecule topology (bonds, angles).
- Atom types: O (LJ and charge) and H (charge only, no LJ), as in SPC/Fw.
- Exclusions: atoms in the same molecule do not interact through LJ or Coulomb terms. Their interaction comes only from the bond and angle terms. The exact rule must match the SPC/Fw paper.
- **Tests:** unit round trips; topology builds the expected bond and angle counts for a known molecule.

### W2: Coulomb electrostatics (Ewald)
- Real-space sum with a cutoff, reciprocal-space sum, self-energy and charged-system correction.
- Forces as the negative gradient of the energy.
- **Tests:**
  - Two-ion energy matches the analytic Coulomb result at short range, with the Ewald correction.
  - NaCl crystal energy matches the Madelung constant (≈ 1.7476 for the NaCl structure, to be verified).
  - Forces equal the negative gradient of energy (finite differences), as for the LJ code.
  - Energy is independent of the Ewald splitting parameter, within tolerance.

### W3: Flexible bonded terms (D3)
- Harmonic bond stretch for O–H and harmonic angle bend for H–O–H.
- Timestep about 0.5 fs, since O–H vibrations are fast. The exact value is set by the energy-drift test.
- **Tests:**
  - Forces equal the negative gradient of the bond and angle energies (finite differences).
  - An isolated molecule vibrates with the period predicted from its force constants (normal-mode check).
  - Energy is conserved in NVE at the chosen timestep (the drift criterion from PLAN.md Section 8).
  - Total momentum and angular momentum of an isolated molecule are conserved.

### W4: SPC/Fw water
- All parameters (charges, LJ on O, equilibrium bond length and angle, force constants) come from
  the SPC/Fw paper (Section 7). No values are written into the code until they are checked against it.
- Build a box of water molecules at the target density.
- **Tests:** a single molecule's equilibrium geometry matches the paper; net charge is zero;
  the molecule stays bonded over a long NVE run.

### W5: Thermostat and barostat (NVT, then NPT)
- Nosé–Hoover thermostat (D4), with tests against an analytic harmonic oscillator.
- Barostat for NPT, once the thermostat is verified.
- **Tests:** the mean temperature matches the target within error; the volume fluctuations match the NPT theory.

### W6: Water validation
- Equilibrate at 298 K and 1 bar. Measure:
  - Density (target near 1.0 g/cm³).
  - O–O radial distribution function (first peak near 2.8 Å).
  - Self-diffusion coefficient from the mean squared displacement (target about 2 × 10⁻⁵ cm²/s).
  - Potential energy per molecule.
- Compare each with the literature values in Section 7.
- Run enough steps to get a converged average. The required length will be set by the diffusion measurement.

### W7: Visualization of H and O atoms
- The trajectory format already has a bond list per frame, currently empty. Water fills it in, so
  the format does not need to change (see [docs/trajectory-format.md](docs/trajectory-format.md)).
- The viewer colours atoms by element (H white, O red, as usual) and draws the bonds.
- Hydrogen bonds between molecules are a physical analysis. They are computed by the core and
  written to the output, not calculated in Godot (PLAN.md Section 11 rule).
- **Tests:** the viewer reads the water trajectory bit-exactly (extend the conformance test).

---

## 4. What the water work needs from the existing code

- The LJ potential is reused for O–O interactions. Its cutoff needs a switch or shift for water (to decide in W1).
- `ForceCalculator` is an O(N²) reference. Water runs need neighbour lists (PLAN.md Milestone 3) to be practical.
  Plan: test the neighbour-list version against the O(N²) version, as PLAN.md already requires.
- The integrator does not change. Velocity Verlet stays, and the bonded forces (W3) are added to the force sum.
- The trajectory writer and reader need the bond list filled in. Their format version may need a bump if the
  bond semantics change (PLAN.md Section 11).

---

## 5. Risks

- **Units errors.** Mixing real and reduced units is the most likely silent error. Every conversion gets a test.
- **Ewald parameters.** Accuracy depends on the splitting parameter and the k-space cutoff. These need a convergence test.
- **Insufficient sampling.** Water's diffusion is slow, so runs need to be long. The sampling time must be stated with each result.
- **Parameter errors.** SPC/Fw values have to match the paper exactly. Test against the paper's numbers, not from memory.
- **Cost of flexibility.** A 0.5 fs timestep means about four times as many steps as rigid water for the same simulated time. Neighbour lists (Section 4) become more important.

---

## 6. Order of work

1. D1–D5 answered (2026-11-06).
2. Update PLAN.md with the adopted decisions (dated), including the electrostatics change.
3. W0 → W1 → W2 → W3 → W4 → W5 → W6 → W7, each with its tests.

---

## 7. References to collect (before W6)

No literature values have been checked yet. Before any water result is compared with reality, collect:

- Wu, Tepper & Voth (2006), *Flexible simple point-charge water model with improved liquid-state properties*, J. Chem. Phys. 124, 024503, for SPC/Fw parameters and the model's own reported liquid properties (citation from memory; verify).
- Experimental density, RDF and diffusion data for liquid water at 298 K (for example, the Soper neutron-diffraction work, and the standard diffusion measurements).
- The Madelung constant for NaCl, for the Ewald test in W2.

Record each source here with the value used.
