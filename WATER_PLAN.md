# Water Simulation Plan

This file plans the move from the LJ-only argon engine to a realistic liquid water simulation.
It follows the conventions in [PLAN.md](PLAN.md): each decision has a status (**Decided**,
**Proposed**, or **Open**), and changes are dated.

**Status: Decisions D1–D5 answered (2026-11-06). W0, W1 and W2 done (2026-11-06).** Next: W3 (bonded terms).

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

### W0: Prerequisites (no physics changes) — Done (2026-11-06)
- D1–D5 are answered (Section 2). Record the electrostatics change in PLAN.md with a date. *(Done: PLAN.md Sections 3, 4, 5, phase table and decision log.)*
- Decide where the water code lives. *(Decided: new classes in `MD.Core`, with the LJ-only path unchanged.)*
- Make sure the existing 37 tests still pass before starting. *(Done: 37/37 passed before any change.)*

### W1: Units and topology
- Add a unit system for real units, with conversion functions tested against the existing reduced-unit code.
- Add per-atom charge and type, and per-molecule topology (bonds, angles).
- Atom types: O (LJ and charge) and H (charge only, no LJ), as in SPC/Fw.
- Exclusions: atoms in the same molecule do not interact through LJ or Coulomb terms. Their interaction comes only from the bond and angle terms. The exact rule must match the SPC/Fw paper.
- **Tests:** unit round trips; topology builds the expected bond and angle counts for a known molecule.

**Done (2026-11-06).** What was built:
- `RealUnits` (`src/MD.Core/RealUnits.cs`): Å, kJ/mol, g/mol, e. The core time unit is the one these
  imply, 100 fs exactly, so velocity Verlet needs no conversion factor; fs is used only at input/output
  (**Decided 2026-11-06**, see PLAN.md decision log). Also R in kJ/(mol·K), the Coulomb constant
  1/(4πε₀) = 1389.35 kJ·Å/(mol·e²) for W2, pressure (bar) and density (g/cm³) conversions.
- `AtomType`, `MoleculeTemplate`, `BondTerm`, `AngleTerm`, `Topology` (`src/MD.Core/Topology.cs`).
  Charge, mass and LJ parameters are per type; a type with ε = 0 has no LJ site (H in SPC/Fw).
  Named `BondTerm`/`AngleTerm` because the trajectory format already has a `Bond`. Force
  constants are added in W3.
- Exclusion rule implemented: all pairs in the same molecule are excluded from LJ and Coulomb.
  For three-site water this equals excluding 1-2 and 1-3 pairs. **Still to check against the SPC/Fw paper (before W4).**
- Tests (`tests/MD.Tests/RealUnitsTests.cs`, `TopologyTests.cs`, 8 new, 45/45 pass):
  constants against known values; round trips; argon's τ, ε/σ³, m/σ³ and T through `RealUnits`
  equal the `ArgonUnits` values; the same 8-atom argon liquid run in reduced and in real units follows
  the same path (within 1e-7 σ after 200 steps); topology counts, global indices, charges and
  exclusions for a placeholder water-shaped molecule (not SPC/Fw values).
- Agreement limit: since the 2019 SI, 1 g/mol / N_A differs from the atomic mass unit in `ArgonUnits`
  by about 3 × 10⁻¹⁰ (relative), so the cross-checks use a 10⁻⁹ tolerance.

**Not done in W1, carried forward:**
- LJ cutoff treatment for water (Section 4): the energy-shifted cutoff, a force-shifted one, or
  whatever the SPC/Fw paper used. **Open**; should match the paper, so it is settled with W4's parameters.
- `Thermo` and `ThermoState` still assume k_B = 1. In real units their "temperature" is in kJ/mol and must
  go through `RealUnits.EnergyToKelvin`. To fix when the water run path is built (W4/W5).

### W2: Coulomb electrostatics (Ewald)
- Real-space sum with a cutoff, reciprocal-space sum, self-energy and charged-system correction.
- Forces as the negative gradient of the energy.
- **Tests:**
  - Two-ion energy matches the analytic Coulomb result at short range, with the Ewald correction.
  - NaCl crystal energy matches the Madelung constant (≈ 1.7476 for the NaCl structure, to be verified).
  - Forces equal the negative gradient of energy (finite differences), as for the LJ code.
  - Energy is independent of the Ewald splitting parameter, within tolerance.

**Done (2026-11-06).** What was built:
- `Ewald` (`src/MD.Core/Ewald.cs`): real space over all images within the cutoff (correct for
  small boxes, like `ForceCalculator`), reciprocal sum over |k| ≤ k_c, self term, exclusion
  correction (−k q_i q_j erf(αr)/r for each intramolecular pair, minimum image), and the
  neutralising-background term −kπQ²/(2Vα²). Conducting (tinfoil) boundary conditions.
  `Ewald.ForAccuracy(r_c, δ, k)` sets α = √(−ln δ)/r_c and k_c = 2α√(−ln δ).
  The virial is not computed yet (needed for the barostat, W5).
- `SpecialFunctions.Erf/Erfc`: .NET has no error function. Series below 2, continued fraction from 2.
  Measured against the C library over [0, 6]: absolute error ≤ 1e-15, erfc relative ≤ 4e-15 from 2.
- `Topology.ExcludedPairs()`, `Vec3.Cross`.
- Tests (`tests/MD.Tests/EwaldTests.cs`, 8 new, 53/53 pass):
  - NaCl, 1 and 2³ conventional cells: E = −M·N_pairs/d to 1e-11 relative, M = 1.7475645946331821; lattice forces vanish.
  - +1/−1 pair in a 30 Å cube at r = 0.5, 1, 2: energy −k/r − 2πk r²/(3V) and the matching force, to 1e-5
    (the dipole term is up to 3e-4, so it is resolved). This is the "Ewald correction" for the two-ion test.
  - One molecule alone: exclusions leave only the image terms. Energy equals the dipole term to 1 %, and
    the remainder falls as L⁻⁵ (ratio 32 between L = 30 and 60 Å, within 1 %). Without exclusions it would be ≈ −1300 kJ/mol.
  - Forces equal −∇E by central differences, for neutral and charged systems with molecules straddling the box edge.
  - Total energy independent of α (four values; terms change by thousands of kJ/mol, total by < 1e-9 kJ/mol), neutral and charged.
  - Forces sum to zero.
- The Forces = −∇E test takes about 5 s (O(N²·N_k) reference code, many energy calls).

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
- The Madelung constant for NaCl, for the Ewald test in W2. **Recorded (2026-11-06):** M = 1.7475645946331821 (nearest-neighbour distance as the length), from the CARMA numbers page, https://carmamaths.org/resources/numbers/page/main.html; Wikipedia's Madelung constant table gives ±1.747565 for the same convention.

Record each source here with the value used.
