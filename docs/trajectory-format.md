# Trajectory file format, version 1

This file is the only interface between the simulation core and the Godot viewer. Both sides
depend on this spec, not on each other's code. Any change that would break an existing reader
requires a version bump.

- **File extension:** `.traj`
- **Byte order:** little-endian for every number.
- **Floating point:** IEEE 754 binary64 (`f64`).
- **Integers:** `u32` unsigned 32-bit, `i32` signed 32-bit, `i64` signed 64-bit.
- **Strings:** `u32` byte length, followed by that many bytes of UTF-8. No terminator.
- **Layout:** one header, then zero or more frames until end of file. There is no frame count
  in the header, so a run that stops early still leaves a readable file. A reader should read
  frames until end of file, and treat a file that ends partway through a frame as truncated
  (keep the complete frames and report the problem).

## Units

The `units` string in the header names the unit system. Every physical value in the file (box,
positions, energies, masses, times) is in that system. There are two:

- **`"lj-reduced"`** (argon runs): reduced Lennard-Jones units. Length in σ, energy in ε, mass in m,
  time in τ = σ√(m/ε), and k_B = 1.
- **`"real"`** (molecular runs, such as water; added 2026-11-06): length in Å, energy in kJ/mol,
  mass in g/mol (u), time in units of exactly 100 fs (the unit that makes Å, kJ/mol and g/mol
  consistent, so velocity Verlet needs no conversion factor).

The header carries scale factors that say what one unit of each quantity equals:

| Quantity | Field | `"lj-reduced"` (argon) | `"real"` |
|---|---|---|---|
| length | `sigma_angstrom` Å | 3.405 | 1 |
| energy | `epsilon_kelvin` × k_B (Kelvin) | 119.8 | 1/R ≈ 120.272 (1 kJ/mol in K) |
| mass | `mass_amu` u | 39.948 | 1 |
| time | `tau_fs` fs | ≈ 2156.3 | 100 |

The field names come from the argon case; for `"real"` they hold the values above.
Temperature in K = (temperature in file units) × `epsilon_kelvin` in both systems.

A reader that does not know a `units` string should still read the file (the layout is the same)
but should not label values with units it cannot name.

## Header

| Field | Type | Notes |
|---|---|---|
| magic | 8 bytes | ASCII `LJMDTRAJ` |
| version | u32 | `1` |
| units | string | `"lj-reduced"` or `"real"` (see Units) |
| sigma_angstrom | f64 | σ in Å (3.405 for argon) |
| epsilon_kelvin | f64 | ε / k_B in K (119.8 for argon) |
| mass_amu | f64 | m in u (39.948 for argon) |
| tau_fs | f64 | τ in fs (≈ 2156.3 for argon) |
| timestep | f64 | integration timestep, time units |
| initial_kinetic_energy | f64 | KE at step 0, energy units |
| initial_potential_energy | f64 | PE at step 0, energy units |
| initial_total_energy | f64 | E₀ = KE + PE at step 0, energy units. Reference for energy drift, so drift can be measured even if frame 0 is missing. |
| atom_count | u32 | N |
| atoms | N × (string element, f64 mass) | element label and mass (mass units) for atom 0 … N−1 |
| settings_json | string | the run settings, as JSON, so every result records how it was made |

In `"lj-reduced"` files element labels are display labels only: all atoms use the same argon-like
LJ parameters, whatever their label. In `"real"` files they are the chemical elements of the
model's atoms (for water: O, H, H per molecule).

## Frame

| Field | Type | Notes |
|---|---|---|
| marker | 4 bytes | ASCII `FRAM` |
| step | i64 | integration step number |
| time | f64 | simulation time, time units |
| kinetic_energy | f64 | energy units, computed by the core at this step |
| potential_energy | f64 | energy units, computed by the core at this step |
| total_energy | f64 | energy units, KE + PE |
| box | 3 × f64 | box lengths Lx, Ly, Lz, length units. The box spans [0, L) on each axis. |
| positions | N × 3 × f64 | x, y, z of each atom, wrapped into [0, L) |
| images | N × 3 × i32 | how many times each atom has crossed each boundary |
| bond_count | u32 | B |
| bonds | B × (u32, u32) | zero-based atom index pairs. Empty (B = 0) for argon runs. For molecular runs, the bonds of the model's topology (for water: O–H1 and O–H2 of every molecule), written by the core in every frame. Never inferred from distances. |

Frame size in bytes: `4 + 8 + 8 + 24 + 24 + 24·N + 12·N + 4 + 8·B`.

### Wrapped and unwrapped positions

Positions are stored wrapped, so atoms that leave one side of the box appear on the opposite side.
The unwrapped (continuous) position is:

```
unwrapped = position + image * L      (per axis)
```

A viewer showing atoms in the box uses the wrapped position. A viewer drawing paths or
interpolating between frames should interpolate unwrapped positions, otherwise an atom crossing
the boundary appears to jump across the whole box.

Wrapped positions can split a molecule across the box: one H on the far side from its O. To draw
molecules whole, a viewer may move bonded atoms by whole box lengths so each is drawn at the
periodic copy nearest its bonded partner. This is display arithmetic; it does not change any value.
The core keeps molecules whole in unwrapped positions (it builds them whole and then only
integrates), but readers should not rely on that.

## Rules for readers (including the Godot viewer)

- Check the magic and version. Reject versions you do not support.
- Display energies as stored. Never recompute energies, forces, or bonds from positions.
- Energy drift `E(t) − E₀` uses `total_energy` from the frame and `initial_total_energy` from the
  header. This subtraction is display arithmetic on values from the core.

## Debugging

`md dump <file.traj> [--all]` prints the header and frames as text.

## Version history

| Version | Date | Change |
|---|---|---|
| 1 | 2026-11-05 | First version |
| 1 | 2026-11-06 | Added the `"real"` units value and bonds for molecular runs. Same layout, so no version bump: version-1 readers still read these files. |
