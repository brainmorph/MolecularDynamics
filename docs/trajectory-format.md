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

All physical values are in **reduced Lennard-Jones units**: length in σ, energy in ε, mass in m,
time in τ = σ√(m/ε), and k_B = 1. The `units` string in the header is `"lj-reduced"`.

The header carries the scale factors needed to show real argon units:

| Quantity | Reduced value 1 equals |
|---|---|
| length | `sigma_angstrom` Å |
| energy | `epsilon_kelvin` × k_B (Kelvin) |
| mass | `mass_amu` atomic mass units |
| time | `tau_fs` femtoseconds |

Temperature in K = reduced temperature × `epsilon_kelvin`.

## Header

| Field | Type | Notes |
|---|---|---|
| magic | 8 bytes | ASCII `LJMDTRAJ` |
| version | u32 | `1` |
| units | string | `"lj-reduced"` |
| sigma_angstrom | f64 | σ in Å (3.405 for argon) |
| epsilon_kelvin | f64 | ε / k_B in K (119.8 for argon) |
| mass_amu | f64 | m in u (39.948 for argon) |
| tau_fs | f64 | τ in fs (≈ 2156.3 for argon) |
| timestep | f64 | integration timestep, τ |
| initial_kinetic_energy | f64 | KE at step 0, ε |
| initial_potential_energy | f64 | PE at step 0, ε |
| initial_total_energy | f64 | E₀ = KE + PE at step 0, ε. Reference for energy drift, so drift can be measured even if frame 0 is missing. |
| atom_count | u32 | N |
| atoms | N × (string element, f64 mass) | element label and mass (reduced) for atom 0 … N−1 |
| settings_json | string | the run settings, as JSON, so every result records how it was made |

Element labels are display labels. In the first phase all atoms use the same argon-like LJ
parameters, whatever their label.

## Frame

| Field | Type | Notes |
|---|---|---|
| marker | 4 bytes | ASCII `FRAM` |
| step | i64 | integration step number |
| time | f64 | simulation time, τ |
| kinetic_energy | f64 | ε, computed by the core at this step |
| potential_energy | f64 | ε, computed by the core at this step |
| total_energy | f64 | ε, KE + PE |
| box | 3 × f64 | box lengths Lx, Ly, Lz, σ. The box spans [0, L) on each axis. |
| positions | N × 3 × f64 | x, y, z of each atom, wrapped into [0, L) |
| images | N × 3 × i32 | how many times each atom has crossed each boundary |
| bond_count | u32 | B |
| bonds | B × (u32, u32) | zero-based atom index pairs. Always empty (B = 0) in the LJ-only phase. |

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
