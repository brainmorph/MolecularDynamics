using System;
using System.Collections.Generic;

namespace MD.Core;

/// <summary>
/// Lennard-Jones forces for a molecular system: only atoms whose type has an LJ site take part,
/// and excluded pairs (same molecule, see <see cref="Topology.IsExcluded"/>) are skipped. The
/// potential is truncated and energy-shifted at the cutoff, as for the argon code (PLAN.md Section 5).
///
/// Only pairs of identical LJ parameters are supported: no mixing rule has been decided, and SPC/Fw
/// has a single LJ type (O). A pair of different LJ types is reported as an error.
///
/// Sums over every periodic image within the cutoff (like <see cref="ForceCalculator"/>); for an
/// excluded pair only the minimum-image copy is left out.
/// </summary>
public sealed class MolecularLennardJones : IForceField
{
    private readonly int[] _ljAtoms;
    private readonly LennardJones _potential;

    public MolecularLennardJones(Topology topology, double cutoff)
    {
        Topology = topology ?? throw new ArgumentNullException(nameof(topology));
        var atoms = new List<int>();
        AtomType? ljType = null;
        for (int i = 0; i < topology.AtomCount; i++)
        {
            AtomType t = topology.TypeOf(i);
            if (!t.HasLj) continue;
            if (ljType is null) ljType = t;
            else if (t.LjSigma != ljType.LjSigma || t.LjEpsilon != ljType.LjEpsilon)
            {
                throw new NotSupportedException($"LJ types {ljType.Name} and {t.Name} differ; no mixing rule has been decided.");
            }
            atoms.Add(i);
        }
        if (ljType is null) throw new ArgumentException("No atom type has an LJ site.", nameof(topology));
        _ljAtoms = atoms.ToArray();
        _potential = new LennardJones(ljType.LjSigma, ljType.LjEpsilon, cutoff);
    }

    public Topology Topology { get; }
    public LennardJones Potential => _potential;

    public ForceResult Compute(ReadOnlySpan<Vec3> positions, Box box, Span<Vec3> forces)
    {
        if (positions.Length != Topology.AtomCount) throw new ArgumentException("Topology and positions have different atom counts.");
        if (forces.Length != positions.Length) throw new ArgumentException("Force and position arrays must have the same length.");
        forces.Clear();

        Vec3 L = box.Lengths;
        double rc = _potential.Cutoff;
        double rc2 = _potential.CutoffSquared;
        int nx = (int)Math.Ceiling(rc / L.X + 0.5);
        int ny = (int)Math.Ceiling(rc / L.Y + 0.5);
        int nz = (int)Math.Ceiling(rc / L.Z + 0.5);
        double energy = 0;
        double virial = 0;

        for (int a = 0; a < _ljAtoms.Length; a++)
        {
            int i = _ljAtoms[a];
            Vec3 fi = Vec3.Zero;
            for (int b = a; b < _ljAtoms.Length; b++)
            {
                int j = _ljAtoms[b];
                bool self = i == j;
                bool excluded = !self && Topology.IsExcluded(i, j);
                double weight = self ? 0.5 : 1.0;
                Vec3 d0 = box.MinimumImage(positions[i] - positions[j]);
                Vec3 fij = Vec3.Zero;
                for (int sx = -nx; sx <= nx; sx++)
                {
                    double dx = d0.X + sx * L.X;
                    double dx2 = dx * dx;
                    if (dx2 >= rc2) continue;
                    for (int sy = -ny; sy <= ny; sy++)
                    {
                        double dy = d0.Y + sy * L.Y;
                        double dxy2 = dx2 + dy * dy;
                        if (dxy2 >= rc2) continue;
                        for (int sz = -nz; sz <= nz; sz++)
                        {
                            if (sx == 0 && sy == 0 && sz == 0 && (self || excluded)) continue;
                            double dz = d0.Z + sz * L.Z;
                            double r2 = dxy2 + dz * dz;
                            if (!_potential.Evaluate(r2, out double u, out double fOverR)) continue;
                            energy += weight * u;
                            virial += weight * fOverR * r2;
                            if (!self) fij += fOverR * new Vec3(dx, dy, dz);
                        }
                    }
                }
                if (!self)
                {
                    fi += fij;
                    forces[j] -= fij;
                }
            }
            forces[i] += fi;
        }
        return new ForceResult(energy, virial);
    }
}
