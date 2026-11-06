using System;

namespace MD.Core;

/// <summary>Potential energy and virial from one force evaluation.</summary>
/// <param name="PotentialEnergy">Total potential energy.</param>
/// <param name="Virial">W = Σ over pairs of r_ij · F_ij (used for pressure).</param>
public readonly record struct ForceResult(double PotentialEnergy, double Virial);

/// <summary>
/// Reference O(N²) force calculation that sums over every periodic image within the cutoff.
///
/// This is correct for any box size, including boxes smaller than twice the cutoff, where the
/// minimum-image convention fails (an atom then interacts with several copies of the same
/// neighbour, and with copies of itself). The first-phase 8-atom box is in that regime.
///
/// It is deliberately simple. Faster methods (minimum image, cell lists, neighbour lists) must be
/// tested against it.
/// </summary>
public sealed class ForceCalculator : IForceField
{
    public ForceCalculator(LennardJones potential)
    {
        Potential = potential ?? throw new ArgumentNullException(nameof(potential));
    }

    public LennardJones Potential { get; }

    /// <summary>Computes forces into <paramref name="forces"/> (overwritten) and returns energy and virial.</summary>
    public ForceResult Compute(ReadOnlySpan<Vec3> positions, Box box, Span<Vec3> forces)
    {
        if (forces.Length != positions.Length)
        {
            throw new ArgumentException("Force and position arrays must have the same length.");
        }
        forces.Clear();

        Vec3 L = box.Lengths;
        double rc = Potential.Cutoff;
        double rc2 = Potential.CutoffSquared;
        // After the minimum-image step each component of d0 lies in [-L/2, L/2], so the images
        // that can be within r_c are those with |shift| ≤ ceil(r_c/L + 1/2).
        int nx = (int)Math.Ceiling(rc / L.X + 0.5);
        int ny = (int)Math.Ceiling(rc / L.Y + 0.5);
        int nz = (int)Math.Ceiling(rc / L.Z + 0.5);

        double energy = 0;
        double virial = 0;
        int n = positions.Length;

        for (int i = 0; i < n; i++)
        {
            Vec3 fi = Vec3.Zero;
            for (int j = i; j < n; j++)
            {
                bool self = i == j;
                // Interactions of an atom with its own images are shared by the atom with itself,
                // so each counts half towards energy and virial. Their forces cancel in ±shift pairs.
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
                            if (self && sx == 0 && sy == 0 && sz == 0) continue;
                            double dz = d0.Z + sz * L.Z;
                            double r2 = dxy2 + dz * dz;
                            if (!Potential.Evaluate(r2, out double u, out double fOverR)) continue;
                            energy += weight * u;
                            virial += weight * fOverR * r2;
                            if (!self)
                            {
                                fij += fOverR * new Vec3(dx, dy, dz);
                            }
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
