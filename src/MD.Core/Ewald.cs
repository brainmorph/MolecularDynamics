using System;
using System.Collections.Generic;

namespace MD.Core;

/// <summary>Coulomb energy terms from one Ewald evaluation (energy units of the Coulomb constant).</summary>
public readonly record struct EwaldResult(double Real, double Reciprocal, double Self, double Exclusion, double ChargedSystem)
{
    public double Total => Real + Reciprocal + Self + Exclusion + ChargedSystem;
}

/// <summary>
/// Ewald summation of the Coulomb energy and forces of a periodic system (WATER_PLAN.md D2, W2),
/// with conducting ("tinfoil") boundary conditions:
///
///   E_real  = k Σ_{i&lt;j} Σ_n q_i q_j erfc(α|r_ij + n|)/|r_ij + n|  (images within r_c; for an
///             excluded pair the n = 0 term is left out; an atom's own images count half)
///   E_recip = k (2π/V) Σ_{k≠0, |k|≤k_c} e^(−k²/4α²)/k² |S(k)|²,  S(k) = Σ_j q_j e^(ik·r_j)
///   E_self  = −k α/√π Σ_i q_i²
///   E_excl  = −k Σ_{excluded i&lt;j} q_i q_j erf(α r_ij)/r_ij   (removes what E_recip adds for them)
///   E_net   = −k π Q² / (2Vα²)   (uniform neutralising background when the net charge Q ≠ 0)
///
/// Real space sums over every periodic image within the cutoff, like <see cref="ForceCalculator"/>,
/// so it is correct for boxes smaller than twice the cutoff. Excluded pairs use the minimum image,
/// so a molecule must be smaller than half the box. This is the O(N²·N_k) reference version; PME
/// and neighbour lists come later and are to be tested against it.
///
/// The virial is not computed yet; it is added with the barostat (W5).
/// </summary>
public sealed class Ewald
{
    public Ewald(double alpha, double realCutoff, double reciprocalCutoff, double coulombConstant)
    {
        if (!(alpha > 0) || double.IsInfinity(alpha)) throw new ArgumentOutOfRangeException(nameof(alpha));
        if (!(realCutoff > 0) || double.IsInfinity(realCutoff)) throw new ArgumentOutOfRangeException(nameof(realCutoff));
        if (!(reciprocalCutoff > 0) || double.IsInfinity(reciprocalCutoff)) throw new ArgumentOutOfRangeException(nameof(reciprocalCutoff));
        if (!(coulombConstant > 0) || double.IsInfinity(coulombConstant)) throw new ArgumentOutOfRangeException(nameof(coulombConstant));
        Alpha = alpha;
        RealCutoff = realCutoff;
        ReciprocalCutoff = reciprocalCutoff;
        CoulombConstant = coulombConstant;
    }

    /// <summary>
    /// Parameters for a target relative accuracy δ: α = √(−ln δ)/r_c, so erfc(α r_c) ≈ δ, and
    /// k_c = 2α√(−ln δ), so the largest omitted Gaussian factor e^(−k_c²/4α²) = δ.
    /// </summary>
    public static Ewald ForAccuracy(double realCutoff, double accuracy, double coulombConstant)
    {
        if (!(accuracy > 0 && accuracy < 1)) throw new ArgumentOutOfRangeException(nameof(accuracy));
        double s = Math.Sqrt(-Math.Log(accuracy));
        double alpha = s / realCutoff;
        return new Ewald(alpha, realCutoff, 2.0 * alpha * s, coulombConstant);
    }

    public double Alpha { get; }
    public double RealCutoff { get; }
    public double ReciprocalCutoff { get; }
    public double CoulombConstant { get; }

    /// <summary>
    /// Computes Coulomb forces into <paramref name="forces"/> (overwritten) and returns the energy terms.
    /// Charges and exclusions come from <paramref name="topology"/>.
    /// </summary>
    public EwaldResult Compute(ReadOnlySpan<Vec3> positions, Box box, Topology topology, Span<Vec3> forces)
    {
        ArgumentNullException.ThrowIfNull(box);
        ArgumentNullException.ThrowIfNull(topology);
        int n = positions.Length;
        if (topology.AtomCount != n) throw new ArgumentException("Topology and positions have different atom counts.");
        if (forces.Length != n) throw new ArgumentException("Force and position arrays must have the same length.");
        forces.Clear();

        double[] q = topology.Charges();
        double real = RealSpace(positions, box, topology, q, forces);
        double recip = ReciprocalSpace(positions, box, q, forces);
        double excl = ExclusionCorrection(positions, box, topology, q, forces);

        double sumQ = 0, sumQ2 = 0;
        foreach (double qi in q)
        {
            sumQ += qi;
            sumQ2 += qi * qi;
        }
        double self = -CoulombConstant * Alpha / Math.Sqrt(Math.PI) * sumQ2;
        double net = -CoulombConstant * Math.PI * sumQ * sumQ / (2.0 * box.Volume * Alpha * Alpha);
        return new EwaldResult(real, recip, self, excl, net);
    }

    private double RealSpace(ReadOnlySpan<Vec3> positions, Box box, Topology topology, double[] q, Span<Vec3> forces)
    {
        Vec3 L = box.Lengths;
        double rc2 = RealCutoff * RealCutoff;
        int nx = (int)Math.Ceiling(RealCutoff / L.X + 0.5);
        int ny = (int)Math.Ceiling(RealCutoff / L.Y + 0.5);
        int nz = (int)Math.Ceiling(RealCutoff / L.Z + 0.5);
        double a = Alpha;
        double gaussCoeff = SpecialFunctions.TwoOverSqrtPi * a;
        double energy = 0;
        int count = positions.Length;

        for (int i = 0; i < count; i++)
        {
            if (q[i] == 0) continue;
            Vec3 fi = Vec3.Zero;
            for (int j = i; j < count; j++)
            {
                if (q[j] == 0) continue;
                bool self = i == j;
                bool excluded = !self && topology.IsExcluded(i, j);
                double weight = self ? 0.5 : 1.0;
                double kqq = CoulombConstant * q[i] * q[j];
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
                            if (r2 >= rc2) continue;
                            double r = Math.Sqrt(r2);
                            double erfcTerm = SpecialFunctions.Erfc(a * r) / r;
                            energy += weight * kqq * erfcTerm;
                            if (!self)
                            {
                                // F_i = k q_i q_j [erfc(αr)/r + (2α/√π) e^(−α²r²)] r_vec / r²
                                double fOverR = kqq * (erfcTerm + gaussCoeff * Math.Exp(-a * a * r2)) / r2;
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
        return energy;
    }

    private double ReciprocalSpace(ReadOnlySpan<Vec3> positions, Box box, double[] q, Span<Vec3> forces)
    {
        Vec3 L = box.Lengths;
        double V = box.Volume;
        double kc2 = ReciprocalCutoff * ReciprocalCutoff;
        double twoPi = 2.0 * Math.PI;
        int mx = (int)Math.Floor(ReciprocalCutoff * L.X / twoPi);
        int my = (int)Math.Floor(ReciprocalCutoff * L.Y / twoPi);
        int mz = (int)Math.Floor(ReciprocalCutoff * L.Z / twoPi);
        double inv4a2 = 1.0 / (4.0 * Alpha * Alpha);
        int count = positions.Length;
        var cos = new double[count];
        var sin = new double[count];
        double energy = 0;
        double energyPrefactor = CoulombConstant * twoPi / V;
        double forcePrefactor = 2.0 * energyPrefactor; // 4πk/V

        for (int hx = -mx; hx <= mx; hx++)
        {
            for (int hy = -my; hy <= my; hy++)
            {
                for (int hz = -mz; hz <= mz; hz++)
                {
                    if (hx == 0 && hy == 0 && hz == 0) continue;
                    var k = new Vec3(twoPi * hx / L.X, twoPi * hy / L.Y, twoPi * hz / L.Z);
                    double k2 = k.LengthSquared;
                    if (k2 > kc2) continue;
                    double ak = Math.Exp(-k2 * inv4a2) / k2;

                    double sRe = 0, sIm = 0;
                    for (int j = 0; j < count; j++)
                    {
                        double phase = Vec3.Dot(k, positions[j]);
                        cos[j] = Math.Cos(phase);
                        sin[j] = Math.Sin(phase);
                        sRe += q[j] * cos[j];
                        sIm += q[j] * sin[j];
                    }
                    energy += energyPrefactor * ak * (sRe * sRe + sIm * sIm);

                    // F_i = (4πk/V) q_i Σ_k A(k) k [Re S sin(k·r_i) − Im S cos(k·r_i)]
                    for (int i = 0; i < count; i++)
                    {
                        if (q[i] == 0) continue;
                        double s = forcePrefactor * ak * q[i] * (sRe * sin[i] - sIm * cos[i]);
                        forces[i] += s * k;
                    }
                }
            }
        }
        return energy;
    }

    private double ExclusionCorrection(ReadOnlySpan<Vec3> positions, Box box, Topology topology, double[] q, Span<Vec3> forces)
    {
        double a = Alpha;
        double gaussCoeff = SpecialFunctions.TwoOverSqrtPi * a;
        double energy = 0;
        foreach (var (i, j) in topology.ExcludedPairs())
        {
            if (q[i] == 0 || q[j] == 0) continue;
            double kqq = CoulombConstant * q[i] * q[j];
            Vec3 d = box.MinimumImage(positions[i] - positions[j]);
            double r2 = d.LengthSquared;
            double r = Math.Sqrt(r2);
            double erfTerm = SpecialFunctions.Erf(a * r) / r;
            energy -= kqq * erfTerm;
            // E = −k q_i q_j erf(αr)/r, so F_i = k q_i q_j [(2α/√π) e^(−α²r²) − erf(αr)/r] r_vec / r².
            double fOverR = kqq * (gaussCoeff * Math.Exp(-a * a * r2) - erfTerm) / r2;
            Vec3 fij = fOverR * d;
            forces[i] += fij;
            forces[j] -= fij;
        }
        return energy;
    }
}
