using System;
using MD.Core;

namespace MD.Tests;

public static class ForceTests
{
    private static readonly LennardJones Lj = new();

    /// <summary>Two atoms: energy and force against the closed-form LJ expressions.</summary>
    [Test]
    public static void TwoParticle_MatchesAnalytic()
    {
        double rc = 2.5;
        double shift = 4 * (Math.Pow(rc, -12) - Math.Pow(rc, -6));
        var calc = new ForceCalculator(Lj);

        foreach (double r in new[] { 0.95, 1.0, Math.Pow(2, 1.0 / 6), 1.5, 2.0, 2.49 })
        {
            var s = Helpers.Pair(r);
            var result = calc.Compute(s.Positions, s.Box, s.Forces);
            double u = 4 * (Math.Pow(r, -12) - Math.Pow(r, -6)) - shift;
            double f = 24 * (2 * Math.Pow(r, -13) - Math.Pow(r, -7)); // −dU/dr, positive = repulsive
            Assert.Close(u, result.PotentialEnergy, 1e-13, 1e-12, $"energy at r={r}");
            // Atom 0 is at smaller x, so a repulsive force pushes it towards −x.
            Assert.Close(-f, s.Forces[0].X, 1e-12, 1e-12, $"force x on atom 0 at r={r}");
            Assert.Close(f, s.Forces[1].X, 1e-12, 1e-12, $"force x on atom 1 at r={r}");
            Assert.Close(0, s.Forces[0].Y, 1e-15, 0, "force y");
            Assert.Close(0, s.Forces[0].Z, 1e-15, 0, "force z");
            Assert.Close(f * r, result.Virial, 1e-12, 1e-12, $"virial at r={r}");
        }
    }

    [Test]
    public static void TwoParticle_ForceIsZeroAtMinimum()
    {
        var s = Helpers.Pair(Math.Pow(2, 1.0 / 6));
        new ForceCalculator(Lj).Compute(s.Positions, s.Box, s.Forces);
        Assert.Close(0, s.Forces[0].X, 1e-12, 0, "force at r = 2^(1/6)");
    }

    /// <summary>The shifted energy goes to zero at the cutoff and is exactly zero beyond it.</summary>
    [Test]
    public static void ShiftedEnergy_IsContinuousAtCutoff()
    {
        Lj.Evaluate(2.5 * 2.5 * (1 - 1e-12), out double inside, out _);
        Assert.Close(0, inside, 1e-12, 0, "energy just inside cutoff");
        bool hit = Lj.Evaluate(2.5 * 2.5, out double at, out double f);
        Assert.True(!hit && at == 0 && f == 0, "pair at exactly r_c must not interact");
        var calc = new ForceCalculator(Lj);
        var s = Helpers.Pair(2.6);
        Assert.Close(0, Helpers.Energy(s, calc), 0, 0, "energy beyond cutoff");
    }

    /// <summary>
    /// Forces equal the negative gradient of the energy (central finite differences), in the 8-atom
    /// box where atoms interact with several images of each neighbour and with images of themselves.
    /// </summary>
    [Test]
    public static void Forces_AreNegativeEnergyGradient_SmallBox()
    {
        CheckFiniteDifference(Helpers.Fcc(2, 1, 1, jitter: 0.2, seed: 3));
    }

    [Test]
    public static void Forces_AreNegativeEnergyGradient_LargeBox()
    {
        CheckFiniteDifference(Helpers.Fcc(3, 3, 3, density: 0.7, jitter: 0.2, seed: 5));
    }

    private static void CheckFiniteDifference(ParticleSystem s)
    {
        var calc = new ForceCalculator(Lj);
        calc.Compute(s.Positions, s.Box, s.Forces);
        var analytic = (Vec3[])s.Forces.Clone();
        double h = 1e-6;
        var original = (Vec3[])s.Positions.Clone();

        for (int i = 0; i < s.Count; i++)
        {
            for (int axis = 0; axis < 3; axis++)
            {
                Vec3 e = axis switch { 0 => new Vec3(h, 0, 0), 1 => new Vec3(0, h, 0), _ => new Vec3(0, 0, h) };
                s.Positions[i] = original[i] + e;
                double up = Helpers.Energy(s, calc);
                s.Positions[i] = original[i] - e;
                double down = Helpers.Energy(s, calc);
                s.Positions[i] = original[i];
                double numeric = -(up - down) / (2 * h);
                Assert.Close(numeric, analytic[i][axis], 1e-6, 1e-6, $"force atom {i} axis {axis}");
            }
        }
    }

    /// <summary>Newton's third law: total force is zero.</summary>
    [Test]
    public static void Forces_SumToZero()
    {
        var s = Helpers.Fcc(2, 1, 1, jitter: 0.3, seed: 11);
        new ForceCalculator(Lj).Compute(s.Positions, s.Box, s.Forces);
        Vec3 sum = Vec3.Zero;
        foreach (var f in s.Forces) sum += f;
        Assert.Close(0, sum.Length, 1e-11, 0, "|Σ F|");
    }

    /// <summary>
    /// A perfect FCC crystal is the same infinite structure whichever supercell is used, so the
    /// energy per atom must not depend on the supercell. Also compared against a direct lattice sum
    /// that does not use the force calculator at all.
    /// </summary>
    [Test]
    public static void PerfectFcc_EnergyPerAtom_IndependentOfSupercell()
    {
        double rho = Helpers.DefaultDensity;
        double expected = DirectFccLatticeSum(rho, 2.5);
        var calc = new ForceCalculator(Lj);
        foreach (var (nx, ny, nz) in new[] { (1, 1, 1), (2, 1, 1), (1, 2, 3), (2, 2, 2), (4, 4, 4) })
        {
            var s = Helpers.Fcc(nx, ny, nz, rho);
            double perAtom = Helpers.Energy(s, calc) / s.Count;
            Assert.Close(expected, perAtom, 1e-12, 1e-12, $"PE/N for {nx}x{ny}x{nz} cells");
            // In a perfect crystal every atom sits at a symmetric site, so each force is zero.
            calc.Compute(s.Positions, s.Box, s.Forces);
            foreach (var f in s.Forces) Assert.Close(0, f.Length, 1e-10, 0, $"force in perfect crystal {nx}x{ny}x{nz}");
        }
    }

    /// <summary>½ Σ over FCC lattice points within r_c of an atom at the origin.</summary>
    private static double DirectFccLatticeSum(double density, double rc)
    {
        double a = Math.Cbrt(4 / density);
        int m = (int)Math.Ceiling(rc / a) + 1;
        double shift = 4 * (Math.Pow(rc, -12) - Math.Pow(rc, -6));
        var basis = new[] { (0.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5), (0.0, 0.5, 0.5) };
        double sum = 0;
        for (int i = -m; i <= m; i++)
        for (int j = -m; j <= m; j++)
        for (int k = -m; k <= m; k++)
        foreach (var (bx, by, bz) in basis)
        {
            double x = a * (i + bx), y = a * (j + by), z = a * (k + bz);
            double r2 = x * x + y * y + z * z;
            if (r2 == 0 || r2 >= rc * rc) continue;
            double s6 = 1 / (r2 * r2 * r2);
            sum += 0.5 * (4 * (s6 * s6 - s6) - shift);
        }
        return sum;
    }

    /// <summary>
    /// When every box length exceeds 2·r_c, the minimum-image convention is valid, so a plain
    /// minimum-image calculation must agree with the image sum.
    /// </summary>
    [Test]
    public static void ImageSum_MatchesMinimumImage_WhenBoxIsLarge()
    {
        var s = Helpers.Fcc(5, 5, 5, jitter: 0.3, seed: 17);
        Assert.True(s.Box.MinLength > 2 * Lj.Cutoff, "test setup: box must exceed 2·r_c");
        var calc = new ForceCalculator(Lj);
        var result = calc.Compute(s.Positions, s.Box, s.Forces);

        var f = new Vec3[s.Count];
        double energy = 0, virial = 0;
        for (int i = 0; i < s.Count; i++)
        for (int j = i + 1; j < s.Count; j++)
        {
            Vec3 d = s.Box.MinimumImage(s.Positions[i] - s.Positions[j]);
            if (!Lj.Evaluate(d.LengthSquared, out double u, out double fr)) continue;
            energy += u;
            virial += fr * d.LengthSquared;
            f[i] += fr * d;
            f[j] -= fr * d;
        }
        Assert.Close(energy, result.PotentialEnergy, 1e-10, 1e-12, "energy");
        Assert.Close(virial, result.Virial, 1e-10, 1e-12, "virial");
        for (int i = 0; i < s.Count; i++) Assert.Close(0, (f[i] - s.Forces[i]).Length, 1e-10, 0, $"force {i}");
    }
}
