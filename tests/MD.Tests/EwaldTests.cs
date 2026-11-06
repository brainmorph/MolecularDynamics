using System;
using MD.Core;

namespace MD.Tests;

public static class EwaldTests
{
    // Placeholder charges with water's shape (not SPC/Fw values; those come in W4).
    private static readonly AtomType Ox = new("OX", "O", mass: 16.0, charge: -0.8, ljSigma: 3.0, ljEpsilon: 0.6);
    private static readonly AtomType Hx = new("HX", "H", mass: 1.0, charge: 0.4);
    private static readonly MoleculeTemplate WaterShaped = new(
        "X3", new[] { Ox, Hx, Hx },
        new[] { new BondTerm(0, 1), new BondTerm(0, 2) },
        new[] { new AngleTerm(1, 0, 2) });

    private static MoleculeTemplate Ion(string name, double charge) =>
        new(name, new[] { new AtomType(name, name, mass: 1.0, charge: charge) }, Array.Empty<BondTerm>(), Array.Empty<AngleTerm>());

    private static readonly MoleculeTemplate Cation = Ion("P", +1.0);
    private static readonly MoleculeTemplate Anion = Ion("M", -1.0);

    // ---- Error function -------------------------------------------------------------------------

    /// <summary>Reference values from the C library (Python math.erf / math.erfc).</summary>
    [Test]
    public static void Erf_MatchesReferenceValues()
    {
        var table = new (double X, double Erf, double Erfc)[]
        {
            (0.0, 0.0, 1.0),
            (0.1, 0.1124629160182849, 0.887537083981715),
            (0.5, 0.5204998778130465, 0.4795001221869535),
            (1.0, 0.8427007929497149, 0.15729920705028513),
            (1.5, 0.9661051464753108, 0.033894853524689274),
            (2.0, 0.9953222650189527, 0.004677734981047266),
            (2.5, 0.999593047982555, 0.0004069520174449589),
            (2.9, 0.9999589021219005, 4.109787809945886e-05),
            (2.99, 0.9999764743969194, 2.3525603080640195e-05),
            (3.0, 0.9999779095030014, 2.209049699858544e-05),
            (3.01, 0.9999792610363629, 2.0738963637132634e-05),
            (3.5, 0.9999992569016276, 7.430983723414128e-07),
            (4.0, 0.9999999845827421, 1.541725790028002e-08),
            (5.0, 0.9999999999984626, 1.537459794428035e-12),
            (6.0, 1.0, 2.1519736712498913e-17),
            (8.0, 1.0, 1.1224297172982926e-29),
            (10.0, 1.0, 2.088487583762545e-45),
        };
        foreach (var (x, erf, erfc) in table)
        {
            // Stated accuracy (SpecialFunctions): absolute 1e-15; erfc relative 1e-14 from x = 2.
            Assert.Close(erf, SpecialFunctions.Erf(x), 1e-15, 0, $"erf({x})");
            Assert.Close(-erf, SpecialFunctions.Erf(-x), 1e-15, 0, $"erf({-x})");
            if (x < 2) Assert.Close(erfc, SpecialFunctions.Erfc(x), 1e-15, 0, $"erfc({x})");
            else Assert.Close(erfc, SpecialFunctions.Erfc(x), 0, 1e-14, $"erfc({x})");
            Assert.Close(2.0 - erfc, SpecialFunctions.Erfc(-x), 1e-15, 0, $"erfc({-x})");
        }
    }

    /// <summary>erf′(x) = (2/√π)e^(−x²), checked by central differences on both sides of the method switch.</summary>
    [Test]
    public static void Erf_DerivativeMatchesGaussian()
    {
        // Below 2, erfc = 1 − erf carries up to 1e-15 of absolute rounding; divided by 2h that is
        // up to about 1e-10 in the difference quotient.
        const double h = 1e-5;
        for (double x = 0.05; x < 5.5; x += 0.1)
        {
            double expected = SpecialFunctions.TwoOverSqrtPi * Math.Exp(-x * x);
            double fdErfc = (SpecialFunctions.Erfc(x + h) - SpecialFunctions.Erfc(x - h)) / (2 * h);
            Assert.Close(-expected, fdErfc, 1e-10, 1e-8, $"erfc'({x:F2})");
        }
        // No jump where the method changes.
        // Compare the adjacent double below 2, so the true slope (≈ 0.02) contributes only ~1e-17.
        Assert.Close(SpecialFunctions.Erfc(Math.BitDecrement(2.0)), SpecialFunctions.Erfc(2.0), 1e-15, 0, "erfc continuous at the switch (x = 2)");
    }

    // ---- Ewald ----------------------------------------------------------------------------------

    /// <summary>
    /// Rock salt, nearest-neighbour distance d = 1, k = 1: E = −M·N_pairs/d with the Madelung
    /// constant M = 1.7475645946331821 (published value, see WATER_PLAN.md Section 7).
    /// Checked for one and for 2×2×2 conventional cells. Forces vanish by symmetry.
    /// </summary>
    [Test]
    public static void NaCl_EnergyMatchesMadelungConstant()
    {
        const double madelung = 1.7475645946331821;
        foreach (int cells in new[] { 1, 2 })
        {
            const double a = 2.0; // conventional cell, so d = a/2 = 1
            var top = new Topology();
            var positions = new System.Collections.Generic.List<Vec3>();
            var cation = new[] { (0, 0, 0), (0, 1, 1), (1, 0, 1), (1, 1, 0) };
            var anion = new[] { (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1) };
            for (int cx = 0; cx < cells; cx++)
                for (int cy = 0; cy < cells; cy++)
                    for (int cz = 0; cz < cells; cz++)
                    {
                        var origin = a * new Vec3(cx, cy, cz);
                        foreach (var (x, y, z) in cation)
                        {
                            top.AddMolecule(Cation);
                            positions.Add(origin + 0.5 * a * new Vec3(x, y, z));
                        }
                        foreach (var (x, y, z) in anion)
                        {
                            top.AddMolecule(Anion);
                            positions.Add(origin + 0.5 * a * new Vec3(x, y, z));
                        }
                    }
            var box = new Box(new Vec3(cells * a, cells * a, cells * a));
            var ewald = Ewald.ForAccuracy(realCutoff: 3.0, accuracy: 1e-14, coulombConstant: 1.0);
            var forces = new Vec3[positions.Count];
            var result = ewald.Compute(positions.ToArray(), box, top, forces);

            int pairs = positions.Count / 2;
            Assert.Close(-madelung * pairs, result.Total, 0, 1e-11, $"NaCl energy, {cells}³ cells");
            Assert.Close(0.0, result.ChargedSystem, 0, 0, "neutral system has no background term");
            foreach (var f in forces) Assert.True(f.Length < 1e-10, $"lattice force {f.Length:G3} should vanish");
        }
    }

    /// <summary>
    /// A +1/−1 pair at separation r in a cubic box of side L ≫ r. With conducting boundaries, the
    /// periodic energy is the isolated pair energy −k/r plus the dipole term −2πk·M²/(3V), M = r
    /// (the image dipole lattice sum vanishes for a cube). Higher terms are of order r⁴/L⁵ ≈ 1e-6 here.
    /// The dipole term (up to ≈ 3e-4) is resolved by the 1e-5 tolerance.
    /// </summary>
    [Test]
    public static void TwoIons_MatchCoulombWithDipoleCorrection()
    {
        const double L = 30.0;
        var box = new Box(new Vec3(L, L, L));
        var dir = new Vec3(1, 2, 2) / 3.0;
        var ewald = Ewald.ForAccuracy(realCutoff: 12.0, accuracy: 1e-12, coulombConstant: 1.0);
        foreach (double r in new[] { 0.5, 1.0, 2.0 })
        {
            var top = new Topology();
            top.AddMolecule(Cation);
            top.AddMolecule(Anion);
            var c = new Vec3(10, 11, 12);
            var positions = new[] { c, c + r * dir };
            var forces = new Vec3[2];
            double e = ewald.Compute(positions, box, top, forces).Total;
            double dipole = -2.0 * Math.PI * r * r / (3.0 * L * L * L);
            Assert.Close(-1.0 / r + dipole, e, 1e-5, 0, $"pair energy at r = {r}");
            // Attractive: force on the cation points towards the anion, magnitude ≈ 1/r² − 4πr/(3V).
            double expectedF = 1.0 / (r * r) - 4.0 * Math.PI * r / (3.0 * L * L * L);
            Assert.Close(expectedF, Vec3.Dot(forces[0], dir), 1e-5, 0, $"pair force at r = {r}");
        }
    }

    /// <summary>
    /// One molecule alone: its own pairs are excluded, so only the image interactions remain. The
    /// leading one is the dipole term −2πk·M²/(3V); the rest is the next multipole order, which for a
    /// cube falls as L⁻⁵. Checked by comparing L = 30 and L = 60 Å: the remainder must shrink by
    /// 2⁵ = 32 (within 1 %). Without exclusions the energy would be about −1300 kJ/mol.
    /// </summary>
    [Test]
    public static void IsolatedMolecule_OnlyImageTermsRemain()
    {
        double k = RealUnits.CoulombConstant;
        double Remainder(double L)
        {
            var box = new Box(new Vec3(L, L, L));
            var top = new Topology();
            top.AddMolecule(WaterShaped);
            var o = new Vec3(0.5 * L, 0.5 * L, 0.5 * L);
            var positions = new[] { o, o + new Vec3(0.8, 0.6, 0), o + new Vec3(-0.8, 0.6, 0) };
            var result = Ewald.ForAccuracy(realCutoff: 12.0, accuracy: 1e-13, coulombConstant: k)
                .Compute(positions, box, top, new Vec3[3]);
            Vec3 m = Vec3.Zero;
            double[] q = top.Charges();
            for (int i = 0; i < 3; i++) m += q[i] * positions[i];
            double dipole = -2.0 * Math.PI * k * m.LengthSquared / (3.0 * L * L * L);
            Assert.Close(dipole, result.Total, 0, 0.01, $"isolated molecule energy is the dipole term to 1 % (L = {L})");
            return result.Total - dipole;
        }
        double r30 = Remainder(30.0);
        double r60 = Remainder(60.0);
        Assert.Close(32.0, r30 / r60, 0, 0.01, "remainder ratio between L = 30 and L = 60 (L⁻⁵ scaling)");
    }

    /// <summary>Forces equal −∇E (central differences), for a neutral and for a charged system.</summary>
    [Test]
    public static void Forces_AreNegativeEnergyGradient()
    {
        foreach (bool charged in new[] { false, true })
        {
            var (positions, box, top) = RandomWaterShapedSystem(molecules: 4, charged, seed: 11);
            var ewald = Ewald.ForAccuracy(realCutoff: 4.0, accuracy: 1e-14, coulombConstant: RealUnits.CoulombConstant);
            var forces = new Vec3[positions.Length];
            ewald.Compute(positions, box, top, forces);
            var scratch = new Vec3[positions.Length];
            const double h = 1e-5;
            for (int i = 0; i < positions.Length; i++)
            {
                for (int c = 0; c < 3; c++)
                {
                    var p = (Vec3[])positions.Clone();
                    p[i] = Shift(positions[i], c, h);
                    double ePlus = ewald.Compute(p, box, top, scratch).Total;
                    p[i] = Shift(positions[i], c, -h);
                    double eMinus = ewald.Compute(p, box, top, scratch).Total;
                    double fd = -(ePlus - eMinus) / (2 * h);
                    Assert.Close(fd, Component(forces[i], c), 2e-5, 1e-7, $"force on atom {i}, component {c} (charged: {charged})");
                }
            }
        }
    }

    /// <summary>The total energy must not depend on the splitting parameter α once converged.</summary>
    [Test]
    public static void Energy_IsIndependentOfAlpha()
    {
        foreach (bool charged in new[] { false, true })
        {
            var (positions, box, top) = RandomWaterShapedSystem(molecules: 4, charged, seed: 5);
            double k = RealUnits.CoulombConstant;
            var forces = new Vec3[positions.Length];
            var reference = Ewald.ForAccuracy(4.0, 1e-14, k).Compute(positions, box, top, forces);
            foreach (double rc in new[] { 3.0, 5.0, 7.0 })
            {
                var other = Ewald.ForAccuracy(rc, 1e-14, k).Compute(positions, box, top, new Vec3[positions.Length]);
                // The individual terms change by hundreds to thousands of kJ/mol between these α values.
                Assert.True(Math.Abs(other.Reciprocal - reference.Reciprocal) > 100.0, $"reciprocal term should change with α (rc = {rc})");
                Assert.Close(reference.Total, other.Total, 1e-7, 0, $"total energy at rc = {rc} (charged: {charged})");
            }
        }
    }

    /// <summary>Internal forces sum to zero (no net force on the system), including the reciprocal part.</summary>
    [Test]
    public static void Forces_SumToZero()
    {
        foreach (bool charged in new[] { false, true })
        {
            var (positions, box, top) = RandomWaterShapedSystem(molecules: 5, charged, seed: 3);
            var forces = new Vec3[positions.Length];
            Ewald.ForAccuracy(4.0, 1e-12, RealUnits.CoulombConstant).Compute(positions, box, top, forces);
            Vec3 sum = Vec3.Zero;
            double scale = 0;
            foreach (var f in forces)
            {
                sum += f;
                scale = Math.Max(scale, f.Length);
            }
            Assert.True(sum.Length < 1e-10 * scale, $"net force {sum.Length:G3} (largest {scale:G3}, charged: {charged})");
        }
    }

    /// <summary>
    /// Water-shaped molecules (O–H 1 Å, H–O–H ≈ 106°) at random positions and orientations in a
    /// non-cubic box, with an extra +1 ion when <paramref name="charged"/>. Some molecules straddle
    /// the box boundary after wrapping, which exercises minimum-image exclusions.
    /// </summary>
    private static (Vec3[] Positions, Box Box, Topology Topology) RandomWaterShapedSystem(int molecules, bool charged, int seed)
    {
        var rng = new Random(seed);
        var box = new Box(new Vec3(8.0, 9.0, 10.0));
        var top = new Topology();
        var list = new System.Collections.Generic.List<Vec3>();
        for (int m = 0; m < molecules; m++)
        {
            top.AddMolecule(WaterShaped);
            var o = new Vec3(rng.NextDouble() * 8.0, rng.NextDouble() * 9.0, rng.NextDouble() * 10.0);
            var u = RandomUnit(rng);
            var w = Vec3.Cross(u, RandomUnit(rng));
            w /= w.Length;
            double half = 0.5 * 106.0 * Math.PI / 180.0;
            list.Add(o);
            list.Add(o + Math.Cos(half) * u + Math.Sin(half) * w);
            list.Add(o + Math.Cos(half) * u - Math.Sin(half) * w);
        }
        if (charged)
        {
            top.AddMolecule(Cation);
            list.Add(new Vec3(rng.NextDouble() * 8.0, rng.NextDouble() * 9.0, rng.NextDouble() * 10.0));
        }
        var positions = list.ToArray();
        for (int i = 0; i < positions.Length; i++) positions[i] = box.Wrap(positions[i], Int3.Zero).Position;
        return (positions, box, top);
    }

    private static Vec3 RandomUnit(Random rng)
    {
        while (true)
        {
            var v = new Vec3(2 * rng.NextDouble() - 1, 2 * rng.NextDouble() - 1, 2 * rng.NextDouble() - 1);
            double l = v.Length;
            if (l > 0.1 && l <= 1) return v / l;
        }
    }

    private static Vec3 Shift(Vec3 v, int c, double h) =>
        c == 0 ? new Vec3(v.X + h, v.Y, v.Z) : c == 1 ? new Vec3(v.X, v.Y + h, v.Z) : new Vec3(v.X, v.Y, v.Z + h);

    private static double Component(Vec3 v, int c) => c == 0 ? v.X : c == 1 ? v.Y : v.Z;
}
