using System;
using MD.Core;

namespace MD.Tests;

public static class SetupTests
{
    [Test]
    public static void Fcc_HasRequestedCountAndDensity()
    {
        var (box, positions) = Lattice.Fcc(2, 1, 1, 0.8442);
        Assert.Equal(8, positions.Length, "atom count for 2x1x1 cells");
        Assert.Close(0.8442, positions.Length / box.Volume, 0, 1e-14, "number density");
        foreach (var p in positions)
        {
            Assert.True(p.X >= 0 && p.X < box.Lengths.X && p.Y >= 0 && p.Y < box.Lengths.Y && p.Z >= 0 && p.Z < box.Lengths.Z,
                $"lattice site {p} inside box");
        }
    }

    [Test]
    public static void Fcc_NearestNeighbourDistance()
    {
        // Nearest-neighbour distance in FCC is a/√2.
        var s = Helpers.Fcc(3, 3, 3);
        double a = Math.Cbrt(4 / Helpers.DefaultDensity);
        double min = double.MaxValue;
        for (int i = 0; i < s.Count; i++)
        for (int j = i + 1; j < s.Count; j++)
            min = Math.Min(min, s.Box.MinimumImage(s.Positions[i] - s.Positions[j]).Length);
        Assert.Close(a / Math.Sqrt(2), min, 1e-12, 0, "nearest-neighbour distance");
    }

    [Test]
    public static void Velocities_HaveZeroMomentumAndExactTemperature()
    {
        var s = Helpers.Fcc(2, 1, 1);
        VelocityInitializer.Assign(s, 0.722, seed: 42);
        Assert.Close(0, s.TotalMomentum().Length, 1e-14, 0, "|total momentum|");
        Assert.Close(0.722, Thermo.Temperature(Thermo.KineticEnergy(s), s.Count), 1e-14, 1e-14, "temperature");
    }

    [Test]
    public static void Velocities_AreDeterministicForSeed()
    {
        var a = Helpers.Fcc(2, 1, 1);
        var b = Helpers.Fcc(2, 1, 1);
        var c = Helpers.Fcc(2, 1, 1);
        VelocityInitializer.Assign(a, 1.0, seed: 9);
        VelocityInitializer.Assign(b, 1.0, seed: 9);
        VelocityInitializer.Assign(c, 1.0, seed: 10);
        bool same = true, differs = false;
        for (int i = 0; i < a.Count; i++)
        {
            same &= a.Velocities[i] == b.Velocities[i];
            differs |= a.Velocities[i] != c.Velocities[i];
        }
        Assert.True(same, "same seed gives same velocities");
        Assert.True(differs, "different seed gives different velocities");
    }

    [Test]
    public static void Velocities_ZeroTemperatureGivesRest()
    {
        var s = Helpers.Fcc(2, 1, 1);
        VelocityInitializer.Assign(s, 0, seed: 1);
        foreach (var v in s.Velocities) Assert.Equal(Vec3.Zero, v, "velocity at T=0");
    }

    /// <summary>Velocities follow Maxwell–Boltzmann: per-component variance ≈ T (unit mass), mean ≈ 0.</summary>
    [Test]
    public static void Velocities_HaveGaussianSpread()
    {
        var s = Helpers.Fcc(10, 10, 10); // 4000 atoms
        VelocityInitializer.Assign(s, 1.5, seed: 3);
        double sum4 = 0, sum2 = 0;
        int n = 0;
        foreach (var v in s.Velocities)
        {
            for (int k = 0; k < 3; k++)
            {
                sum2 += v[k] * v[k];
                sum4 += Math.Pow(v[k], 4);
                n++;
            }
        }
        double var = sum2 / n;
        // For a Gaussian, <v^4> = 3 <v^2>^2. Checks the shape, not just the scale (which is set exactly).
        Assert.Close(3.0, sum4 / n / (var * var), 0.15, 0, "kurtosis");
    }

    [Test]
    public static void Box_WrapsAndCountsImages()
    {
        var box = new Box(new Vec3(2, 3, 4));
        var (p, im) = box.Wrap(new Vec3(2.5, -0.5, 9.0), Int3.Zero);
        Assert.Close(0.5, p.X, 1e-15, 0, "x");
        Assert.Close(2.5, p.Y, 1e-15, 0, "y");
        Assert.Close(1.0, p.Z, 1e-15, 0, "z");
        Assert.Equal(new Int3(1, -1, 2), im, "images");
        var u = box.Unwrap(p, im);
        Assert.Close(2.5, u.X, 1e-15, 0, "unwrapped x");
        Assert.Close(-0.5, u.Y, 1e-15, 0, "unwrapped y");
        Assert.Close(9.0, u.Z, 1e-15, 0, "unwrapped z");
    }

    [Test]
    public static void Box_WrapTinyNegativeStaysInside()
    {
        var box = new Box(new Vec3(1.7, 1.7, 1.7));
        var (p, im) = box.Wrap(new Vec3(-1e-17, 0, 0), Int3.Zero);
        Assert.True(p.X >= 0 && p.X < 1.7, $"wrapped x {p.X:G17} in [0, L)");
        Assert.True(im.X == -1 || im.X == 0, "image count");
    }

    [Test]
    public static void Box_MinimumImage()
    {
        var box = new Box(new Vec3(10, 10, 10));
        var d = box.MinimumImage(new Vec3(9, -6, 4));
        Assert.Equal(new Vec3(-1, 4, 4), d, "minimum image");
    }

    [Test]
    public static void Thermo_TemperatureAndPressure()
    {
        // Two atoms of unit mass, velocities ±(1,0,0): KE = 1, dof = 3, T = 2/3.
        var s = Helpers.Pair(3.0);
        s.Velocities[0] = new Vec3(1, 0, 0);
        s.Velocities[1] = new Vec3(-1, 0, 0);
        double ke = Thermo.KineticEnergy(s);
        Assert.Close(1.0, ke, 0, 1e-15, "KE");
        Assert.Close(2.0 / 3.0, Thermo.Temperature(ke, 2), 0, 1e-15, "T");
        // Ideal-gas part only (no virial): P = N T / V.
        Assert.Close(2 * (2.0 / 3.0) / 1000.0, Thermo.Pressure(ke, 0, 2, 1000.0), 0, 1e-15, "P");
    }
}
