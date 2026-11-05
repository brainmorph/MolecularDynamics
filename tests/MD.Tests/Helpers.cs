using System;
using MD.Core;

namespace MD.Tests;

internal static class Helpers
{
    public const double DefaultDensity = 0.8442;

    /// <summary>FCC system with unit masses, optionally with positions randomly displaced.</summary>
    public static ParticleSystem Fcc(int nx, int ny, int nz, double density = DefaultDensity, double jitter = 0, int seed = 1)
    {
        var (box, positions) = Lattice.Fcc(nx, ny, nz, density);
        var rng = new Random(seed);
        for (int i = 0; i < positions.Length; i++)
        {
            positions[i] += jitter * new Vec3(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5, rng.NextDouble() - 0.5);
        }
        var masses = new double[positions.Length];
        Array.Fill(masses, 1.0);
        var elements = new string[positions.Length];
        Array.Fill(elements, "Ar");
        return new ParticleSystem(box, positions, masses, elements);
    }

    public static ParticleSystem Pair(double r, double boxLength = 10.0)
    {
        var box = new Box(new Vec3(boxLength, boxLength, boxLength));
        var c = 0.5 * boxLength;
        var positions = new[] { new Vec3(c - 0.5 * r, c, c), new Vec3(c + 0.5 * r, c, c) };
        return new ParticleSystem(box, positions, new[] { 1.0, 1.0 }, new[] { "Ar", "Ar" });
    }

    public static Simulation NewSimulation(ParticleSystem s, double temperature, double dt, int seed = 7, double cutoff = 2.5)
    {
        VelocityInitializer.Assign(s, temperature, seed);
        return new Simulation(s, new VelocityVerlet(new ForceCalculator(new LennardJones(cutoff: cutoff)), dt));
    }

    public static double Energy(ParticleSystem s, ForceCalculator calc) =>
        calc.Compute(s.Positions, s.Box, new Vec3[s.Count]).PotentialEnergy;
}
