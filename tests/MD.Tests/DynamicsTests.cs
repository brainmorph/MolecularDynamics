using System;
using MD.Core;

namespace MD.Tests;

public static class DynamicsTests
{
    /// <summary>
    /// A bound LJ dimer (r stays in [1.05, 1.26], far inside the cutoff) is a smooth system, so the
    /// velocity Verlet energy error must scale as dt². Halving dt should divide the error by 4.
    /// This isolates the integrator from the force jump at the cutoff.
    /// </summary>
    [Test]
    public static void Dimer_EnergyError_ScalesAsDtSquared()
    {
        double previous = double.NaN;
        foreach (double dt in new[] { 0.004, 0.002, 0.001 })
        {
            var s = Helpers.Pair(1.05);
            var sim = new Simulation(s, new VelocityVerlet(new ForceCalculator(new LennardJones()), dt));
            double e0 = sim.Initial.TotalEnergy, maxDev = 0;
            for (long k = 0; k < (long)Math.Round(10 / dt); k++)
            {
                sim.Step();
                maxDev = Math.Max(maxDev, Math.Abs(sim.State.TotalEnergy - e0));
                double r = (s.UnwrappedPosition(1) - s.UnwrappedPosition(0)).Length;
                Assert.True(r < 2.0, "dimer must stay well inside the cutoff for this test to be valid");
            }
            if (!double.IsNaN(previous))
            {
                Assert.Close(4.0, previous / maxDev, 0.1, 0, $"error ratio for dt {2 * dt} → {dt}");
            }
            previous = maxDev;
        }
    }

    /// <summary>
    /// 8-atom LJ fluid (the first-phase run), NVE, dt = 0.005, 100 τ. Total energy must stay within
    /// 10⁻³ of its starting value (observed ≈ 3×10⁻⁴; the residual comes from the force jump at the
    /// cutoff, since the potential is energy-shifted but not force-shifted).
    /// </summary>
    [Test]
    public static void LjFluid8_EnergyConserved()
    {
        CheckEnergyConservation(Helpers.Fcc(2, 1, 1), steps: 20_000, maxRelativeDeviation: 1e-3);
    }

    /// <summary>
    /// 8 atoms at 1/125 of the liquid density (box sides 5x longer, larger than 2·r_c): a dilute gas.
    /// E0 is almost all kinetic and small, so the relative tolerance is looser; the absolute deviation
    /// (≈ 0.01 ε, from close collisions) is similar to the liquid's.
    /// </summary>
    [Test]
    public static void LjGas8_EnergyConserved()
    {
        CheckEnergyConservation(Helpers.Fcc(2, 1, 1, density: 0.8442 / 125), steps: 20_000, maxRelativeDeviation: 3e-3);
    }

    [Test]
    public static void LjFluid32_EnergyConserved()
    {
        CheckEnergyConservation(Helpers.Fcc(2, 2, 2), steps: 5_000, maxRelativeDeviation: 1e-3);
    }

    private static void CheckEnergyConservation(ParticleSystem s, int steps, double maxRelativeDeviation)
    {
        var sim = Helpers.NewSimulation(s, 0.722, 0.005);
        double e0 = sim.Initial.TotalEnergy;
        double firstMean = 0, lastMean = 0;
        int window = steps / 10;
        for (int k = 0; k < steps; k++)
        {
            sim.Step();
            double e = sim.State.TotalEnergy;
            Assert.Close(e0, e, 0, maxRelativeDeviation, $"total energy at step {k + 1}");
            if (k < window) firstMean += e / window;
            if (k >= steps - window) lastMean += e / window;
        }
        // No systematic drift: early and late averages agree to the same tolerance.
        Assert.Close(firstMean, lastMean, 0, maxRelativeDeviation, "mean energy, first vs last 10% of run");
    }

    /// <summary>Total momentum starts at zero and stays there (pairwise forces cancel exactly).</summary>
    [Test]
    public static void Momentum_IsConserved()
    {
        var s = Helpers.Fcc(2, 1, 1);
        var sim = Helpers.NewSimulation(s, 1.2, 0.005);
        for (int k = 0; k < 5000; k++) sim.Step();
        Assert.Close(0, s.TotalMomentum().Length, 1e-11, 0, "|total momentum|");
    }

    /// <summary>
    /// Run forward, reverse velocities, run back the same number of steps: the start state must be
    /// recovered. Uses unwrapped positions, so atoms that crossed the boundary are compared correctly.
    /// 1000 steps (5 τ) is short enough that chaotic growth of rounding error stays below 1e-9.
    /// </summary>
    [Test]
    public static void TimeReversible()
    {
        var s = Helpers.Fcc(2, 1, 1);
        var sim = Helpers.NewSimulation(s, 0.722, 0.005);
        var x0 = new Vec3[s.Count];
        var v0 = (Vec3[])s.Velocities.Clone();
        for (int i = 0; i < s.Count; i++) x0[i] = s.UnwrappedPosition(i);

        bool crossed = false;
        for (int k = 0; k < 1000; k++)
        {
            sim.Step();
            foreach (var im in s.Images) crossed |= im != Int3.Zero;
        }
        Assert.True(crossed, "test setup: at least one atom should cross a boundary");

        sim.ReverseVelocities();
        for (int k = 0; k < 1000; k++) sim.Step();

        for (int i = 0; i < s.Count; i++)
        {
            Assert.Close(0, (s.UnwrappedPosition(i) - x0[i]).Length, 1e-9, 0, $"position of atom {i}");
            Assert.Close(0, (s.Velocities[i] + v0[i]).Length, 1e-9, 0, $"velocity of atom {i}");
        }
    }

    /// <summary>Same seed and settings give bit-identical trajectories.</summary>
    [Test]
    public static void Deterministic()
    {
        var a = Helpers.NewSimulation(Helpers.Fcc(2, 1, 1), 0.722, 0.005, seed: 99);
        var b = Helpers.NewSimulation(Helpers.Fcc(2, 1, 1), 0.722, 0.005, seed: 99);
        for (int k = 0; k < 500; k++) { a.Step(); b.Step(); }
        for (int i = 0; i < a.System.Count; i++)
        {
            Assert.Equal(a.System.Positions[i], b.System.Positions[i], $"position {i}");
        }
        Assert.Equal(a.State.TotalEnergy, b.State.TotalEnergy, "energy");
    }
}
