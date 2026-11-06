using System;
using MD.Core;

namespace MD.Tests;

public static class RealUnitsTests
{
    // The 2019 SI fixes N_A but not the atomic mass unit, so 1 g/mol / N_A differs from
    // ArgonUnits.AtomicMassUnitKilogram by about 3e-10 (relative). Cross-checks allow for that.
    private const double CrossCheckRelTol = 1e-9;

    private static readonly double ArgonEpsilonKj = RealUnits.KelvinToEnergy(ArgonUnits.EpsilonOverBoltzmannKelvin);

    [Test]
    public static void DerivedConstants_MatchKnownValues()
    {
        Assert.Close(100.0, RealUnits.TimeUnitFemtoseconds, 0, 1e-14, "core time unit in fs");
        Assert.Close(0.00831446261815324, RealUnits.GasConstantKilojoulesPerMoleKelvin, 0, 1e-11, "R in kJ/(mol K)");
        // GROMACS ONE_4PI_EPS0 = 138.935458 kJ nm/(mol e²), i.e. 1389.35458 in Å.
        Assert.Close(1389.35458, RealUnits.CoulombConstant, 0, 1e-8, "1/(4πε₀) in kJ Å/(mol e²)");
        // 1 kJ/(mol Å³) = 1e3 / N_A J / 1e-30 m³ ≈ 16605.39 bar.
        Assert.Close(16605.39, RealUnits.PressureToBar(1.0), 0.01, 0, "kJ/(mol Å³) in bar");
        // 1 (g/mol)/Å³ = 1 / N_A g / 1e-24 cm³ ≈ 1.66054 g/cm³.
        Assert.Close(1.66054, RealUnits.DensityToGramsPerCubicCentimetre(1.0), 1e-5, 0, "(g/mol)/Å³ in g/cm³");
    }

    [Test]
    public static void Conversions_RoundTrip()
    {
        foreach (double x in new[] { 0.0, 0.005, 0.5, 1.0, 298.15, 123.456 })
        {
            Assert.Close(x, RealUnits.TimeToFemtoseconds(RealUnits.FemtosecondsToTime(x)), 1e-15, 1e-14, "time");
            Assert.Close(x, RealUnits.VelocityToAngstromPerFemtosecond(RealUnits.AngstromPerFemtosecondToVelocity(x)), 1e-15, 1e-14, "velocity");
            Assert.Close(x, RealUnits.EnergyToKelvin(RealUnits.KelvinToEnergy(x)), 1e-15, 1e-14, "temperature");
            Assert.Close(x, RealUnits.BarToPressure(RealUnits.PressureToBar(x)), 1e-15, 1e-14, "pressure");
            Assert.Close(x, RealUnits.GramsPerCubicCentimetreToDensity(RealUnits.DensityToGramsPerCubicCentimetre(x)), 1e-15, 1e-14, "density");
        }
    }

    /// <summary>Argon's reduced units, expressed through RealUnits, must equal the ArgonUnits values.</summary>
    [Test]
    public static void ArgonReducedUnits_AgreeWithArgonUnits()
    {
        double sigma = ArgonUnits.SigmaAngstrom;
        double mass = ArgonUnits.MassAtomicUnits;
        double eps = ArgonEpsilonKj;

        Assert.Close(ArgonUnits.EnergyToKilojoulesPerMole(1.0), eps, 0, CrossCheckRelTol, "ε in kJ/mol");
        double tau = sigma * Math.Sqrt(mass / eps);
        Assert.Close(ArgonUnits.TauFemtoseconds, RealUnits.TimeToFemtoseconds(tau), 0, CrossCheckRelTol, "τ in fs");
        double sigma3 = sigma * sigma * sigma;
        Assert.Close(ArgonUnits.PressureToBar(1.0), RealUnits.PressureToBar(eps / sigma3), 0, CrossCheckRelTol, "ε/σ³ in bar");
        Assert.Close(ArgonUnits.DensityToGramsPerCubicCentimetre(1.0),
            RealUnits.DensityToGramsPerCubicCentimetre(mass / sigma3), 0, CrossCheckRelTol, "m/σ³ in g/cm³");
        Assert.Close(ArgonUnits.TemperatureToKelvin(0.722), RealUnits.EnergyToKelvin(0.722 * eps), 0, CrossCheckRelTol, "T* = 0.722 in K");
    }

    /// <summary>
    /// The same argon liquid run in reduced units and in real units (positions scaled by σ, LJ in
    /// Å and kJ/mol, masses in g/mol, timestep converted through fs) must follow the same path,
    /// with no change to the integrator.
    /// </summary>
    [Test]
    public static void ArgonDynamics_SameInReducedAndRealUnits()
    {
        const int steps = 200;
        const double dtReduced = 0.005;
        double sigma = ArgonUnits.SigmaAngstrom;
        double mass = ArgonUnits.MassAtomicUnits;
        double eps = ArgonEpsilonKj;
        double tauCore = RealUnits.FemtosecondsToTime(ArgonUnits.TauFemtoseconds);

        var reduced = Helpers.Fcc(2, 1, 1, jitter: 0.05);
        VelocityInitializer.Assign(reduced, 0.722, seed: 3);

        var realPositions = new Vec3[reduced.Count];
        for (int i = 0; i < reduced.Count; i++) realPositions[i] = sigma * reduced.UnwrappedPosition(i);
        var realMasses = new double[reduced.Count];
        Array.Fill(realMasses, mass);
        var real = new ParticleSystem(new Box(sigma * reduced.Box.Lengths), realPositions, realMasses, reduced.Elements);
        for (int i = 0; i < reduced.Count; i++) real.Velocities[i] = (sigma / tauCore) * reduced.Velocities[i];

        var simReduced = new Simulation(reduced, new VelocityVerlet(new ForceCalculator(new LennardJones()), dtReduced));
        var simReal = new Simulation(real, new VelocityVerlet(
            new ForceCalculator(new LennardJones(sigma, eps, 2.5 * sigma)), dtReduced * tauCore));

        Assert.Close(simReduced.State.PotentialEnergy * eps, simReal.State.PotentialEnergy, 0, CrossCheckRelTol, "initial PE");
        for (int s = 0; s < steps; s++)
        {
            simReduced.Step();
            simReal.Step();
        }
        for (int i = 0; i < reduced.Count; i++)
        {
            Vec3 d = real.UnwrappedPosition(i) / sigma - reduced.UnwrappedPosition(i);
            Assert.True(d.Length < 1e-7, $"atom {i} position differs by {d.Length:G3} σ after {steps} steps");
        }
        Assert.Close(simReduced.State.TotalEnergy * eps, simReal.State.TotalEnergy, 0, 1e-7, "total energy");
    }
}
