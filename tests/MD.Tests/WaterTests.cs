using System;
using System.Linq;
using MD.Core;

namespace MD.Tests;

public static class WaterTests
{
    /// <summary>Parameters as recorded from the sources (SpcFw doc comment; WATER_PLAN.md Section 7), converted to kJ/mol.</summary>
    [Test]
    public static void SpcFwParameters_AreTheRecordedValues()
    {
        Assert.Close(1059.162 * 4.184, SpcFw.Bond.ForceConstant, 0, 1e-15, "k_b in kJ/(mol Å²)");
        Assert.Close(1.012, SpcFw.Bond.Length, 0, 0, "r0 in Å");
        Assert.Close(75.90 * 4.184, SpcFw.Angle.ForceConstant, 0, 1e-15, "k_a in kJ/(mol rad²)");
        Assert.Close(113.24, SpcFw.Angle.Angle * 180.0 / Math.PI, 0, 1e-14, "θ0 in degrees");
        Assert.Close(-0.82, SpcFw.Oxygen.Charge, 0, 0, "q_O");
        Assert.Close(0.41, SpcFw.Hydrogen.Charge, 0, 0, "q_H");
        Assert.Close(0.0, SpcFw.Molecule.TotalCharge(), 1e-15, 0, "molecule net charge");
        Assert.Close(0.1554 * 4.184, SpcFw.Oxygen.LjEpsilon, 0, 1e-15, "ε_OO in kJ/mol");
        Assert.Close(3.1655, SpcFw.Oxygen.LjSigma, 0, 0, "σ_OO in Å");
        Assert.True(!SpcFw.Hydrogen.HasLj, "H has no LJ site");
        Assert.Close(18.0154, SpcFw.MolarMass, 0, 1e-12, "molar mass");
    }

    /// <summary>
    /// A single molecule relaxed under the full force field (bonded + LJ + Ewald) in a 30 Å box must
    /// reach the model's equilibrium geometry: all its own non-bonded pairs are excluded, and the
    /// image forces (≈ 0.03 kJ/(mol·Å)) move it by less than 1e-4 Å and 0.01°. If the exclusions
    /// leaked, the intramolecular Coulomb term (hundreds of kJ/mol) would distort it grossly.
    /// </summary>
    [Test]
    public static void SingleMolecule_RelaxesToModelGeometry()
    {
        var top = new Topology();
        top.AddMolecule(SpcFw.Molecule);
        var box = new Box(new Vec3(30, 30, 30));
        var ff = SpcFw.ForceField(top, ljCutoff: 9.0, ewaldRealCutoff: 9.0, ewaldAccuracy: 1e-8);
        var o = new Vec3(15, 15, 15);
        var x = new[] { o, o + new Vec3(1.05, 0.0, 0.0), o + 0.98 * new Vec3(Math.Cos(1.85), Math.Sin(1.85), 0.0) };
        var f = new Vec3[3];
        // Steepest descent. Stability needs step × (largest curvature) < 2; O sits in two bonds, so the
        // largest curvature is about 3k_b ≈ 13300 and the step must be below 1.5e-4.
        const double step = 1e-4;
        for (int it = 0; it < 1500; it++)
        {
            ff.Compute(x, box, f);
            for (int i = 0; i < 3; i++) x[i] += step * f[i];
        }
        ff.Compute(x, box, f);
        double fmax = f.Max(v => v.Length);
        Assert.True(fmax < 0.1, $"residual force {fmax:G3} kJ/(mol Å)");
        double r1 = (x[1] - x[0]).Length, r2 = (x[2] - x[0]).Length;
        double theta = Math.Acos(Vec3.Dot(x[1] - x[0], x[2] - x[0]) / (r1 * r2)) * 180.0 / Math.PI;
        Assert.Close(SpcFw.BondLengthAngstrom, r1, 1e-4, 0, "O–H1 length");
        Assert.Close(SpcFw.BondLengthAngstrom, r2, 1e-4, 0, "O–H2 length");
        Assert.Close(SpcFw.AngleDegrees, theta, 0.01, 0, "H–O–H angle");
    }

    [Test]
    public static void WaterBox_HasTargetDensityAndIntactMolecules()
    {
        var (box, x, top) = WaterBox.Cubic(3, 1.0, seed: 1);
        Assert.Equal(27, top.MoleculeCount, "molecules");
        Assert.Equal(81, top.AtomCount, "atoms");
        Assert.Close(0.0, top.TotalCharge(), 1e-13, 0, "net charge");
        double density = RealUnits.DensityToGramsPerCubicCentimetre(top.Masses().Sum() / box.Volume);
        Assert.Close(1.0, density, 0, 1e-12, "mass density");
        for (int m = 0; m < 27; m++)
        {
            var (r1, r2, theta) = Geometry(x, box, m);
            Assert.Close(SpcFw.BondLengthAngstrom, r1, 1e-12, 0, "O–H1");
            Assert.Close(SpcFw.BondLengthAngstrom, r2, 1e-12, 0, "O–H2");
            Assert.Close(SpcFw.AngleDegrees, theta, 1e-9, 0, "H–O–H");
        }
    }

    /// <summary>
    /// 27 SPC/Fw molecules at 1 g/cm³, full force field, NVE for 1 ps at Δt = 0.5 fs from a 300 K start.
    /// The grid start releases potential energy, so the temperature rises (no thermostat until W5).
    /// Criteria:
    /// - the temperature stays below 1000 K;
    /// - every O–H length and H–O–H angle stays within 6 thermal standard deviations of r0 and θ0 at
    ///   1000 K (σ = √(RT/k): 0.043 Å and 9.3°), i.e. a molecule is never stretched or bent beyond what
    ///   heat alone could do;
    /// - each H stays closer to its own O than to any other (no dissociation or H transfer);
    /// - max |E − E0| &lt; 1 % of KE0; drift (mean of the last tenth minus the first tenth) &lt; 0.1 % of KE0.
    /// </summary>
    [Test]
    public static void Water_StaysBondedAndConservesEnergy_1ps()
    {
        var (sim, box, s) = WaterSimulation(dtFs: 0.5, seed: 1);
        double e0 = sim.Initial.TotalEnergy, ke0 = sim.Initial.KineticEnergy;
        const int steps = 2000;
        double maxDev = 0, first = 0, last = 0;
        for (int n = 1; n <= steps; n++)
        {
            sim.Step();
            double e = sim.State.TotalEnergy;
            maxDev = Math.Max(maxDev, Math.Abs(e - e0));
            if (n <= steps / 10) first += e;
            if (n > steps - steps / 10) last += e;
            if (n % 20 == 0)
            {
                CheckIntact(s, box, n);
                double t = RealUnits.EnergyToKelvin(sim.State.Temperature);
                Assert.True(t < 1000.0, $"step {n}: temperature {t:F0} K");
            }
        }
        double drift = (last - first) / (steps / 10);
        Assert.True(maxDev < 0.01 * ke0, $"max |E − E0| = {maxDev:G3} kJ/mol vs KE0 {ke0:G4}");
        Assert.True(Math.Abs(drift) < 1e-3 * ke0, $"drift {drift:G3} kJ/mol vs KE0 {ke0:G4}");
    }

    /// <summary>With the full SPC/Fw force field the integration error is still second order: halving Δt (0.5 → 0.25 fs) cuts max |E − E0| over 200 fs by 3.5 to 4.5.</summary>
    [Test]
    public static void Water_EnergyErrorIsSecondOrder()
    {
        double MaxDeviation(double dtFs)
        {
            var (sim, _, _) = WaterSimulation(dtFs, seed: 1);
            double e0 = sim.Initial.TotalEnergy, maxDev = 0;
            int steps = (int)Math.Round(200.0 / dtFs);
            for (int n = 0; n < steps; n++)
            {
                sim.Step();
                maxDev = Math.Max(maxDev, Math.Abs(sim.State.TotalEnergy - e0));
            }
            return maxDev;
        }
        double ratio = MaxDeviation(0.5) / MaxDeviation(0.25);
        Assert.True(ratio > 3.5 && ratio < 4.5, $"ratio {ratio:F2}, expected ≈ 4");
    }

    // ---- helpers --------------------------------------------------------------------------------

    private static (Simulation Sim, Box Box, ParticleSystem System) WaterSimulation(double dtFs, int seed)
    {
        var (box, x, top) = WaterBox.Cubic(3, 1.0, seed);
        var s = new ParticleSystem(box, x, top.Masses(), top.Elements());
        VelocityInitializer.Assign(s, RealUnits.KelvinToEnergy(300.0), seed);
        // Box side 9.32 Å: cutoffs just under half the box. Test settings, not production ones.
        double rc = 0.49 * box.Lengths.X;
        var ff = SpcFw.ForceField(top, ljCutoff: rc, ewaldRealCutoff: rc, ewaldAccuracy: 1e-8);
        return (new Simulation(s, new VelocityVerlet(ff, RealUnits.FemtosecondsToTime(dtFs))), box, s);
    }

    private static (double R1, double R2, double Theta) Geometry(Vec3[] x, Box box, int m)
    {
        Vec3 a = box.MinimumImage(x[3 * m + 1] - x[3 * m]);
        Vec3 b = box.MinimumImage(x[3 * m + 2] - x[3 * m]);
        double theta = Math.Acos(Vec3.Dot(a, b) / (a.Length * b.Length)) * 180.0 / Math.PI;
        return (a.Length, b.Length, theta);
    }

    private static void CheckIntact(ParticleSystem s, Box box, int step)
    {
        double rt = RealUnits.KelvinToEnergy(1000.0);
        double sigmaR = Math.Sqrt(rt / SpcFw.Bond.ForceConstant);
        double sigmaTheta = Math.Sqrt(rt / SpcFw.Angle.ForceConstant) * 180.0 / Math.PI;
        int molecules = s.Count / 3;
        for (int m = 0; m < molecules; m++)
        {
            var (r1, r2, theta) = Geometry(s.Positions, box, m);
            double r0 = SpcFw.BondLengthAngstrom, t0 = SpcFw.AngleDegrees;
            Assert.True(Math.Abs(r1 - r0) < 6 * sigmaR && Math.Abs(r2 - r0) < 6 * sigmaR, $"step {step}, molecule {m}: O–H {r1:F3}, {r2:F3} Å");
            Assert.True(Math.Abs(theta - t0) < 6 * sigmaTheta, $"step {step}, molecule {m}: angle {theta:F1}°");
            for (int h = 1; h <= 2; h++)
            {
                double own = box.MinimumImage(s.Positions[3 * m + h] - s.Positions[3 * m]).Length;
                for (int other = 0; other < molecules; other++)
                {
                    if (other == m) continue;
                    double d = box.MinimumImage(s.Positions[3 * m + h] - s.Positions[3 * other]).Length;
                    Assert.True(d > own, $"step {step}: H {3 * m + h} is closer to O of molecule {other} ({d:F3} Å) than to its own ({own:F3} Å)");
                }
            }
        }
    }
}
