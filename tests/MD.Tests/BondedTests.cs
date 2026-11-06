using System;
using System.Linq;
using MD.Core;

namespace MD.Tests;

public static class BondedTests
{
    // Placeholder parameters of realistic size for an O–H stretch and H–O–H bend. NOT SPC/Fw values;
    // those come from the paper in W4. Real units: kJ/mol, Å, rad, g/mol.
    private const double Kb = 4000.0;
    private const double R0 = 1.0;
    private const double Ka = 300.0;
    private static readonly double Theta0 = 109.47 * Math.PI / 180.0;
    private const double MassO = 16.0;
    private const double MassH = 1.0;

    private static readonly AtomType Ox = new("OX", "O", MassO, charge: -0.8, ljSigma: 3.0, ljEpsilon: 0.6);
    private static readonly AtomType Hx = new("HX", "H", MassH, charge: 0.4);
    private static readonly HarmonicBond Bond = new(Kb, R0);
    private static readonly HarmonicAngle Bend = new(Ka, Theta0);
    private static readonly MoleculeTemplate Flexible = new(
        "X3F", new[] { Ox, Hx, Hx },
        new[] { new BondTerm(0, 1, Bond), new BondTerm(0, 2, Bond) },
        new[] { new AngleTerm(1, 0, 2, Bend) });

    /// <summary>Equilibrium geometry with O at <paramref name="o"/>, in the xy plane.</summary>
    private static Vec3[] Equilibrium(Vec3 o, double r = R0, double theta = double.NaN)
    {
        if (double.IsNaN(theta)) theta = Theta0;
        double h = 0.5 * theta;
        return new[] { o, o + r * new Vec3(Math.Sin(h), Math.Cos(h), 0), o + r * new Vec3(-Math.Sin(h), Math.Cos(h), 0) };
    }

    private static Topology Molecules(int n)
    {
        var top = new Topology();
        for (int m = 0; m < n; m++) top.AddMolecule(Flexible);
        return top;
    }

    [Test]
    public static void Forces_AreNegativeEnergyGradient()
    {
        var rng = new Random(17);
        var box = new Box(new Vec3(6.0, 7.0, 8.0));
        var top = Molecules(4);
        var positions = new Vec3[top.AtomCount];
        for (int m = 0; m < 4; m++)
        {
            // Distorted geometries, placed near the box edges so some molecules straddle the boundary.
            var o = new Vec3(rng.NextDouble() < 0.5 ? 0.2 : 5.8, rng.NextDouble() * 7.0, rng.NextDouble() * 8.0);
            var p = Equilibrium(o, R0 + 0.1 * (rng.NextDouble() - 0.5), Theta0 + 0.4 * (rng.NextDouble() - 0.5));
            for (int k = 0; k < 3; k++)
            {
                var jitter = 0.05 * new Vec3(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5, rng.NextDouble() - 0.5);
                positions[3 * m + k] = box.Wrap(p[k] + jitter, Int3.Zero).Position;
            }
        }
        var bonded = new BondedForces(top);
        var forces = new Vec3[positions.Length];
        double e0 = bonded.Compute(positions, box, forces).PotentialEnergy;
        Assert.True(e0 > 1.0, $"test geometry should be away from equilibrium (E = {e0:G3} kJ/mol)");
        var scratch = new Vec3[positions.Length];
        const double h = 1e-6;
        for (int i = 0; i < positions.Length; i++)
        {
            for (int c = 0; c < 3; c++)
            {
                var p = (Vec3[])positions.Clone();
                p[i] = Shift(positions[i], c, h);
                double ePlus = bonded.Compute(p, box, scratch).PotentialEnergy;
                p[i] = Shift(positions[i], c, -h);
                double eMinus = bonded.Compute(p, box, scratch).PotentialEnergy;
                Assert.Close(-(ePlus - eMinus) / (2 * h), Component(forces[i], c), 1e-5, 1e-7, $"force on atom {i}, component {c}");
            }
        }
        Vec3 sum = Vec3.Zero;
        foreach (var f in forces) sum += f;
        Assert.True(sum.Length < 1e-10, $"bonded forces sum to {sum.Length:G3}");
    }

    [Test]
    public static void Equilibrium_HasZeroEnergyAndForce()
    {
        var box = new Box(new Vec3(20, 20, 20));
        var top = Molecules(1);
        var forces = new Vec3[3];
        var r = new BondedForces(top).Compute(Equilibrium(new Vec3(10, 10, 10)), box, forces);
        Assert.Close(0.0, r.PotentialEnergy, 1e-20, 0, "energy at equilibrium");
        foreach (var f in forces) Assert.True(f.Length < 1e-10, $"force at equilibrium {f.Length:G3}");
    }

    [Test]
    public static void MissingParameters_AreRejected()
    {
        var bare = new MoleculeTemplate("X3", new[] { Ox, Hx, Hx },
            new[] { new BondTerm(0, 1), new BondTerm(0, 2, Bond) }, new[] { new AngleTerm(1, 0, 2, Bend) });
        var top = new Topology();
        top.AddMolecule(bare);
        Assert.Throws<ArgumentException>(() => new BondedForces(top), "bond without parameters");
    }

    /// <summary>
    /// Harmonic frequencies of a bent XY₂ molecule (bond k_b, angle k_θ, no cross terms), from the
    /// Wilson GF method with μ = 1/m:
    ///   antisymmetric  λ3 = k_b(μ_Y + μ_X(1 − cos θ))
    ///   symmetric      λ1 + λ2 = k_b(μ_Y + μ_X(1 + cos θ)) + (2k_θ/r²)(μ_Y + μ_X(1 − cos θ))
    ///                  λ1·λ2   = (2 k_b k_θ / r²)·μ_Y(μ_Y + 2μ_X)
    /// (these agree with Herzberg's valence-force formulas for XY₂, with 2α = θ and k_δ = k_θ).
    /// The engine's mass-weighted Hessian (finite differences of its forces) must have these three
    /// eigenvalues, plus six zero modes (translation and rotation).
    /// </summary>
    [Test]
    public static void NormalModes_MatchWilsonGF()
    {
        var (l1, l2, l3) = PredictedEigenvalues();
        var box = new Box(new Vec3(20, 20, 20));
        var top = Molecules(1);
        var bonded = new BondedForces(top);
        var x0 = Equilibrium(new Vec3(10, 10, 10));
        double[] masses = top.Masses();

        var hessian = new double[9, 9];
        var fPlus = new Vec3[3];
        var fMinus = new Vec3[3];
        const double h = 1e-5;
        for (int j = 0; j < 9; j++)
        {
            var p = (Vec3[])x0.Clone();
            p[j / 3] = Shift(x0[j / 3], j % 3, h);
            bonded.Compute(p, box, fPlus);
            p[j / 3] = Shift(x0[j / 3], j % 3, -h);
            bonded.Compute(p, box, fMinus);
            for (int i = 0; i < 9; i++)
            {
                double dF = Component(fPlus[i / 3], i % 3) - Component(fMinus[i / 3], i % 3);
                hessian[i, j] = -dF / (2 * h) / Math.Sqrt(masses[i / 3] * masses[j / 3]);
            }
        }
        for (int i = 0; i < 9; i++)
            for (int j = 0; j < i; j++)
                hessian[i, j] = hessian[j, i] = 0.5 * (hessian[i, j] + hessian[j, i]);

        double[] eig = SymmetricEigenvalues(hessian).OrderBy(v => v).ToArray();
        for (int k = 0; k < 6; k++) Assert.True(Math.Abs(eig[k]) < 1e-4 * l3, $"zero mode {k}: {eig[k]:G3}");
        var expected = new[] { l1, l2, l3 }.OrderBy(v => v).ToArray();
        for (int k = 0; k < 3; k++) Assert.Close(expected[k], eig[6 + k], 0, 1e-6, $"vibrational eigenvalue {k}");
    }

    /// <summary>
    /// Dynamics check: displace the two bonds antisymmetrically by ±1e-3 Å and run NVE at a small
    /// timestep (0.05 fs). The bond-length difference must oscillate with period 2π/√λ3. Expected
    /// shifts: Verlet (ωΔt)²/24 ≈ 5e-5, anharmonic ≈ (1e-3)² — both below the 2e-4 tolerance.
    /// </summary>
    [Test]
    public static void AntisymmetricStretch_HasPredictedPeriod()
    {
        var (_, _, l3) = PredictedEigenvalues();
        double expectedPeriod = 2 * Math.PI / Math.Sqrt(l3);
        var box = new Box(new Vec3(20, 20, 20));
        var top = Molecules(1);
        var x = Equilibrium(new Vec3(10, 10, 10));
        const double eps = 1e-3;
        x[1] += eps * (x[1] - x[0]) / R0;
        x[2] -= eps * (x[2] - x[0]) / R0;
        var s = new ParticleSystem(box, x, top.Masses(), top.Elements());
        double dt = RealUnits.FemtosecondsToTime(0.05);
        var sim = new Simulation(s, new VelocityVerlet(new BondedForces(top), dt));

        var crossings = new System.Collections.Generic.List<double>();
        double prev = Asymmetry(s);
        int steps = (int)(25 * expectedPeriod / dt);
        for (int n = 1; n <= steps; n++)
        {
            sim.Step();
            double cur = Asymmetry(s);
            if (prev < 0 && cur >= 0) crossings.Add((n - 1 + prev / (prev - cur)) * dt); // upward zero crossing
            prev = cur;
        }
        Assert.True(crossings.Count >= 20, $"expected ≥ 20 periods, found {crossings.Count} crossings");
        double measured = (crossings[^1] - crossings[0]) / (crossings.Count - 1);
        Assert.Close(expectedPeriod, measured, 0, 2e-4,
            $"antisymmetric stretch period (expected {RealUnits.TimeToFemtoseconds(expectedPeriod):F4} fs)");
    }

    /// <summary>
    /// NVE energy conservation for 8 flexible molecules at about 300 K (bonded forces only; no
    /// intermolecular forces yet). Stated criteria (PLAN.md Section 8):
    /// - second order: the largest |E − E0| over 5 ps shrinks by 4 (3.5 to 4.5) when Δt halves;
    /// - at 0.5 fs, the largest |E − E0| is below 2 % of the starting kinetic energy. Velocity Verlet
    ///   on an oscillator of frequency ω fluctuates by about (ωΔt)²/8 of its energy, ≈ 1.4 % here;
    /// - no drift: the mean of E over the last tenth differs from the first tenth by &lt; 0.1 % of KE0.
    /// </summary>
    [Test]
    public static void Energy_IsConservedAtHalfFemtosecond()
    {
        var at05 = RunNve(0.5, 5000.0);
        var at025 = RunNve(0.25, 5000.0);
        double ratio = at05.MaxDeviation / at025.MaxDeviation;
        Assert.True(ratio > 3.5 && ratio < 4.5, $"fluctuation ratio for Δt 0.5 → 0.25 fs is {ratio:F2}, expected ≈ 4");
        Assert.True(at05.MaxDeviation < 0.02 * at05.Ke0, $"max |E − E0| = {at05.MaxDeviation:G3} kJ/mol vs KE0 {at05.Ke0:G3}");
        Assert.True(Math.Abs(at05.Drift) < 1e-3 * at05.Ke0, $"drift {at05.Drift:G3} kJ/mol vs KE0 {at05.Ke0:G3}");
    }

    /// <summary>Run forward, reverse velocities, run back: the start must be recovered (PLAN.md Section 8).</summary>
    [Test]
    public static void Dynamics_AreTimeReversible()
    {
        var (s, sim) = MoleculeGas(molecules: 8, dtFs: 0.5, seed: 4);
        var start = Enumerable.Range(0, s.Count).Select(s.UnwrappedPosition).ToArray();
        for (int n = 0; n < 2000; n++) sim.Step();
        sim.ReverseVelocities();
        for (int n = 0; n < 2000; n++) sim.Step();
        for (int i = 0; i < s.Count; i++)
        {
            double d = (s.UnwrappedPosition(i) - start[i]).Length;
            Assert.True(d < 1e-8, $"atom {i} returns within {d:G3} Å");
        }
    }

    /// <summary>An isolated molecule with random velocities (including rotation) keeps its total momentum and angular momentum.</summary>
    [Test]
    public static void IsolatedMolecule_ConservesMomentumAndAngularMomentum()
    {
        var box = new Box(new Vec3(30, 30, 30));
        var top = Molecules(1);
        var x = Equilibrium(new Vec3(15, 15, 15), R0 + 0.03, Theta0 - 0.1);
        var s = new ParticleSystem(box, x, top.Masses(), top.Elements());
        var rng = new Random(9);
        for (int i = 0; i < 3; i++) s.Velocities[i] = 2.0 * new Vec3(rng.NextDouble() - 0.5, rng.NextDouble() - 0.5, rng.NextDouble() - 0.5);
        var sim = new Simulation(s, new VelocityVerlet(new BondedForces(top), RealUnits.FemtosecondsToTime(0.5)));
        Vec3 p0 = s.TotalMomentum();
        Vec3 l0 = AngularMomentum(s);
        Assert.True(l0.Length > 0.1, "test needs a rotating molecule");
        for (int n = 0; n < 20000; n++) sim.Step();
        Assert.True((s.TotalMomentum() - p0).Length < 1e-12, $"momentum change {(s.TotalMomentum() - p0).Length:G3}");
        Vec3 dl = AngularMomentum(s) - l0;
        Assert.True(dl.Length < 1e-11 * l0.Length, $"angular momentum change {dl.Length:G3} of {l0.Length:G3}");
    }

    // ---- helpers --------------------------------------------------------------------------------

    private static (double L1, double L2, double L3) PredictedEigenvalues()
    {
        double muX = 1.0 / MassO, muY = 1.0 / MassH, c = Math.Cos(Theta0), r2 = R0 * R0;
        double l3 = Kb * (muY + muX * (1 - c));
        double sum = Kb * (muY + muX * (1 + c)) + 2 * Ka / r2 * (muY + muX * (1 - c));
        double prod = 2 * Kb * Ka / r2 * muY * (muY + 2 * muX);
        double disc = Math.Sqrt(sum * sum - 4 * prod);
        return (0.5 * (sum - disc), 0.5 * (sum + disc), l3);
    }

    private static double Asymmetry(ParticleSystem s)
    {
        Vec3 o = s.UnwrappedPosition(0);
        return (s.UnwrappedPosition(1) - o).Length - (s.UnwrappedPosition(2) - o).Length;
    }

    private static Vec3 AngularMomentum(ParticleSystem s)
    {
        Vec3 com = Vec3.Zero;
        for (int i = 0; i < s.Count; i++) com += s.Masses[i] * s.UnwrappedPosition(i);
        com /= s.TotalMass();
        Vec3 l = Vec3.Zero;
        for (int i = 0; i < s.Count; i++) l += s.Masses[i] * Vec3.Cross(s.UnwrappedPosition(i) - com, s.Velocities[i]);
        return l;
    }

    /// <summary>Molecules on a coarse grid (well apart), slightly distorted, with 300 K velocities.</summary>
    private static (ParticleSystem System, Simulation Sim) MoleculeGas(int molecules, double dtFs, int seed)
    {
        var rng = new Random(seed);
        var box = new Box(new Vec3(20, 20, 20));
        var top = Molecules(molecules);
        var x = new Vec3[top.AtomCount];
        for (int m = 0; m < molecules; m++)
        {
            var o = new Vec3(5 + 10 * (m % 2), 5 + 10 * (m / 2 % 2), 5 + 10 * (m / 4 % 2));
            var p = Equilibrium(o, R0 + 0.02 * (rng.NextDouble() - 0.5), Theta0 + 0.1 * (rng.NextDouble() - 0.5));
            for (int k = 0; k < 3; k++) x[3 * m + k] = p[k];
        }
        var s = new ParticleSystem(box, x, top.Masses(), top.Elements());
        VelocityInitializer.Assign(s, RealUnits.KelvinToEnergy(300.0), seed);
        var sim = new Simulation(s, new VelocityVerlet(new BondedForces(top), RealUnits.FemtosecondsToTime(dtFs)));
        return (s, sim);
    }

    private static (double MaxDeviation, double Drift, double Ke0) RunNve(double dtFs, double durationFs)
    {
        var (_, sim) = MoleculeGas(molecules: 8, dtFs, seed: 2);
        double e0 = sim.Initial.TotalEnergy;
        double ke0 = sim.Initial.KineticEnergy;
        int steps = (int)Math.Round(durationFs / dtFs);
        int tenth = steps / 10;
        double maxDev = 0, firstSum = 0, lastSum = 0;
        for (int n = 1; n <= steps; n++)
        {
            sim.Step();
            double e = sim.State.TotalEnergy;
            maxDev = Math.Max(maxDev, Math.Abs(e - e0));
            if (n <= tenth) firstSum += e;
            if (n > steps - tenth) lastSum += e;
        }
        return (maxDev, (lastSum - firstSum) / tenth, ke0);
    }

    /// <summary>Eigenvalues of a small symmetric matrix by cyclic Jacobi rotations.</summary>
    private static double[] SymmetricEigenvalues(double[,] matrix)
    {
        int n = matrix.GetLength(0);
        var a = (double[,])matrix.Clone();
        for (int sweep = 0; sweep < 100; sweep++)
        {
            double off = 0;
            for (int p = 0; p < n; p++)
                for (int q = p + 1; q < n; q++)
                    off += a[p, q] * a[p, q];
            if (off < 1e-30) break;
            for (int p = 0; p < n; p++)
            {
                for (int q = p + 1; q < n; q++)
                {
                    if (a[p, q] == 0) continue;
                    double theta = 0.5 * (a[q, q] - a[p, p]) / a[p, q];
                    double t = Math.Sign(theta == 0 ? 1 : theta) / (Math.Abs(theta) + Math.Sqrt(theta * theta + 1));
                    double c = 1 / Math.Sqrt(t * t + 1), s = t * c;
                    for (int k = 0; k < n; k++)
                    {
                        double akp = a[k, p], akq = a[k, q];
                        a[k, p] = c * akp - s * akq;
                        a[k, q] = s * akp + c * akq;
                    }
                    for (int k = 0; k < n; k++)
                    {
                        double apk = a[p, k], aqk = a[q, k];
                        a[p, k] = c * apk - s * aqk;
                        a[q, k] = s * apk + c * aqk;
                    }
                }
            }
        }
        return Enumerable.Range(0, n).Select(i => a[i, i]).ToArray();
    }

    private static Vec3 Shift(Vec3 v, int c, double h) =>
        c == 0 ? new Vec3(v.X + h, v.Y, v.Z) : c == 1 ? new Vec3(v.X, v.Y + h, v.Z) : new Vec3(v.X, v.Y, v.Z + h);

    private static double Component(Vec3 v, int c) => c == 0 ? v.X : c == 1 ? v.Y : v.Z;
}
