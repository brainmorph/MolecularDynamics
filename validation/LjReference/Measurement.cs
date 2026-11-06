using System;
using System.Linq;
using MD.Core;

// Shared measurement protocol for the Lennard-Jones validation runs.
//
// Each run: FCC start, melt at a high temperature, cool to the target, hold, then NVE production.
// Velocity rescaling is used only during preparation. Measured quantities come from NVE dynamics,
// and energy drift is measured over production only.
public static class Measurement
{
    public const int CellsPerAxis = 5;     // 5 × 5 × 5 FCC cells = 500 atoms
    public const double Timestep = 0.004;  // τ, as in Johnson et al.
    public const double HotTemperature = 2.0;
    public const int MeltSteps = 10000;
    public const int CoolSteps = 10000;
    public const int HoldSteps = 10000;
    public const int Blocks = 10;
    public const int StepsPerBlock = 10000;
    public const int SampleEvery = 10;

    public sealed record Result(
        double Density,
        double Temperature,
        int Atoms,
        double MeanTemperature,
        double Pressure,
        double PressureError,
        double PotentialPerAtom,
        double PotentialPerAtomError,
        double Drift);

    // Long-range corrections for a truncated, unshifted-force LJ fluid (g(r) = 1 beyond r_c).
    public static double PressureTail(double rho, double rc) =>
        16.0 * Math.PI / 3.0 * rho * rho * (2.0 / 3.0 * Math.Pow(rc, -9) - Math.Pow(rc, -3));

    public static double EnergyTail(double rho, double rc) =>
        8.0 * Math.PI / 3.0 * rho * (1.0 / 3.0 * Math.Pow(rc, -9) - Math.Pow(rc, -3));

    // Energy per atom added back for the shift: (ρ/2)·(4πr_c³/3)·U(r_c), with U(r_c) < 0.
    public static double ShiftEnergy(double rho, double rc)
    {
        double uRc = 4.0 * (Math.Pow(rc, -12) - Math.Pow(rc, -6));
        return 0.5 * rho * (4.0 * Math.PI * rc * rc * rc / 3.0) * uRc;
    }

    // With holdTemperature, production also rescales velocities to the target every SampleEvery
    // steps. That makes this an NVT average at exactly the target temperature, and the energy drift
    // figure then means nothing. Without it, production is NVE and the measured temperature can
    // drift from the target.
    public static Result Run(double density, double temperature, double cutoff, int seed = 1, bool holdTemperature = false)
    {
        var (box, positions) = Lattice.Fcc(CellsPerAxis, CellsPerAxis, CellsPerAxis, density);
        int n = positions.Length;
        var system = new ParticleSystem(
            box,
            positions,
            Enumerable.Repeat(1.0, n).ToArray(),
            Enumerable.Repeat("Ar", n).ToArray());
        VelocityInitializer.Assign(system, temperature, seed);
        var sim = new Simulation(system, new VelocityVerlet(new ForceCalculator(new LennardJones(cutoff: cutoff)), Timestep));

        Prepare(sim, system, temperature);
        ThermoState productionStart = sim.State;

        double[] blockP = new double[Blocks];
        double[] blockU = new double[Blocks];
        double sumT = 0, sumP = 0, sumU = 0;
        long samples = 0;
        for (int b = 0; b < Blocks; b++)
        {
            double bp = 0, bu = 0;
            int bs = 0;
            for (int k = 1; k <= StepsPerBlock; k++)
            {
                sim.Step();
                if (holdTemperature && k % SampleEvery == 0) RescaleToTemperature(system, temperature);
                if (k % SampleEvery != 0) continue;
                ThermoState s = sim.State;
                bp += s.Pressure;
                bu += s.PotentialEnergy / n;
                sumT += s.Temperature;
                bs++;
            }
            blockP[b] = bp / bs;
            blockU[b] = bu / bs;
            sumP += bp;
            sumU += bu;
            samples += bs;
        }

        double drift = (sim.State.TotalEnergy - productionStart.TotalEnergy) / Math.Abs(productionStart.TotalEnergy);
        return new Result(
            density,
            temperature,
            n,
            sumT / samples,
            sumP / samples,
            StandardError(blockP),
            sumU / samples,
            StandardError(blockU),
            drift);
    }

    private static void Prepare(Simulation sim, ParticleSystem s, double target)
    {
        for (int step = 1; step <= MeltSteps; step++)
        {
            sim.Step();
            if (step % 10 == 0) RescaleToTemperature(s, HotTemperature);
        }
        for (int step = 1; step <= CoolSteps; step++)
        {
            sim.Step();
            if (step % 10 == 0)
            {
                double t = HotTemperature + (target - HotTemperature) * step / CoolSteps;
                RescaleToTemperature(s, t);
            }
        }
        for (int step = 1; step <= HoldSteps; step++)
        {
            sim.Step();
            if (step % 10 == 0) RescaleToTemperature(s, target);
        }
    }

    private static void RescaleToTemperature(ParticleSystem s, double target)
    {
        double current = Thermo.Temperature(Thermo.KineticEnergy(s), s.Count);
        if (current <= 0) return;
        double scale = Math.Sqrt(target / current);
        for (int i = 0; i < s.Count; i++) s.Velocities[i] *= scale;
    }

    private static double StandardError(double[] blockMeans)
    {
        double mean = blockMeans.Average();
        double var = blockMeans.Sum(x => (x - mean) * (x - mean)) / (blockMeans.Length - 1);
        return Math.Sqrt(var / blockMeans.Length);
    }
}
