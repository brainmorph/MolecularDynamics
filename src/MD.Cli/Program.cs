using System;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text.Json;
using MD.Core;
using MD.Core.Trajectory;

namespace MD.Cli;

public static class Program
{
    private const string Usage =
        "Usage:\n" +
        "  md run <settings.json>          Run an NVE simulation (argon, or \"system\": \"water\") and write a trajectory and log\n" +
        "  md dump <file.traj> [--all]     Print a trajectory's header and frames (first/last 3 unless --all)\n";

    public static int Main(string[] args)
    {
        CultureInfo.DefaultThreadCurrentCulture = CultureInfo.InvariantCulture;
        CultureInfo.CurrentCulture = CultureInfo.InvariantCulture;
        try
        {
            return args switch
            {
                ["run", var path] => Run(path),
                ["dump", var path] => Dump(path, all: false),
                ["dump", var path, "--all"] => Dump(path, all: true),
                _ => Fail(Usage),
            };
        }
        catch (Exception e) when (e is SettingsException or JsonException or IOException or TrajectoryFormatException or UnauthorizedAccessException)
        {
            return Fail($"error: {e.Message}");
        }
    }

    private static int Fail(string message)
    {
        Console.Error.WriteLine(message);
        return 1;
    }

    /// <summary>Everything the shared run loop needs, built per system (argon or water).</summary>
    private sealed record RunPlan(
        Simulation Sim,
        TrajectoryHeader Header,
        Bond[] Bonds,
        string[] Preamble,
        string ColumnHeader,
        Func<ThermoState, ThermoState, string> FormatState,
        string EnergyUnit,
        long Steps,
        int TrajectoryEvery,
        int LogEvery,
        string Output,
        string? Log);

    private static int Run(string settingsPath)
    {
        string json = File.ReadAllText(settingsPath);
        RunPlan plan = RunSettings.SystemOf(json) switch
        {
            RunSettings.SystemName => ArgonPlan(RunSettings.Parse(json)),
            WaterRunSettings.SystemName => WaterPlan(WaterRunSettings.Parse(json)),
            var other => throw new SettingsException($"Unknown system \"{other}\" (expected \"{RunSettings.SystemName}\" or \"{WaterRunSettings.SystemName}\")."),
        };

        string baseDir = Path.GetDirectoryName(Path.GetFullPath(settingsPath))!;
        string outputPath = Path.GetFullPath(Path.Combine(baseDir, plan.Output));
        string logPath = plan.Log is null
            ? Path.ChangeExtension(outputPath, ".log")
            : Path.GetFullPath(Path.Combine(baseDir, plan.Log));
        Directory.CreateDirectory(Path.GetDirectoryName(outputPath)!);
        Directory.CreateDirectory(Path.GetDirectoryName(logPath)!);

        Simulation sim = plan.Sim;
        ThermoState start = sim.Initial;

        using var log = new StreamWriter(logPath);
        void Log(string line)
        {
            Console.WriteLine(line);
            log.WriteLine(line);
        }

        foreach (string line in plan.Preamble) Log(line);
        Log($"# trajectory {outputPath}");
        Log(plan.ColumnHeader);

        using var traj = new TrajectoryWriter(File.Create(outputPath));
        traj.WriteHeader(plan.Header);
        traj.WriteFrame(TrajectoryWriter.FrameFrom(sim, plan.Bonds));
        Log(plan.FormatState(start, start));

        double maxAbsDrift = 0;
        var clock = System.Diagnostics.Stopwatch.StartNew();
        for (long k = 1; k <= plan.Steps; k++)
        {
            sim.Step();
            ThermoState t = sim.State;
            maxAbsDrift = Math.Max(maxAbsDrift, Math.Abs(t.TotalEnergy - start.TotalEnergy));
            if (!double.IsFinite(t.TotalEnergy))
            {
                Log($"# ERROR: energy is not finite at step {k}; stopping. Try a smaller timestep.");
                return 2;
            }
            if (k % plan.TrajectoryEvery == 0) traj.WriteFrame(TrajectoryWriter.FrameFrom(sim, plan.Bonds));
            if (k % plan.LogEvery == 0 || k == plan.Steps) Log(plan.FormatState(t, start));
        }
        double seconds = clock.Elapsed.TotalSeconds;

        ThermoState end = sim.State;
        Log($"# done: {plan.Steps} steps in {seconds:F2} s ({(seconds > 0 ? plan.Steps / seconds : 0):F0} steps/s)");
        Log($"# energy ({plan.EnergyUnit}): E0 = {start.TotalEnergy:F8}, final drift = {end.TotalEnergy - start.TotalEnergy:E3}, " +
            $"max |E - E0| = {maxAbsDrift:E3} ({(start.TotalEnergy != 0 ? 100 * maxAbsDrift / Math.Abs(start.TotalEnergy) : 0):E2} % of |E0|, " +
            $"{(start.KineticEnergy != 0 ? 100 * maxAbsDrift / start.KineticEnergy : 0):E2} % of KE0)");
        Log($"# |total momentum| = {sim.System.TotalMomentum().Length:E2}");
        return 0;
    }

    private static RunPlan ArgonPlan(RunSettings settings)
    {
        var (box, positions) = settings.Box is null
            ? Lattice.Fcc(settings.Cells[0], settings.Cells[1], settings.Cells[2], settings.Density)
            : Lattice.Chain(settings.Cells[0], settings.Cells[1], settings.Cells[2], settings.Density,
                settings.ClusterCount, settings.ClusterGap, new Vec3(settings.Box[0], settings.Box[1], settings.Box[2]));
        int n = positions.Length;
        var masses = Enumerable.Repeat(1.0, n).ToArray();
        var elements = Enumerable.Range(0, n).Select(i => settings.Elements[i % settings.Elements.Length]).ToArray();
        var system = new ParticleSystem(box, positions, masses, elements);
        VelocityInitializer.Assign(system, settings.Temperature, settings.Seed);
        var potential = new LennardJones(cutoff: settings.Cutoff);
        var sim = new Simulation(system, new VelocityVerlet(new ForceCalculator(potential), settings.Timestep));
        ThermoState start = sim.Initial;

        var preamble = new System.Collections.Generic.List<string>
        {
            "# LJ molecular dynamics, NVE, velocity Verlet. Reduced LJ units (sigma, epsilon, m, k_B = 1).",
        };
        double actualDensity = n / box.Volume;
        string densityNote = settings.Box is null ? "" : $" (lattice spacing from density {settings.Density})";
        string clusterNote = settings.ClusterCount > 1 ? $" x {settings.ClusterCount} clusters, gap {settings.ClusterGap}" : "";
        preamble.Add($"# atoms {n}  cells {string.Join("x", settings.Cells)}{clusterNote}  density {actualDensity:G6}{densityNote}  " +
            $"box {box.Lengths.X:F6} x {box.Lengths.Y:F6} x {box.Lengths.Z:F6}");
        preamble.Add($"# timestep {settings.Timestep}  steps {settings.Steps}  cutoff {settings.Cutoff} (energy-shifted)  seed {settings.Seed}");
        preamble.Add($"# argon scale: dt = {ArgonUnits.TimeToFemtoseconds(settings.Timestep):F3} fs, " +
            $"T0 = {ArgonUnits.TemperatureToKelvin(settings.Temperature):F2} K, " +
            $"density = {ArgonUnits.DensityToGramsPerCubicCentimetre(actualDensity):F4} g/cm3");
        if (box.MinLength < 2 * settings.Cutoff)
        {
            preamble.Add("# note: box is smaller than 2 x cutoff; forces sum over all periodic images within the cutoff.");
        }
        if (n < 100)
        {
            preamble.Add($"# note: {n} atoms is a correctness check only; too few for meaningful thermodynamic averages.");
        }

        var header = new TrajectoryHeader(
            TrajectoryFormat.Version,
            TrajectoryFormat.ReducedLjUnits,
            new UnitScales(ArgonUnits.SigmaAngstrom, ArgonUnits.EpsilonOverBoltzmannKelvin, ArgonUnits.MassAtomicUnits, ArgonUnits.TauFemtoseconds),
            elements,
            masses,
            settings.Timestep,
            start.KineticEnergy,
            start.PotentialEnergy,
            start.TotalEnergy,
            settings.ToJson());

        return new RunPlan(sim, header, Array.Empty<Bond>(), preamble.ToArray(),
            "#     step         time           KE             PE              E          drift(E-E0)        T             P",
            (t, s0) => string.Create(CultureInfo.InvariantCulture,
                $"{t.Step,10} {t.Time,12:F4} {t.KineticEnergy,14:F8} {t.PotentialEnergy,14:F8} {t.TotalEnergy,14:F8} {t.TotalEnergy - s0.TotalEnergy,14:E4} {t.Temperature,12:F6} {t.Pressure,12:F6}"),
            "epsilon", settings.Steps, settings.TrajectoryEvery, settings.LogEvery, settings.Output, settings.Log);
    }

    private static RunPlan WaterPlan(WaterRunSettings settings)
    {
        var (box, positions, topology) = WaterBox.Cubic(settings.MoleculesPerSide, settings.DensityGramsPerCm3, settings.Seed);
        var system = new ParticleSystem(box, positions, topology.Masses(), topology.Elements());
        VelocityInitializer.Assign(system, RealUnits.KelvinToEnergy(settings.TemperatureKelvin), settings.Seed);
        var forceField = SpcFw.ForceField(topology, settings.LjCutoffAngstrom, settings.EwaldCutoffAngstrom, settings.EwaldAccuracy);
        double dt = RealUnits.FemtosecondsToTime(settings.TimestepFemtoseconds);
        var sim = new Simulation(system, new VelocityVerlet(forceField, dt));
        ThermoState start = sim.Initial;
        var ewald = Ewald.ForAccuracy(settings.EwaldCutoffAngstrom, settings.EwaldAccuracy, RealUnits.CoulombConstant);

        var preamble = new System.Collections.Generic.List<string>
        {
            "# SPC/Fw flexible water, NVE, velocity Verlet. Real units: A, kJ/mol, g/mol, time unit 100 fs; T from R.",
            $"# molecules {topology.MoleculeCount} ({settings.MoleculesPerSide}^3 grid, random orientations)  atoms {topology.AtomCount}  " +
            $"density {RealUnits.DensityToGramsPerCubicCentimetre(topology.Masses().Sum() / box.Volume):F4} g/cm3  box {box.Lengths.X:F4} A",
            $"# timestep {settings.TimestepFemtoseconds} fs  steps {settings.Steps} ({settings.Steps * settings.TimestepFemtoseconds / 1000.0:G4} ps)  " +
            $"T0 {settings.TemperatureKelvin} K  seed {settings.Seed}",
            $"# LJ cutoff {settings.LjCutoffAngstrom} A (energy-shifted)  Ewald: real cutoff {settings.EwaldCutoffAngstrom} A, " +
            $"accuracy {settings.EwaldAccuracy:G3}, alpha {ewald.Alpha:F4} 1/A, k cutoff {ewald.ReciprocalCutoff:F4} 1/A",
            "# note: no thermostat yet (WATER_PLAN W5). The grid start releases energy and the water heats up.",
        };
        if (box.MinLength < 2 * Math.Max(settings.LjCutoffAngstrom, settings.EwaldCutoffAngstrom))
        {
            preamble.Add("# note: box is smaller than 2 x cutoff; pair terms sum over all periodic images within the cutoff.");
        }

        var header = new TrajectoryHeader(
            TrajectoryFormat.Version,
            TrajectoryFormat.MolecularRealUnits,
            UnitScales.Real,
            topology.Elements(),
            topology.Masses(),
            dt,
            start.KineticEnergy,
            start.PotentialEnergy,
            start.TotalEnergy,
            settings.ToJson());

        return new RunPlan(sim, header, TrajectoryWriter.BondsFrom(topology), preamble.ToArray(),
            "#     step     time(ps)      KE(kJ/mol)     PE(kJ/mol)      E(kJ/mol)     drift(E-E0)       T(K)",
            (t, s0) => string.Create(CultureInfo.InvariantCulture,
                $"{t.Step,10} {RealUnits.TimeToFemtoseconds(t.Time) / 1000.0,12:F4} {t.KineticEnergy,14:F4} {t.PotentialEnergy,14:F4} {t.TotalEnergy,14:F4} {t.TotalEnergy - s0.TotalEnergy,14:E4} {RealUnits.EnergyToKelvin(t.Temperature),10:F2}"),
            "kJ/mol", settings.Steps, settings.TrajectoryEvery, settings.LogEvery, settings.Output, settings.Log);
    }

    private static int Dump(string path, bool all)
    {
        using var reader = new TrajectoryReader(File.OpenRead(path));
        TrajectoryHeader h = reader.Header;
        Console.WriteLine($"format version : {h.Version}");
        Console.WriteLine($"units          : {h.Units} (1 length = {h.Scales.SigmaAngstrom} A, 1 energy = {h.Scales.EpsilonKelvin} K x k_B, 1 mass = {h.Scales.MassAmu} u, 1 time = {h.Scales.TauFemtoseconds:F3} fs)");
        Console.WriteLine($"atoms          : {h.AtomCount}");
        Console.WriteLine($"elements       : {string.Join(" ", h.Elements)}");
        Console.WriteLine($"masses         : {string.Join(" ", h.Masses)}");
        Console.WriteLine($"timestep       : {h.Timestep}");
        Console.WriteLine($"initial KE/PE/E: {h.InitialKineticEnergy:F8} {h.InitialPotentialEnergy:F8} {h.InitialTotalEnergy:F8}");
        Console.WriteLine("settings       :");
        Console.WriteLine(h.SettingsJson);

        int count = 0;
        var tail = new System.Collections.Generic.Queue<TrajectoryFrame>();
        foreach (TrajectoryFrame f in reader.ReadFrames())
        {
            if (all || count < 3) PrintFrame(f, h);
            else
            {
                tail.Enqueue(f);
                if (tail.Count > 3) tail.Dequeue();
            }
            count++;
        }
        if (!all && count > 6) Console.WriteLine($"\n... {count - 6} frames omitted ...");
        foreach (var f in tail) PrintFrame(f, h);
        Console.WriteLine($"\nframes: {count}");
        return 0;
    }

    private static void PrintFrame(TrajectoryFrame f, TrajectoryHeader h)
    {
        Console.WriteLine();
        Console.WriteLine($"frame step {f.Step}  time {f.Time:F4}  KE {f.KineticEnergy:F8}  PE {f.PotentialEnergy:F8}  E {f.TotalEnergy:F8}  (E - E0 = {f.TotalEnergy - h.InitialTotalEnergy:E3})");
        Console.WriteLine($"  box {f.BoxLengths.X:F6} {f.BoxLengths.Y:F6} {f.BoxLengths.Z:F6}  bonds {f.Bonds.Length}");
        for (int i = 0; i < f.Positions.Length; i++)
        {
            Vec3 p = f.Positions[i];
            Int3 im = f.Images[i];
            Console.WriteLine($"  {i,5} {h.Elements[i],-3} {p.X,12:F6} {p.Y,12:F6} {p.Z,12:F6}   image {im.X,3} {im.Y,3} {im.Z,3}");
        }
        foreach (Bond b in f.Bonds) Console.WriteLine($"  bond {b.I} {b.J}");
    }
}
