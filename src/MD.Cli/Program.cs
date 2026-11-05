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
        "  md run <settings.json>          Run an NVE simulation and write a trajectory and log\n" +
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

    private static int Run(string settingsPath)
    {
        RunSettings settings = RunSettings.Parse(File.ReadAllText(settingsPath));
        string baseDir = Path.GetDirectoryName(Path.GetFullPath(settingsPath))!;
        string outputPath = Path.GetFullPath(Path.Combine(baseDir, settings.Output));
        string logPath = settings.Log is null
            ? Path.ChangeExtension(outputPath, ".log")
            : Path.GetFullPath(Path.Combine(baseDir, settings.Log));
        Directory.CreateDirectory(Path.GetDirectoryName(outputPath)!);
        Directory.CreateDirectory(Path.GetDirectoryName(logPath)!);

        // Build the system.
        var (box, positions) = Lattice.Fcc(settings.Cells[0], settings.Cells[1], settings.Cells[2], settings.Density);
        int n = positions.Length;
        var masses = Enumerable.Repeat(1.0, n).ToArray();
        var elements = Enumerable.Range(0, n).Select(i => settings.Elements[i % settings.Elements.Length]).ToArray();
        var system = new ParticleSystem(box, positions, masses, elements);
        VelocityInitializer.Assign(system, settings.Temperature, settings.Seed);
        var potential = new LennardJones(cutoff: settings.Cutoff);
        var sim = new Simulation(system, new VelocityVerlet(new ForceCalculator(potential), settings.Timestep));
        ThermoState start = sim.Initial;

        using var log = new StreamWriter(logPath);
        void Log(string line)
        {
            Console.WriteLine(line);
            log.WriteLine(line);
        }

        Log("# LJ molecular dynamics, NVE, velocity Verlet. Reduced LJ units (sigma, epsilon, m, k_B = 1).");
        Log($"# atoms {n}  cells {string.Join("x", settings.Cells)}  density {settings.Density}  " +
            $"box {box.Lengths.X:F6} x {box.Lengths.Y:F6} x {box.Lengths.Z:F6}");
        Log($"# timestep {settings.Timestep}  steps {settings.Steps}  cutoff {settings.Cutoff} (energy-shifted)  seed {settings.Seed}");
        Log($"# argon scale: dt = {ArgonUnits.TimeToFemtoseconds(settings.Timestep):F3} fs, " +
            $"T0 = {ArgonUnits.TemperatureToKelvin(settings.Temperature):F2} K, " +
            $"density = {ArgonUnits.DensityToGramsPerCubicCentimetre(settings.Density):F4} g/cm3");
        if (box.MinLength < 2 * settings.Cutoff)
        {
            Log($"# note: box is smaller than 2 x cutoff; forces sum over all periodic images within the cutoff.");
        }
        if (n < 100)
        {
            Log($"# note: {n} atoms is a correctness check only; too few for meaningful thermodynamic averages.");
        }
        Log($"# trajectory {outputPath}");
        Log("#     step         time           KE             PE              E          drift(E-E0)        T             P");

        void LogState(ThermoState t) => Log(string.Create(CultureInfo.InvariantCulture,
            $"{t.Step,10} {t.Time,12:F4} {t.KineticEnergy,14:F8} {t.PotentialEnergy,14:F8} {t.TotalEnergy,14:F8} {t.TotalEnergy - start.TotalEnergy,14:E4} {t.Temperature,12:F6} {t.Pressure,12:F6}"));

        using var traj = new TrajectoryWriter(File.Create(outputPath));
        traj.WriteHeader(new TrajectoryHeader(
            TrajectoryFormat.Version,
            TrajectoryFormat.ReducedLjUnits,
            new UnitScales(ArgonUnits.SigmaAngstrom, ArgonUnits.EpsilonOverBoltzmannKelvin, ArgonUnits.MassAtomicUnits, ArgonUnits.TauFemtoseconds),
            elements,
            masses,
            settings.Timestep,
            start.KineticEnergy,
            start.PotentialEnergy,
            start.TotalEnergy,
            settings.ToJson()));
        traj.WriteFrame(TrajectoryWriter.FrameFrom(sim));
        LogState(start);

        double maxAbsDrift = 0;
        var clock = System.Diagnostics.Stopwatch.StartNew();
        for (long k = 1; k <= settings.Steps; k++)
        {
            sim.Step();
            ThermoState t = sim.State;
            maxAbsDrift = Math.Max(maxAbsDrift, Math.Abs(t.TotalEnergy - start.TotalEnergy));
            if (!double.IsFinite(t.TotalEnergy))
            {
                Log($"# ERROR: energy is not finite at step {k}; stopping. Try a smaller timestep.");
                return 2;
            }
            if (k % settings.TrajectoryEvery == 0) traj.WriteFrame(TrajectoryWriter.FrameFrom(sim));
            if (k % settings.LogEvery == 0 || k == settings.Steps) LogState(t);
        }
        double seconds = clock.Elapsed.TotalSeconds;

        ThermoState end = sim.State;
        Log($"# done: {settings.Steps} steps in {seconds:F2} s ({(seconds > 0 ? settings.Steps / seconds : 0):F0} steps/s)");
        Log($"# energy: E0 = {start.TotalEnergy:F8}, final drift = {end.TotalEnergy - start.TotalEnergy:E3}, " +
            $"max |E - E0| = {maxAbsDrift:E3} ({(start.TotalEnergy != 0 ? 100 * maxAbsDrift / Math.Abs(start.TotalEnergy) : 0):E2} % of |E0|)");
        Log($"# |total momentum| = {system.TotalMomentum().Length:E2}");
        return 0;
    }

    private static int Dump(string path, bool all)
    {
        using var reader = new TrajectoryReader(File.OpenRead(path));
        TrajectoryHeader h = reader.Header;
        Console.WriteLine($"format version : {h.Version}");
        Console.WriteLine($"units          : {h.Units} (sigma = {h.Scales.SigmaAngstrom} A, epsilon/k_B = {h.Scales.EpsilonKelvin} K, m = {h.Scales.MassAmu} u, tau = {h.Scales.TauFemtoseconds:F3} fs)");
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
