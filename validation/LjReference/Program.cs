using System;
using System.Globalization;
using System.Linq;

// Compares the core's Lennard-Jones NVE engine against published simulation data:
//   J. K. Johnson, J. A. Zollweg, K. E. Gubbins, "The Lennard-Jones equation of state revisited",
//   Mol. Phys. 78, 591 (1993), Table 2.
//
// Johnson et al. ran 864 atoms with a cut-and-shifted potential at r_c = 4σ. Their table gives
// two versions of each quantity:
//   P*, U*   full LJ potential (shift added back, standard tail corrections applied)
//   Pc*, Uc* cut-and-shifted potential, no corrections (what a simulation like this one measures)
// This program compares its raw pressure and energy with Pc* and Uc*. It also checks that the
// tail and shift formulas turn Pc*, Uc* into P*, U*, which confirms how the columns are read.
//
// Usage:
//   dotnet run -c Release --project validation/LjReference                  all Johnson points
//   dotnet run -c Release --project validation/LjReference -- <index>       one Johnson point (0..2)
//   dotnet run -c Release --project validation/LjReference -- argon         real-argon saturation check
public static class Program
{
    private const double Cutoff = 4.0;      // r_c in σ, as in Johnson et al.
    private const double RelativeTolerance = 0.01; // 1% of the reference value, in addition to 3σ error

    private sealed record Row(double Density, double Temperature, double Pc, double Uc, double PFull, double UFull);
    private sealed record Target(double Density, double Temperature);

    // Table 2 rows, read from the paper's text layer (the columns are as labelled in the paper).
    private static readonly Row[] Rows =
    {
        new(0.80, 0.70, -0.357, -5.606, -0.525, -5.815),
        new(0.80, 0.75, -0.088, -5.556, -0.256, -5.765),
        new(0.90, 0.70,  1.350, -6.193,  1.140, -6.429),
        new(0.90, 0.75,  1.715, -6.129,  1.503, -6.365),
    };

    // State points to simulate. Each density has tabulated rows at T = 0.70 and 0.75. The reference
    // is interpolated linearly to the temperature each run actually measures, because the NVE
    // temperature sits a few thousandths away from the target.
    private static readonly Target[] Targets =
    {
        new(0.80, 0.70),
        new(0.80, 0.75),
        new(0.90, 0.75),
    };

    public static int Main(string[] args)
    {
        CultureInfo.CurrentCulture = CultureInfo.InvariantCulture;

        if (args.Length == 1 && args[0] == "argon")
        {
            return ArgonSaturation.Run() ? 0 : 1;
        }

        int[] indices = args.Length == 0
            ? Enumerable.Range(0, Targets.Length).ToArray()
            : new[] { int.Parse(args[0], CultureInfo.InvariantCulture) };

        CheckTailIdentities();

        bool allPass = true;
        foreach (int i in indices)
        {
            allPass &= Run(Targets[i]);
        }
        Console.WriteLine(allPass ? "\nOVERALL: PASS" : "\nOVERALL: FAIL");
        return allPass ? 0 : 1;
    }

    private static void CheckTailIdentities()
    {
        Console.WriteLine($"Reading check (r_c = {Cutoff}): full = cut-shifted + corrections");
        Console.WriteLine("  rho    T      P*(table)  Pc+tail     U*(table)  Uc+shift+tail");
        foreach (Row p in Rows)
        {
            double pCalc = p.Pc + Measurement.PressureTail(p.Density, Cutoff);
            double uCalc = p.Uc + Measurement.ShiftEnergy(p.Density, Cutoff) + Measurement.EnergyTail(p.Density, Cutoff);
            Console.WriteLine($"  {p.Density,-5} {p.Temperature,-5} {p.PFull,9:F3}  {pCalc,9:F3}   {p.UFull,9:F3}  {uCalc,9:F3}");
        }
        Console.WriteLine();
    }

    // Linear interpolation in T between the two tabulated rows at this density. Returns (Pc, Uc).
    private static (double Pc, double Uc) Reference(double density, double temperature)
    {
        Row[] rows = Rows.Where(r => r.Density == density).OrderBy(r => r.Temperature).ToArray();
        Row lo = rows[0], hi = rows[^1];
        double f = (temperature - lo.Temperature) / (hi.Temperature - lo.Temperature);
        return (lo.Pc + f * (hi.Pc - lo.Pc), lo.Uc + f * (hi.Uc - lo.Uc));
    }

    private static bool Run(Target p)
    {
        Console.WriteLine($"rho = {p.Density}, T = {p.Temperature}: r_c = {Cutoff}");
        Measurement.Result r = Measurement.Run(p.Density, p.Temperature, Cutoff);
        int total = Measurement.Blocks * Measurement.StepsPerBlock;

        (double refP, double refU) = Reference(p.Density, r.MeanTemperature);
        bool passP = Within(r.Pressure, r.PressureError, refP);
        bool passU = Within(r.PotentialPerAtom, r.PotentialPerAtomError, refU);

        Console.WriteLine($"  {r.Atoms} atoms, measured T {r.MeanTemperature:F4} (target {p.Temperature})");
        Console.WriteLine($"  energy drift   {r.Drift:P3} over {total} NVE steps");
        Console.WriteLine($"  quantity       simulated        reference at T={r.MeanTemperature:F4} (Pc*/Uc*)   difference   pass");
        Console.WriteLine($"  P (cut)        {r.Pressure,8:F4} ± {r.PressureError:F4}   {refP,12:F4}        {r.Pressure - refP,+8:F4}   {Flag(passP)}");
        Console.WriteLine($"  U/N (cut)      {r.PotentialPerAtom,8:F4} ± {r.PotentialPerAtomError:F4}   {refU,12:F4}        {r.PotentialPerAtom - refU,+8:F4}   {Flag(passU)}");
        Console.WriteLine();
        return passP && passU;
    }

    private static bool Within(double simulated, double standardError, double reference)
    {
        double tolerance = 3.0 * standardError + RelativeTolerance * Math.Abs(reference);
        return Math.Abs(simulated - reference) <= tolerance;
    }

    private static string Flag(bool pass) => pass ? "PASS" : "FAIL";
}
