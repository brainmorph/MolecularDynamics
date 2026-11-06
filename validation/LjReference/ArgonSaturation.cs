using System;
using System.Threading.Tasks;
using MD.Core;

// Compares the LJ fluid with real argon at saturation, using NIST WebBook data for argon
// (saturation at 86.485 K: liquid density 35.056 mol/l, vapour pressure 0.092862 MPa).
//
// At saturation the liquid and vapour pressures are equal, so the LJ liquid density at which the
// simulated pressure equals P_sat is the model's saturated liquid density. That density is
// compared with argon's real saturated liquid density at the same reduced temperature.
//
// Argon parameters (ArgonUnits): σ = 3.405 Å, ε/k_B = 119.8 K, m = 39.948 u. Real argon is not
// exactly LJ, so disagreement here is partly model error, not only simulation error.
public static class ArgonSaturation
{
    private const double SaturationKelvin = 86.485;
    private const double LiquidMolPerLitre = 35.056;   // NIST WebBook, saturated liquid
    private const double VaporPressureMPa = 0.092862;   // NIST WebBook, saturated vapour pressure
    private const double MolarMassGramsPerMole = 39.948;
    private const double Cutoff = 4.0;

    private static readonly double[] Densities = { 0.82, 0.83, 0.835, 0.84, 0.8442, 0.85 };

    public static bool Run()
    {
        double tStar = SaturationKelvin / ArgonUnits.EpsilonOverBoltzmannKelvin;
        double rhoArgonGramsPerCm3 = LiquidMolPerLitre * MolarMassGramsPerMole / 1000.0;
        double rhoArgonStar = ArgonUnits.GramsPerCubicCentimetreToDensity(rhoArgonGramsPerCm3);
        double pSatMegaPascal = VaporPressureMPa;
        double pSatStar = ArgonUnits.BarToPressure(pSatMegaPascal * 10.0);

        Console.WriteLine($"Argon at {SaturationKelvin} K (T* = {tStar:F4}):");
        Console.WriteLine($"  saturated liquid density  {rhoArgonGramsPerCm3:F4} g/cm³  = {rhoArgonStar:F4} σ⁻³");
        Console.WriteLine($"  vapour pressure           {pSatMegaPascal:F6} MPa   = {pSatStar:F6} ε/σ³");
        Console.WriteLine($"Simulating {Densities.Length} densities at T* = {tStar:F4}, r_c = {Cutoff}, in parallel");
        Console.WriteLine("(temperature held at the target during production, so this is an NVT average)");
        Console.WriteLine();

        var results = new Measurement.Result[Densities.Length];
        Parallel.For(0, Densities.Length, i =>
        {
            results[i] = Measurement.Run(Densities[i], tStar, Cutoff, holdTemperature: true);
        });

        // Full-LJ pressure (tail added, as in the paper's P* column). Pressure is increasing in density
        // on the liquid branch, so find where it crosses the saturation pressure.
        Console.WriteLine("  rho*     P (raw, cut)     P (full LJ)    P (full, MPa)   measured T");
        double[] pFull = new double[results.Length];
        for (int i = 0; i < results.Length; i++)
        {
            var r = results[i];
            pFull[i] = r.Pressure + Measurement.PressureTail(r.Density, Cutoff);
            double pMPa = ArgonUnits.PressureToBar(pFull[i]) / 10.0;
            Console.WriteLine($"  {r.Density,-7} {r.Pressure,9:F4} ± {r.PressureError:F4}   {pFull[i],9:F4}       {pMPa,9:F2}     {r.MeanTemperature,8:F4}");
        }
        Console.WriteLine();

        bool found = TryFindCrossing(Densities, pFull, pSatStar, out double rhoLj);
        if (!found)
        {
            Console.WriteLine("No density in the scanned range brackets the saturation pressure.");
            return false;
        }

        double pctDiff = 100.0 * (rhoLj - rhoArgonStar) / rhoArgonStar;
        Console.WriteLine($"LJ density at P = P_sat      {rhoLj:F4} σ⁻³  = {rhoLj * ArgonUnits.DensityToGramsPerCubicCentimetre(1.0):F4} g/cm³ (using argon's σ, ε, m)");
        Console.WriteLine($"Argon saturated liquid       {rhoArgonStar:F4} σ⁻³  = {rhoArgonGramsPerCm3:F4} g/cm³");
        Console.WriteLine($"Difference                   {pctDiff:+0.00;-0.00;0.00}% (model vs argon)");
        Console.WriteLine();
        Console.WriteLine("Note: this is a comparison, not a pass/fail test. The LJ model is not exact for argon.");
        return true;
    }

    // Linear interpolation of density where the (increasing) pressure curve crosses pSat.
    private static bool TryFindCrossing(double[] rho, double[] p, double pSat, out double rhoAtSat)
    {
        for (int i = 0; i + 1 < rho.Length; i++)
        {
            if ((p[i] - pSat) * (p[i + 1] - pSat) <= 0 && p[i + 1] != p[i])
            {
                double f = (pSat - p[i]) / (p[i + 1] - p[i]);
                rhoAtSat = rho[i] + f * (rho[i + 1] - rho[i]);
                return true;
            }
        }
        rhoAtSat = 0;
        return false;
    }
}
