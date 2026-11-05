using MD.Core;

namespace MD.Tests;

public static class UnitsTests
{
    /// <summary>Published argon values: τ ≈ 2.156 ps, ε/σ³ ≈ 419 bar, m/σ³ ≈ 1.680 g/cm³.</summary>
    [Test]
    public static void DerivedUnits_MatchKnownArgonValues()
    {
        Assert.Close(2156.0, ArgonUnits.TauFemtoseconds, 1.0, 0, "τ in fs");
        Assert.Close(419.0, ArgonUnits.PressureToBar(1.0), 1.0, 0, "ε/σ³ in bar");
        Assert.Close(1.680, ArgonUnits.DensityToGramsPerCubicCentimetre(1.0), 0.002, 0, "m/σ³ in g/cm³");
        Assert.Close(0.9961, ArgonUnits.EnergyToKilojoulesPerMole(1.0), 0.0005, 0, "ε in kJ/mol");
        Assert.Close(86.49, ArgonUnits.TemperatureToKelvin(0.722), 0.01, 0, "T* = 0.722 in K");
        Assert.Close(3.405, ArgonUnits.LengthToAngstrom(1.0), 0, 0, "σ in Å");
        // Default timestep 0.005 τ ≈ 10.8 fs.
        Assert.Close(10.78, ArgonUnits.TimeToFemtoseconds(0.005), 0.01, 0, "0.005 τ in fs");
    }

    [Test]
    public static void Conversions_RoundTrip()
    {
        foreach (double x in new[] { 0.0, 0.005, 0.8442, 1.0, 2.5, 123.456 })
        {
            Assert.Close(x, ArgonUnits.AngstromToLength(ArgonUnits.LengthToAngstrom(x)), 1e-15, 1e-14, "length");
            Assert.Close(x, ArgonUnits.FemtosecondsToTime(ArgonUnits.TimeToFemtoseconds(x)), 1e-15, 1e-14, "time");
            Assert.Close(x, ArgonUnits.KelvinToTemperature(ArgonUnits.TemperatureToKelvin(x)), 1e-15, 1e-14, "temperature");
            Assert.Close(x, ArgonUnits.KilojoulesPerMoleToEnergy(ArgonUnits.EnergyToKilojoulesPerMole(x)), 1e-15, 1e-14, "energy");
            Assert.Close(x, ArgonUnits.BarToPressure(ArgonUnits.PressureToBar(x)), 1e-15, 1e-14, "pressure");
            Assert.Close(x, ArgonUnits.GramsPerCubicCentimetreToDensity(ArgonUnits.DensityToGramsPerCubicCentimetre(x)), 1e-15, 1e-14, "density");
        }
    }
}
