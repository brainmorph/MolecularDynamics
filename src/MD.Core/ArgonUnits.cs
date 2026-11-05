namespace MD.Core;

/// <summary>
/// Conversion between reduced LJ units and real units, using argon parameters
/// (Rahman 1964: σ = 3.405 Å, ε/k_B = 119.8 K, m = 39.948 u).
///
/// The core works only in reduced units (length σ, energy ε, mass m, k_B = 1). These functions are
/// for input and output only.
///
/// Derived units:
///   time τ = σ·√(m/ε) ≈ 2.156 ps
///   pressure ε/σ³ ≈ 419 bar
///   mass density m/σ³ ≈ 1.680 g/cm³
/// </summary>
public static class ArgonUnits
{
    public const double SigmaAngstrom = 3.405;
    public const double EpsilonOverBoltzmannKelvin = 119.8;
    public const double MassAtomicUnits = 39.948;

    // CODATA / SI exact values.
    public const double BoltzmannJoulePerKelvin = 1.380649e-23;
    public const double AvogadroPerMole = 6.02214076e23;
    public const double AtomicMassUnitKilogram = 1.66053906660e-27;

    public static double SigmaMetres => SigmaAngstrom * 1e-10;
    public static double EpsilonJoules => EpsilonOverBoltzmannKelvin * BoltzmannJoulePerKelvin;
    public static double MassKilograms => MassAtomicUnits * AtomicMassUnitKilogram;

    /// <summary>τ = σ·√(m/ε), in seconds.</summary>
    public static double TauSeconds => SigmaMetres * System.Math.Sqrt(MassKilograms / EpsilonJoules);
    public static double TauFemtoseconds => TauSeconds * 1e15;

    public static double LengthToAngstrom(double reduced) => reduced * SigmaAngstrom;
    public static double AngstromToLength(double angstrom) => angstrom / SigmaAngstrom;

    public static double TimeToFemtoseconds(double reduced) => reduced * TauFemtoseconds;
    public static double FemtosecondsToTime(double fs) => fs / TauFemtoseconds;

    public static double TemperatureToKelvin(double reduced) => reduced * EpsilonOverBoltzmannKelvin;
    public static double KelvinToTemperature(double kelvin) => kelvin / EpsilonOverBoltzmannKelvin;

    /// <summary>Energy per atom (or per pair) in kJ/mol.</summary>
    public static double EnergyToKilojoulesPerMole(double reduced) => reduced * EpsilonJoules * AvogadroPerMole / 1000.0;
    public static double KilojoulesPerMoleToEnergy(double kjPerMol) => kjPerMol * 1000.0 / (EpsilonJoules * AvogadroPerMole);

    public static double PressureToBar(double reduced) =>
        reduced * EpsilonJoules / (SigmaMetres * SigmaMetres * SigmaMetres) / 1e5;
    public static double BarToPressure(double bar) =>
        bar * 1e5 * (SigmaMetres * SigmaMetres * SigmaMetres) / EpsilonJoules;

    /// <summary>Number density in reduced units (ρσ³) to mass density in g/cm³.</summary>
    public static double DensityToGramsPerCubicCentimetre(double reduced) =>
        reduced * MassKilograms / (SigmaMetres * SigmaMetres * SigmaMetres) / 1000.0;
    public static double GramsPerCubicCentimetreToDensity(double gPerCm3) =>
        gPerCm3 * 1000.0 * (SigmaMetres * SigmaMetres * SigmaMetres) / MassKilograms;
}
