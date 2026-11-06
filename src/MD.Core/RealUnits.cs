using System;

namespace MD.Core;

/// <summary>
/// Real-unit system for molecular (water) simulations, decision D1 in WATER_PLAN.md.
///
/// The core uses a consistent set, so the integrator needs no conversion factors:
///   length  Å
///   energy  kJ/mol
///   mass    g/mol (u)
///   charge  e
///   time    Å·√((g/mol)/(kJ/mol)) = 1e-10 m · √(1e-3 kg / 1e3 J) = 1e-13 s = 100 fs (exact)
/// Force is kJ/(mol·Å). Temperature enters only through R = k_B·N_A (kJ/mol/K).
///
/// Times are given and reported in fs; these functions convert at input/output only, in the same
/// way <see cref="ArgonUnits"/> does for the reduced LJ system. The two systems are never mixed in
/// one simulation.
/// </summary>
public static class RealUnits
{
    // CODATA 2018 / SI values. The first four are exact by the 2019 SI definition.
    public const double BoltzmannJoulePerKelvin = 1.380649e-23;
    public const double AvogadroPerMole = 6.02214076e23;
    public const double ElementaryChargeCoulomb = 1.602176634e-19;
    public const double VacuumPermittivityFaradPerMetre = 8.8541878128e-12;

    /// <summary>Core time unit in seconds, derived from the length, energy and mass units.</summary>
    public static readonly double TimeUnitSeconds = 1e-10 * Math.Sqrt(1e-3 / 1e3);

    /// <summary>Core time unit in fs (100).</summary>
    public static readonly double TimeUnitFemtoseconds = TimeUnitSeconds * 1e15;

    /// <summary>Molar gas constant R = k_B·N_A in kJ/(mol·K).</summary>
    public static readonly double GasConstantKilojoulesPerMoleKelvin = BoltzmannJoulePerKelvin * AvogadroPerMole / 1000.0;

    /// <summary>
    /// Coulomb constant 1/(4πε₀) in kJ·Å/(mol·e²): the energy of two unit charges 1 Å apart.
    /// </summary>
    public static readonly double CoulombConstant =
        ElementaryChargeCoulomb * ElementaryChargeCoulomb / (4.0 * Math.PI * VacuumPermittivityFaradPerMetre)
        * AvogadroPerMole / 1000.0 / 1e-10;

    public static double FemtosecondsToTime(double fs) => fs / TimeUnitFemtoseconds;
    public static double TimeToFemtoseconds(double time) => time * TimeUnitFemtoseconds;

    /// <summary>Velocity in Å/fs to core units (Å per 100 fs).</summary>
    public static double AngstromPerFemtosecondToVelocity(double v) => v * TimeUnitFemtoseconds;
    public static double VelocityToAngstromPerFemtosecond(double v) => v / TimeUnitFemtoseconds;

    /// <summary>k_B·T in kJ/mol for a temperature in K.</summary>
    public static double KelvinToEnergy(double kelvin) => kelvin * GasConstantKilojoulesPerMoleKelvin;
    public static double EnergyToKelvin(double kjPerMol) => kjPerMol / GasConstantKilojoulesPerMoleKelvin;

    /// <summary>Pressure in kJ/(mol·Å³) to bar.</summary>
    public static double PressureToBar(double p) => p * 1e3 / AvogadroPerMole / 1e-30 / 1e5;
    public static double BarToPressure(double bar) => bar * 1e5 * 1e-30 * AvogadroPerMole / 1e3;

    /// <summary>Mass density in (g/mol)/Å³ to g/cm³.</summary>
    public static double DensityToGramsPerCubicCentimetre(double d) => d / AvogadroPerMole / 1e-24;
    public static double GramsPerCubicCentimetreToDensity(double gPerCm3) => gPerCm3 * AvogadroPerMole * 1e-24;
}
