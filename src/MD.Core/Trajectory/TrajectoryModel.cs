using System;
using System.Collections.Generic;

namespace MD.Core.Trajectory;

/// <summary>Constants for trajectory format version 1. See docs/trajectory-format.md.</summary>
public static class TrajectoryFormat
{
    /// <summary>File magic: ASCII "LJMDTRAJ".</summary>
    public static ReadOnlySpan<byte> Magic => "LJMDTRAJ"u8;

    /// <summary>Frame marker: ASCII "FRAM".</summary>
    public static ReadOnlySpan<byte> FrameMarker => "FRAM"u8;

    public const uint Version = 1;

    /// <summary>Argon runs: σ, ε, m, τ, k_B = 1.</summary>
    public const string ReducedLjUnits = "lj-reduced";

    /// <summary>Molecular runs (<see cref="MD.Core.RealUnits"/>): Å, kJ/mol, g/mol, time unit 100 fs.</summary>
    public const string MolecularRealUnits = "real";
}

/// <summary>
/// What one unit of length, energy, mass and time in the file equals (written in the header). For
/// "lj-reduced" these are argon's σ, ε/k_B, m and τ; for "real" they are given by <see cref="Real"/>.
/// </summary>
public sealed record UnitScales(double SigmaAngstrom, double EpsilonKelvin, double MassAmu, double TauFemtoseconds)
{
    /// <summary>Scales for "real" units: 1 Å, 1 kJ/mol = 1/R K, 1 u, 100 fs.</summary>
    public static UnitScales Real { get; } = new(
        1.0, 1.0 / MD.Core.RealUnits.GasConstantKilojoulesPerMoleKelvin, 1.0, MD.Core.RealUnits.TimeUnitFemtoseconds);
}

public sealed record TrajectoryHeader(
    uint Version,
    string Units,
    UnitScales Scales,
    IReadOnlyList<string> Elements,
    IReadOnlyList<double> Masses,
    double Timestep,
    double InitialKineticEnergy,
    double InitialPotentialEnergy,
    double InitialTotalEnergy,
    string SettingsJson)
{
    public int AtomCount => Elements.Count;
}

/// <summary>A bond between two atoms (zero-based indices). Empty for argon runs; the topology's bonds for molecular runs.</summary>
public readonly record struct Bond(uint I, uint J);

public sealed record TrajectoryFrame(
    long Step,
    double Time,
    double KineticEnergy,
    double PotentialEnergy,
    double TotalEnergy,
    Vec3 BoxLengths,
    Vec3[] Positions,
    Int3[] Images,
    Bond[] Bonds);

public sealed class TrajectoryFormatException : Exception
{
    public TrajectoryFormatException(string message) : base(message) { }
    public TrajectoryFormatException(string message, Exception inner) : base(message, inner) { }
}
