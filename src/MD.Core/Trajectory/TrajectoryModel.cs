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
    public const string ReducedLjUnits = "lj-reduced";
}

/// <summary>Real-unit scale factors for converting reduced quantities (written in the header).</summary>
public sealed record UnitScales(double SigmaAngstrom, double EpsilonKelvin, double MassAmu, double TauFemtoseconds);

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

/// <summary>A bond between two atoms (zero-based indices). Always empty in the LJ-only phase.</summary>
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
