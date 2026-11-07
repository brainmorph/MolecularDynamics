using System;
using System.Collections.Generic;
using System.Text.Json;

namespace MD.Cli;

/// <summary>
/// Settings for an SPC/Fw water run (<c>"system": "water"</c>), read from JSON. Values are in real units,
/// with the unit in each name: Å, K, fs, g/cm³. Unknown keys are rejected, including the argon keys,
/// so an argon setting cannot silently do nothing in a water run.
/// </summary>
public sealed class WaterRunSettings
{
    public const string SystemName = "water";

    /// <summary>Must be "water".</summary>
    public string System { get; set; } = SystemName;

    /// <summary>Molecules along each side of the starting cubic grid. Molecule count = n³.</summary>
    public int MoleculesPerSide { get; set; } = 4;

    /// <summary>Mass density in g/cm³. Sets the box side.</summary>
    public double DensityGramsPerCm3 { get; set; } = 1.0;

    /// <summary>Starting temperature in K, used to draw the velocities.</summary>
    public double TemperatureKelvin { get; set; } = 300.0;

    public double TimestepFemtoseconds { get; set; } = 0.5;
    public long Steps { get; set; } = 2000;
    public int Seed { get; set; } = 12345;

    /// <summary>O–O Lennard-Jones cutoff in Å (truncated and energy-shifted).</summary>
    public double LjCutoffAngstrom { get; set; } = 6.0;

    /// <summary>Ewald real-space cutoff in Å.</summary>
    public double EwaldCutoffAngstrom { get; set; } = 6.0;

    /// <summary>Ewald accuracy δ: erfc(α r_c) = δ sets α, and the reciprocal cutoff follows (see Ewald.ForAccuracy).</summary>
    public double EwaldAccuracy { get; set; } = 1e-6;

    /// <summary>Write a trajectory frame every N steps (frame 0 is always written).</summary>
    public int TrajectoryEvery { get; set; } = 10;

    /// <summary>Write a log line every N steps.</summary>
    public int LogEvery { get; set; } = 100;

    /// <summary>Trajectory path, relative to the settings file.</summary>
    public string Output { get; set; } = "out/water.traj";

    /// <summary>Log path, relative to the settings file. Defaults to the output path with ".log".</summary>
    public string? Log { get; set; }

    public static WaterRunSettings Parse(string json)
    {
        WaterRunSettings s = JsonSerializer.Deserialize<WaterRunSettings>(json, RunSettings.JsonOptions)
            ?? throw new SettingsException("Settings file is empty.");
        s.Validate();
        return s;
    }

    public string ToJson() => JsonSerializer.Serialize(this, RunSettings.JsonOptions);

    public void Validate()
    {
        var errors = new List<string>();
        if (System != SystemName) errors.Add($"system must be \"{SystemName}\"");
        if (MoleculesPerSide < 1 || MoleculesPerSide > 40) errors.Add("moleculesPerSide must be from 1 to 40");
        if (!IsPositive(DensityGramsPerCm3)) errors.Add("densityGramsPerCm3 must be positive");
        if (!(TemperatureKelvin >= 0) || double.IsInfinity(TemperatureKelvin)) errors.Add("temperatureKelvin must be zero or positive");
        if (!IsPositive(TimestepFemtoseconds)) errors.Add("timestepFemtoseconds must be positive");
        if (Steps < 0) errors.Add("steps must be zero or positive");
        if (!IsPositive(LjCutoffAngstrom)) errors.Add("ljCutoffAngstrom must be positive");
        if (!IsPositive(EwaldCutoffAngstrom)) errors.Add("ewaldCutoffAngstrom must be positive");
        if (!(EwaldAccuracy > 0 && EwaldAccuracy < 1)) errors.Add("ewaldAccuracy must be between 0 and 1");
        if (TrajectoryEvery < 1) errors.Add("trajectoryEvery must be at least 1");
        if (LogEvery < 1) errors.Add("logEvery must be at least 1");
        if (string.IsNullOrWhiteSpace(Output)) errors.Add("output must be set");
        if (errors.Count > 0) throw new SettingsException("Invalid settings:\n  " + string.Join("\n  ", errors));
    }

    private static bool IsPositive(double x) => x > 0 && !double.IsInfinity(x);
}
