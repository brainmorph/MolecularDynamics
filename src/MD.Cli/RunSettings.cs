using System;
using System.Collections.Generic;
using System.Text.Json;
using System.Text.Json.Serialization;

namespace MD.Cli;

/// <summary>
/// Run settings, read from a JSON file. All physical values are in reduced LJ units.
/// Unknown keys are rejected so typos do not silently fall back to defaults.
/// </summary>
public sealed class RunSettings
{
    /// <summary>FCC cells along x, y, z. Atom count = 4·x·y·z. [2,1,1] gives 8 atoms.</summary>
    public int[] Cells { get; set; } = { 2, 1, 1 };

    /// <summary>
    /// Number density ρσ³. Without <see cref="Box"/> this is the system's density. With <see cref="Box"/>
    /// it only sets the lattice spacing (a = (4/ρ)^(1/3)); the system's density is then atoms / box volume.
    /// </summary>
    public double Density { get; set; } = 0.8442;

    /// <summary>
    /// Optional box lengths (x, y, z) in σ. When set, the lattice is centred in this box, leaving the rest
    /// empty. Use it to start a compact cluster in a large box. Omit it to fill the box at Density.
    /// </summary>
    public double[]? Box { get; set; }

    /// <summary>Starting temperature (k_B T / ε).</summary>
    public double Temperature { get; set; } = 0.722;

    /// <summary>Timestep in τ.</summary>
    public double Timestep { get; set; } = 0.005;

    public long Steps { get; set; } = 20_000;
    public int Seed { get; set; } = 12345;
    public double Cutoff { get; set; } = 2.5;

    /// <summary>Write a trajectory frame every N steps (frame 0 is always written).</summary>
    public int TrajectoryEvery { get; set; } = 10;

    /// <summary>Write a log line every N steps.</summary>
    public int LogEvery { get; set; } = 100;

    /// <summary>Element labels, assigned to atoms in turn. Labels only: all atoms use the same argon-like LJ parameters.</summary>
    public string[] Elements { get; set; } = { "Ar" };

    /// <summary>Trajectory path, relative to the settings file.</summary>
    public string Output { get; set; } = "out/run.traj";

    /// <summary>Log path, relative to the settings file. Defaults to the output path with ".log".</summary>
    public string? Log { get; set; }

    public static readonly JsonSerializerOptions JsonOptions = new()
    {
        PropertyNamingPolicy = JsonNamingPolicy.CamelCase,
        UnmappedMemberHandling = JsonUnmappedMemberHandling.Disallow,
        ReadCommentHandling = JsonCommentHandling.Skip,
        AllowTrailingCommas = true,
        WriteIndented = true,
    };

    public static RunSettings Parse(string json)
    {
        RunSettings s = JsonSerializer.Deserialize<RunSettings>(json, JsonOptions)
            ?? throw new SettingsException("Settings file is empty.");
        s.Validate();
        return s;
    }

    public string ToJson() => JsonSerializer.Serialize(this, JsonOptions);

    [JsonIgnore]
    public int AtomCount => 4 * Cells[0] * Cells[1] * Cells[2];

    public void Validate()
    {
        var errors = new List<string>();
        if (Cells is not { Length: 3 } || Cells[0] < 1 || Cells[1] < 1 || Cells[2] < 1)
            errors.Add("cells must be three integers, each at least 1");
        else if ((long)Cells[0] * Cells[1] * Cells[2] > 1_000_000)
            errors.Add("cells gives too many atoms");
        if (!IsPositive(Density)) errors.Add("density must be positive");
        if (Box is not null && (Box.Length != 3 || Array.Exists(Box, x => !IsPositive(x))))
            errors.Add("box must be three positive lengths");
        if (!(Temperature >= 0) || double.IsInfinity(Temperature)) errors.Add("temperature must be zero or positive");
        if (!IsPositive(Timestep)) errors.Add("timestep must be positive");
        if (Steps < 0) errors.Add("steps must be zero or positive");
        if (!IsPositive(Cutoff)) errors.Add("cutoff must be positive");
        if (TrajectoryEvery < 1) errors.Add("trajectoryEvery must be at least 1");
        if (LogEvery < 1) errors.Add("logEvery must be at least 1");
        if (Elements is not { Length: > 0 } || Array.Exists(Elements, e => string.IsNullOrWhiteSpace(e)))
            errors.Add("elements must be a non-empty list of non-empty labels");
        if (string.IsNullOrWhiteSpace(Output)) errors.Add("output must be set");
        if (errors.Count > 0) throw new SettingsException("Invalid settings:\n  " + string.Join("\n  ", errors));
    }

    private static bool IsPositive(double x) => x > 0 && !double.IsInfinity(x);
}

public sealed class SettingsException : Exception
{
    public SettingsException(string message) : base(message) { }
}
