// Reader for trajectory format version 1, written from docs/trajectory-format.md only.
// It must not use MD.Core, and it uses no Godot types, so the core test suite can
// compile this file to check that it agrees with the core's writer.
using System;
using System.Collections.Generic;
using System.IO;
using System.Text;

namespace MDViewer.Trajectory;

public sealed class TrajectoryFrame
{
    public long Step;
    public double Time;
    public double KineticEnergy;
    public double PotentialEnergy;
    public double TotalEnergy;
    public double BoxX, BoxY, BoxZ;

    /// <summary>Wrapped positions, x0 y0 z0 x1 y1 z1 ... (σ).</summary>
    public double[] Positions = Array.Empty<double>();

    /// <summary>Image counts, same layout as <see cref="Positions"/>.</summary>
    public int[] Images = Array.Empty<int>();

    /// <summary>Bond atom-index pairs, i0 j0 i1 j1 ...</summary>
    public uint[] Bonds = Array.Empty<uint>();

    /// <summary>Unwrapped coordinate: position + image × box length.</summary>
    public double Unwrapped(int atom, int axis)
    {
        double box = axis == 0 ? BoxX : axis == 1 ? BoxY : BoxZ;
        return Positions[3 * atom + axis] + Images[3 * atom + axis] * box;
    }
}

public sealed class TrajectoryFile
{
    public uint Version;
    public string Units = "";
    public double SigmaAngstrom, EpsilonKelvin, MassAmu, TauFemtoseconds;
    public double Timestep;
    public double InitialKineticEnergy, InitialPotentialEnergy, InitialTotalEnergy;
    public string[] Elements = Array.Empty<string>();
    public double[] Masses = Array.Empty<double>();
    public string SettingsJson = "";
    public List<TrajectoryFrame> Frames = new();

    /// <summary>Set when the file ends partway through a frame. Complete frames are kept.</summary>
    public string? Warning;

    public int AtomCount => Elements.Length;

    public const uint SupportedVersion = 1;

    public static TrajectoryFile Load(string path)
    {
        using var stream = File.OpenRead(path);
        return Load(stream);
    }

    /// <summary>Reads the whole file. Throws <see cref="InvalidDataException"/> for a bad header.</summary>
    public static TrajectoryFile Load(Stream stream)
    {
        // BinaryReader is always little-endian, as the spec requires.
        using var r = new BinaryReader(stream, Encoding.UTF8, leaveOpen: true);
        var t = new TrajectoryFile();
        try
        {
            if (Encoding.ASCII.GetString(r.ReadBytes(8)) != "LJMDTRAJ")
                throw new InvalidDataException("Not a trajectory file (bad magic).");
            t.Version = r.ReadUInt32();
            if (t.Version != SupportedVersion)
                throw new InvalidDataException($"Trajectory version {t.Version} is not supported (viewer supports {SupportedVersion}).");
            t.Units = ReadString(r);
            t.SigmaAngstrom = r.ReadDouble();
            t.EpsilonKelvin = r.ReadDouble();
            t.MassAmu = r.ReadDouble();
            t.TauFemtoseconds = r.ReadDouble();
            t.Timestep = r.ReadDouble();
            t.InitialKineticEnergy = r.ReadDouble();
            t.InitialPotentialEnergy = r.ReadDouble();
            t.InitialTotalEnergy = r.ReadDouble();
            uint n = r.ReadUInt32();
            if (n > 100_000_000) throw new InvalidDataException($"Atom count {n} is not plausible.");
            t.Elements = new string[n];
            t.Masses = new double[n];
            for (int i = 0; i < n; i++)
            {
                t.Elements[i] = ReadString(r);
                t.Masses[i] = r.ReadDouble();
            }
            t.SettingsJson = ReadString(r);
        }
        catch (EndOfStreamException)
        {
            throw new InvalidDataException("File ends inside the header.");
        }

        int atoms = t.AtomCount;
        while (true)
        {
            byte[] marker = r.ReadBytes(4);
            if (marker.Length == 0) break;
            if (Encoding.ASCII.GetString(marker) != "FRAM")
            {
                t.Warning = $"Bad frame marker after frame {t.Frames.Count}; the rest of the file was ignored.";
                break;
            }
            try
            {
                var f = new TrajectoryFrame
                {
                    Step = r.ReadInt64(),
                    Time = r.ReadDouble(),
                    KineticEnergy = r.ReadDouble(),
                    PotentialEnergy = r.ReadDouble(),
                    TotalEnergy = r.ReadDouble(),
                    BoxX = r.ReadDouble(),
                    BoxY = r.ReadDouble(),
                    BoxZ = r.ReadDouble(),
                    Positions = new double[3 * atoms],
                    Images = new int[3 * atoms],
                };
                for (int k = 0; k < 3 * atoms; k++) f.Positions[k] = r.ReadDouble();
                for (int k = 0; k < 3 * atoms; k++) f.Images[k] = r.ReadInt32();
                uint bonds = r.ReadUInt32();
                if (bonds > 100_000_000) throw new InvalidDataException("Bond count is not plausible.");
                f.Bonds = new uint[2 * bonds];
                for (int k = 0; k < 2 * bonds; k++)
                {
                    f.Bonds[k] = r.ReadUInt32();
                    if (f.Bonds[k] >= atoms) throw new InvalidDataException("Bond refers to a missing atom.");
                }
                t.Frames.Add(f);
            }
            catch (Exception e) when (e is EndOfStreamException or InvalidDataException)
            {
                t.Warning = $"File is truncated or damaged after frame {t.Frames.Count} ({e.Message}). Complete frames were kept.";
                break;
            }
        }
        return t;
    }

    private static string ReadString(BinaryReader r)
    {
        uint len = r.ReadUInt32();
        if (len > 16 * 1024 * 1024) throw new InvalidDataException("String length is not plausible.");
        byte[] bytes = r.ReadBytes((int)len);
        if (bytes.Length != len) throw new EndOfStreamException();
        return Encoding.UTF8.GetString(bytes);
    }
}
