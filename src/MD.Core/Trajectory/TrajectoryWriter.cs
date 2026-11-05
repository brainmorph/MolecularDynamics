using System;
using System.IO;
using System.Text;

namespace MD.Core.Trajectory;

/// <summary>Writes trajectory format version 1 to a stream. All numbers are little-endian.</summary>
public sealed class TrajectoryWriter : IDisposable
{
    private readonly BinaryWriter _w;
    private int _atomCount = -1;

    public TrajectoryWriter(Stream stream, bool leaveOpen = false)
    {
        ArgumentNullException.ThrowIfNull(stream);
        _w = new BinaryWriter(stream, new UTF8Encoding(false, true), leaveOpen);
    }

    public void WriteHeader(TrajectoryHeader h)
    {
        ArgumentNullException.ThrowIfNull(h);
        if (_atomCount >= 0) throw new InvalidOperationException("Header already written.");
        if (h.Masses.Count != h.Elements.Count) throw new ArgumentException("Masses and elements must have the same length.");

        _w.Write(TrajectoryFormat.Magic);
        _w.Write(h.Version);
        WriteString(h.Units);
        _w.Write(h.Scales.SigmaAngstrom);
        _w.Write(h.Scales.EpsilonKelvin);
        _w.Write(h.Scales.MassAmu);
        _w.Write(h.Scales.TauFemtoseconds);
        _w.Write(h.Timestep);
        _w.Write(h.InitialKineticEnergy);
        _w.Write(h.InitialPotentialEnergy);
        _w.Write(h.InitialTotalEnergy);
        _w.Write((uint)h.AtomCount);
        for (int i = 0; i < h.AtomCount; i++)
        {
            WriteString(h.Elements[i]);
            _w.Write(h.Masses[i]);
        }
        WriteString(h.SettingsJson);
        _atomCount = h.AtomCount;
    }

    public void WriteFrame(TrajectoryFrame f)
    {
        ArgumentNullException.ThrowIfNull(f);
        if (_atomCount < 0) throw new InvalidOperationException("Write the header before any frame.");
        if (f.Positions.Length != _atomCount || f.Images.Length != _atomCount)
        {
            throw new ArgumentException($"Frame must have {_atomCount} positions and images.");
        }

        _w.Write(TrajectoryFormat.FrameMarker);
        _w.Write(f.Step);
        _w.Write(f.Time);
        _w.Write(f.KineticEnergy);
        _w.Write(f.PotentialEnergy);
        _w.Write(f.TotalEnergy);
        WriteVec3(f.BoxLengths);
        foreach (Vec3 p in f.Positions) WriteVec3(p);
        foreach (Int3 im in f.Images)
        {
            _w.Write(im.X);
            _w.Write(im.Y);
            _w.Write(im.Z);
        }
        _w.Write((uint)f.Bonds.Length);
        foreach (Bond b in f.Bonds)
        {
            if (b.I >= _atomCount || b.J >= _atomCount) throw new ArgumentException("Bond index out of range.");
            _w.Write(b.I);
            _w.Write(b.J);
        }
    }

    /// <summary>Builds a frame from the current state of a simulation.</summary>
    public static TrajectoryFrame FrameFrom(Simulation sim)
    {
        ThermoState t = sim.State;
        ParticleSystem s = sim.System;
        return new TrajectoryFrame(
            t.Step, t.Time, t.KineticEnergy, t.PotentialEnergy, t.TotalEnergy,
            s.Box.Lengths, (Vec3[])s.Positions.Clone(), (Int3[])s.Images.Clone(), Array.Empty<Bond>());
    }

    public void Flush() => _w.Flush();

    public void Dispose() => _w.Dispose();

    private void WriteVec3(Vec3 v)
    {
        _w.Write(v.X);
        _w.Write(v.Y);
        _w.Write(v.Z);
    }

    private void WriteString(string s)
    {
        byte[] bytes = Encoding.UTF8.GetBytes(s);
        _w.Write((uint)bytes.Length);
        _w.Write(bytes);
    }
}
