using System;
using System.Collections.Generic;
using System.IO;
using System.Text;

namespace MD.Core.Trajectory;

/// <summary>
/// Reads trajectory format version 1. Used by the dump tool and tests. (The Godot viewer has its own
/// reader written from the format spec; it does not use this code.)
/// </summary>
public sealed class TrajectoryReader : IDisposable
{
    // Guards against absurd lengths from a corrupt file.
    private const uint MaxStringBytes = 16 * 1024 * 1024;
    private const uint MaxAtoms = 100_000_000;

    private readonly BinaryReader _r;
    private readonly Stream _stream;

    public TrajectoryReader(Stream stream, bool leaveOpen = false)
    {
        _stream = stream ?? throw new ArgumentNullException(nameof(stream));
        _r = new BinaryReader(stream, new UTF8Encoding(false, true), leaveOpen);
        Header = ReadHeader();
    }

    public TrajectoryHeader Header { get; }

    /// <summary>Reads frames until end of file. Throws if the file ends partway through a frame.</summary>
    public IEnumerable<TrajectoryFrame> ReadFrames()
    {
        while (true)
        {
            byte[] marker = _r.ReadBytes(4);
            if (marker.Length == 0) yield break;
            if (marker.Length < 4 || !marker.AsSpan().SequenceEqual(TrajectoryFormat.FrameMarker))
            {
                throw new TrajectoryFormatException($"Expected frame marker at byte {_stream.Position - marker.Length}.");
            }
            yield return ReadFrameBody();
        }
    }

    public void Dispose() => _r.Dispose();

    private TrajectoryHeader ReadHeader()
    {
        try
        {
            byte[] magic = _r.ReadBytes(8);
            if (!magic.AsSpan().SequenceEqual(TrajectoryFormat.Magic))
            {
                throw new TrajectoryFormatException("Not a trajectory file (bad magic).");
            }
            uint version = _r.ReadUInt32();
            if (version != TrajectoryFormat.Version)
            {
                throw new TrajectoryFormatException($"Unsupported trajectory version {version} (reader supports {TrajectoryFormat.Version}).");
            }
            string units = ReadString();
            var scales = new UnitScales(_r.ReadDouble(), _r.ReadDouble(), _r.ReadDouble(), _r.ReadDouble());
            double dt = _r.ReadDouble();
            double ke0 = _r.ReadDouble();
            double pe0 = _r.ReadDouble();
            double e0 = _r.ReadDouble();
            uint n = _r.ReadUInt32();
            if (n > MaxAtoms) throw new TrajectoryFormatException($"Atom count {n} is not plausible.");
            var elements = new string[n];
            var masses = new double[n];
            for (int i = 0; i < n; i++)
            {
                elements[i] = ReadString();
                masses[i] = _r.ReadDouble();
            }
            string settings = ReadString();
            return new TrajectoryHeader(version, units, scales, elements, masses, dt, ke0, pe0, e0, settings);
        }
        catch (EndOfStreamException e)
        {
            throw new TrajectoryFormatException("File ends inside the header.", e);
        }
    }

    private TrajectoryFrame ReadFrameBody()
    {
        try
        {
            int n = Header.AtomCount;
            long step = _r.ReadInt64();
            double time = _r.ReadDouble();
            double ke = _r.ReadDouble();
            double pe = _r.ReadDouble();
            double e = _r.ReadDouble();
            Vec3 box = ReadVec3();
            var pos = new Vec3[n];
            for (int i = 0; i < n; i++) pos[i] = ReadVec3();
            var images = new Int3[n];
            for (int i = 0; i < n; i++) images[i] = new Int3(_r.ReadInt32(), _r.ReadInt32(), _r.ReadInt32());
            uint bondCount = _r.ReadUInt32();
            if (bondCount > MaxAtoms) throw new TrajectoryFormatException($"Bond count {bondCount} is not plausible.");
            var bonds = new Bond[bondCount];
            for (int i = 0; i < bondCount; i++)
            {
                var b = new Bond(_r.ReadUInt32(), _r.ReadUInt32());
                if (b.I >= n || b.J >= n) throw new TrajectoryFormatException($"Bond {i} in step {step} refers to a missing atom.");
                bonds[i] = b;
            }
            return new TrajectoryFrame(step, time, ke, pe, e, box, pos, images, bonds);
        }
        catch (EndOfStreamException ex)
        {
            throw new TrajectoryFormatException("File ends partway through a frame (truncated).", ex);
        }
    }

    private Vec3 ReadVec3() => new(_r.ReadDouble(), _r.ReadDouble(), _r.ReadDouble());

    private string ReadString()
    {
        uint len = _r.ReadUInt32();
        if (len > MaxStringBytes) throw new TrajectoryFormatException($"String length {len} is not plausible.");
        byte[] bytes = _r.ReadBytes((int)len);
        if (bytes.Length != len) throw new EndOfStreamException();
        return Encoding.UTF8.GetString(bytes);
    }
}
