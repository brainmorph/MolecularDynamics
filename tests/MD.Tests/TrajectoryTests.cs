using System;
using System.IO;
using System.Linq;
using MD.Core;
using MD.Core.Trajectory;

namespace MD.Tests;

public static class TrajectoryTests
{
    private static TrajectoryHeader SampleHeader(int n) => new(
        TrajectoryFormat.Version, TrajectoryFormat.ReducedLjUnits,
        new UnitScales(3.405, 119.8, 39.948, 2156.0),
        Enumerable.Range(0, n).Select(i => i % 2 == 0 ? "O" : "H").ToArray(),
        Enumerable.Repeat(1.0, n).ToArray(),
        0.005, 1.25, -3.5, -2.25, "{\"cells\":[2,1,1]}");

    private static TrajectoryFrame SampleFrame(int n, long step, int bonds = 0) => new(
        step, step * 0.005, 1.0 + step, -2.0 - step, -1.0,
        new Vec3(3.4, 1.7, 1.7),
        Enumerable.Range(0, n).Select(i => new Vec3(i * 0.1, i * 0.2 + 1e-17, Math.PI / (i + 1))).ToArray(),
        Enumerable.Range(0, n).Select(i => new Int3(i, -i, int.MaxValue - i)).ToArray(),
        Enumerable.Range(0, bonds).Select(i => new Bond((uint)i, (uint)(n - 1 - i))).ToArray());

    [Test]
    public static void RoundTrip_PreservesEverythingExactly()
    {
        var ms = new MemoryStream();
        var header = SampleHeader(8);
        var frames = new[] { SampleFrame(8, 0), SampleFrame(8, 10, bonds: 2), SampleFrame(8, 20) };
        using (var w = new TrajectoryWriter(ms, leaveOpen: true))
        {
            w.WriteHeader(header);
            foreach (var f in frames) w.WriteFrame(f);
        }

        ms.Position = 0;
        using var r = new TrajectoryReader(ms);
        var h = r.Header;
        Assert.Equal(header.Version, h.Version, "version");
        Assert.Equal(header.Units, h.Units, "units");
        Assert.Equal(header.Scales, h.Scales, "scales");
        Assert.True(header.Elements.SequenceEqual(h.Elements), "elements");
        Assert.True(header.Masses.SequenceEqual(h.Masses), "masses");
        Assert.Equal(header.Timestep, h.Timestep, "timestep");
        Assert.Equal(header.InitialKineticEnergy, h.InitialKineticEnergy, "KE0");
        Assert.Equal(header.InitialPotentialEnergy, h.InitialPotentialEnergy, "PE0");
        Assert.Equal(header.InitialTotalEnergy, h.InitialTotalEnergy, "E0");
        Assert.Equal(header.SettingsJson, h.SettingsJson, "settings");

        var read = r.ReadFrames().ToList();
        Assert.Equal(frames.Length, read.Count, "frame count");
        for (int k = 0; k < frames.Length; k++)
        {
            var a = frames[k];
            var b = read[k];
            Assert.Equal(a.Step, b.Step, "step");
            Assert.Equal(a.Time, b.Time, "time");
            Assert.Equal(a.KineticEnergy, b.KineticEnergy, "KE");
            Assert.Equal(a.PotentialEnergy, b.PotentialEnergy, "PE");
            Assert.Equal(a.TotalEnergy, b.TotalEnergy, "E");
            Assert.Equal(a.BoxLengths, b.BoxLengths, "box");
            Assert.True(a.Positions.SequenceEqual(b.Positions), "positions bit-exact");
            Assert.True(a.Images.SequenceEqual(b.Images), "images");
            Assert.True(a.Bonds.SequenceEqual(b.Bonds), "bonds");
        }
    }

    /// <summary>Pins the byte layout so an accidental format change breaks a test (and needs a version bump).</summary>
    [Test]
    public static void ByteLayout_MatchesSpec()
    {
        var ms = new MemoryStream();
        using (var w = new TrajectoryWriter(ms, leaveOpen: true))
        {
            w.WriteHeader(SampleHeader(2) with { SettingsJson = "{}" });
            w.WriteFrame(SampleFrame(2, 5, bonds: 1));
        }
        byte[] bytes = ms.ToArray();
        // Header: magic 8 + version 4 + units (4+10) + scales 32 + dt 8 + energies 24 + count 4
        //         + 2 atoms × (string 4+1 + mass 8) + settings (4+2) = 120
        int headerSize = 8 + 4 + 14 + 32 + 8 + 24 + 4 + 2 * 13 + 6;
        // Frame: marker 4 + step 8 + time 8 + energies 24 + box 24 + positions 2×24 + images 2×12 + bond count 4 + 1 bond × 8
        int frameSize = 4 + 8 + 8 + 24 + 24 + 48 + 24 + 4 + 8;
        Assert.Equal(headerSize + frameSize, bytes.Length, "file size");
        Assert.True(bytes.AsSpan(0, 8).SequenceEqual("LJMDTRAJ"u8), "magic");
        Assert.Equal(1u, BitConverter.ToUInt32(bytes, 8), "version (little-endian)");
        Assert.True(bytes.AsSpan(headerSize, 4).SequenceEqual("FRAM"u8), "frame marker");
        Assert.Equal(5L, BitConverter.ToInt64(bytes, headerSize + 4), "step");
    }

    [Test]
    public static void BadMagic_IsRejected()
    {
        var ms = new MemoryStream(new byte[64]);
        Assert.Throws<TrajectoryFormatException>(() => new TrajectoryReader(ms), "bad magic");
    }

    [Test]
    public static void WrongVersion_IsRejected()
    {
        var ms = new MemoryStream();
        using (var w = new TrajectoryWriter(ms, leaveOpen: true)) w.WriteHeader(SampleHeader(1) with { Version = 2 });
        ms.Position = 0;
        Assert.Throws<TrajectoryFormatException>(() => new TrajectoryReader(ms), "version 2");
    }

    [Test]
    public static void TruncatedFrame_IsReported()
    {
        var ms = new MemoryStream();
        using (var w = new TrajectoryWriter(ms, leaveOpen: true))
        {
            w.WriteHeader(SampleHeader(4));
            w.WriteFrame(SampleFrame(4, 0));
            w.WriteFrame(SampleFrame(4, 1));
        }
        var cut = new MemoryStream(ms.ToArray().AsSpan(0, (int)ms.Length - 10).ToArray());
        using var r = new TrajectoryReader(cut);
        int count = 0;
        Assert.Throws<TrajectoryFormatException>(() => { foreach (var _ in r.ReadFrames()) count++; }, "truncated");
        Assert.Equal(1, count, "complete frames before the truncated one");
    }

    [Test]
    public static void Writer_RejectsWrongAtomCount()
    {
        using var w = new TrajectoryWriter(new MemoryStream());
        w.WriteHeader(SampleHeader(4));
        Assert.Throws<ArgumentException>(() => w.WriteFrame(SampleFrame(3, 0)), "frame with 3 atoms in a 4-atom file");
    }
}
