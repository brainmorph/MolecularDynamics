using System.IO;
using MD.Core;
using MD.Core.Trajectory;
using Viewer = MDViewer.Trajectory;

namespace MD.Tests;

/// <summary>
/// The Godot viewer has its own reader, written from docs/trajectory-format.md. These tests write
/// with the core and read with the viewer's reader, so a disagreement about the format fails here
/// instead of showing wrong numbers on screen.
/// </summary>
public static class ViewerConformanceTests
{
    private static (MemoryStream Stream, Simulation Sim, TrajectoryFrame[] Frames) WriteRun()
    {
        var s = Helpers.Fcc(2, 1, 1);
        for (int i = 0; i < s.Count; i++) s.Elements[i] = i % 3 == 0 ? "O" : "H";
        var sim = Helpers.NewSimulation(s, 0.722, 0.005);
        var ms = new MemoryStream();
        var frames = new TrajectoryFrame[6];
        using (var w = new TrajectoryWriter(ms, leaveOpen: true))
        {
            w.WriteHeader(new TrajectoryHeader(
                TrajectoryFormat.Version, TrajectoryFormat.ReducedLjUnits,
                new UnitScales(ArgonUnits.SigmaAngstrom, ArgonUnits.EpsilonOverBoltzmannKelvin, ArgonUnits.MassAtomicUnits, ArgonUnits.TauFemtoseconds),
                s.Elements, s.Masses, 0.005,
                sim.Initial.KineticEnergy, sim.Initial.PotentialEnergy, sim.Initial.TotalEnergy, "{\"x\":\"é\"}"));
            for (int k = 0; k < frames.Length; k++)
            {
                frames[k] = TrajectoryWriter.FrameFrom(sim);
                // Bonds are always empty in the LJ phase; put some in one frame to exercise the field.
                if (k == 2) frames[k] = frames[k] with { Bonds = new[] { new Bond(0, 1), new Bond(3, 7) } };
                w.WriteFrame(frames[k]);
                for (int j = 0; j < 300; j++) sim.Step(); // long enough for some atoms to cross the boundary
            }
        }
        ms.Position = 0;
        return (ms, sim, frames);
    }

    [Test]
    public static void ViewerReader_ReadsCoreOutputExactly()
    {
        var (ms, sim, frames) = WriteRun();
        var t = Viewer.TrajectoryFile.Load(ms);

        Assert.True(t.Warning == null, $"no warning expected, got: {t.Warning}");
        Assert.Equal(1u, t.Version, "version");
        Assert.Equal("lj-reduced", t.Units, "units");
        Assert.Equal(ArgonUnits.TauFemtoseconds, t.TauFemtoseconds, "tau");
        Assert.Equal(ArgonUnits.SigmaAngstrom, t.SigmaAngstrom, "sigma");
        Assert.Equal(sim.Initial.TotalEnergy, t.InitialTotalEnergy, "E0");
        Assert.Equal(sim.Initial.KineticEnergy, t.InitialKineticEnergy, "KE0");
        Assert.Equal(sim.Initial.PotentialEnergy, t.InitialPotentialEnergy, "PE0");
        Assert.Equal("{\"x\":\"é\"}", t.SettingsJson, "settings (UTF-8)");
        Assert.Equal(sim.System.Count, t.AtomCount, "atom count");
        for (int i = 0; i < t.AtomCount; i++) Assert.Equal(sim.System.Elements[i], t.Elements[i], $"element {i}");
        Assert.Equal(frames.Length, t.Frames.Count, "frame count");

        bool sawImage = false;
        for (int k = 0; k < frames.Length; k++)
        {
            var a = frames[k];
            var b = t.Frames[k];
            Assert.Equal(a.Step, b.Step, "step");
            Assert.Equal(a.Time, b.Time, "time");
            Assert.Equal(a.KineticEnergy, b.KineticEnergy, "KE");
            Assert.Equal(a.PotentialEnergy, b.PotentialEnergy, "PE");
            Assert.Equal(a.TotalEnergy, b.TotalEnergy, "E");
            Assert.Equal(a.BoxLengths, new Vec3(b.BoxX, b.BoxY, b.BoxZ), "box");
            for (int i = 0; i < t.AtomCount; i++)
            {
                Assert.Equal(a.Positions[i], new Vec3(b.Positions[3 * i], b.Positions[3 * i + 1], b.Positions[3 * i + 2]), $"position {i}");
                Assert.Equal(a.Images[i], new Int3(b.Images[3 * i], b.Images[3 * i + 1], b.Images[3 * i + 2]), $"image {i}");
                var unwrapped = new Box(a.BoxLengths).Unwrap(a.Positions[i], a.Images[i]);
                Assert.Equal(unwrapped, new Vec3(b.Unwrapped(i, 0), b.Unwrapped(i, 1), b.Unwrapped(i, 2)), $"unwrapped {i}");
                sawImage |= a.Images[i] != Int3.Zero;
            }
            Assert.Equal(2 * a.Bonds.Length, b.Bonds.Length, "bond count");
            for (int j = 0; j < a.Bonds.Length; j++)
            {
                Assert.Equal(a.Bonds[j].I, b.Bonds[2 * j], "bond i");
                Assert.Equal(a.Bonds[j].J, b.Bonds[2 * j + 1], "bond j");
            }
        }
        Assert.True(sawImage, "test setup: some frame should have a non-zero image count");
    }

    [Test]
    public static void ViewerReader_KeepsCompleteFramesOfTruncatedFile()
    {
        var (ms, _, frames) = WriteRun();
        byte[] bytes = ms.ToArray();
        var cut = new MemoryStream(bytes, 0, bytes.Length - 5);
        var t = Viewer.TrajectoryFile.Load(cut);
        Assert.Equal(frames.Length - 1, t.Frames.Count, "complete frames kept");
        Assert.True(t.Warning != null, "truncation reported");
    }

    [Test]
    public static void ViewerReader_RejectsBadMagicAndVersion()
    {
        Assert.Throws<InvalidDataException>(() => Viewer.TrajectoryFile.Load(new MemoryStream(new byte[100])), "bad magic");
        var (ms, _, _) = WriteRun();
        byte[] bytes = ms.ToArray();
        bytes[8] = 2; // version field
        Assert.Throws<InvalidDataException>(() => Viewer.TrajectoryFile.Load(new MemoryStream(bytes)), "version 2");
    }
}
