using System;
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

    /// <summary>A water run in real units with the topology's bonds in every frame, written by the core and read by the viewer.</summary>
    private static (MemoryStream Stream, Simulation Sim, Topology Topology, TrajectoryFrame[] Frames) WriteWaterRun()
    {
        var (box, x, top) = WaterBox.Cubic(2, 1.0, seed: 3);
        // Shift everything (a rigid translation) so that some molecules straddle the box boundary.
        for (int i = 0; i < x.Length; i++) x[i] += new Vec3(1.3, 1.3, 1.3);
        var s = new ParticleSystem(box, x, top.Masses(), top.Elements());
        VelocityInitializer.Assign(s, RealUnits.KelvinToEnergy(300.0), 3);
        double rc = 0.49 * box.Lengths.X;
        var sim = new Simulation(s, new VelocityVerlet(SpcFw.ForceField(top, rc, rc, 1e-6), RealUnits.FemtosecondsToTime(0.5)));
        Bond[] bonds = TrajectoryWriter.BondsFrom(top);
        var ms = new MemoryStream();
        var frames = new TrajectoryFrame[5];
        using (var w = new TrajectoryWriter(ms, leaveOpen: true))
        {
            w.WriteHeader(new TrajectoryHeader(
                TrajectoryFormat.Version, TrajectoryFormat.MolecularRealUnits, UnitScales.Real,
                top.Elements(), top.Masses(), RealUnits.FemtosecondsToTime(0.5),
                sim.Initial.KineticEnergy, sim.Initial.PotentialEnergy, sim.Initial.TotalEnergy, "{}"));
            for (int k = 0; k < frames.Length; k++)
            {
                frames[k] = TrajectoryWriter.FrameFrom(sim, bonds);
                w.WriteFrame(frames[k]);
                for (int j = 0; j < 40; j++) sim.Step();
            }
        }
        ms.Position = 0;
        return (ms, sim, top, frames);
    }

    [Test]
    public static void ViewerReader_ReadsWaterRunWithBondsAndRealUnits()
    {
        var (ms, sim, top, frames) = WriteWaterRun();
        var t = Viewer.TrajectoryFile.Load(ms);
        Assert.True(t.Warning == null, $"no warning expected, got: {t.Warning}");
        Assert.Equal("real", t.Units, "units");
        Assert.Equal(1.0, t.SigmaAngstrom, "length scale: 1 Å");
        Assert.Equal(1.0 / RealUnits.GasConstantKilojoulesPerMoleKelvin, t.EpsilonKelvin, "energy scale: 1 kJ/mol in K");
        Assert.Equal(1.0, t.MassAmu, "mass scale: 1 u");
        Assert.Equal(RealUnits.TimeUnitFemtoseconds, t.TauFemtoseconds, "time scale: 100 fs");
        Assert.Close(100.0, t.TauFemtoseconds, 0, 1e-15, "time unit is 100 fs");
        Assert.Equal(sim.Initial.TotalEnergy, t.InitialTotalEnergy, "E0");
        Assert.Equal(top.AtomCount, t.AtomCount, "atom count");
        for (int i = 0; i < t.AtomCount; i++)
        {
            Assert.Equal(top.TypeOf(i).Element, t.Elements[i], $"element {i}");
            Assert.Equal(top.TypeOf(i).Mass, t.Masses[i], $"mass {i}");
        }
        Assert.Equal(frames.Length, t.Frames.Count, "frame count");
        for (int k = 0; k < frames.Length; k++)
        {
            Assert.Equal(frames[k].TotalEnergy, t.Frames[k].TotalEnergy, $"frame {k} E");
            Assert.Equal(2 * top.Bonds.Count, t.Frames[k].Bonds.Length, $"frame {k} bond count");
            for (int j = 0; j < top.Bonds.Count; j++)
            {
                Assert.Equal((uint)top.Bonds[j].I, t.Frames[k].Bonds[2 * j], $"frame {k} bond {j} i");
                Assert.Equal((uint)top.Bonds[j].J, t.Frames[k].Bonds[2 * j + 1], $"frame {k} bond {j} j");
            }
            for (int i = 0; i < t.AtomCount; i++)
            {
                Assert.Equal(frames[k].Positions[i], new Vec3(t.Frames[k].Positions[3 * i], t.Frames[k].Positions[3 * i + 1], t.Frames[k].Positions[3 * i + 2]), $"frame {k} position {i}");
            }
        }
    }

    /// <summary>
    /// The viewer's layout draws each molecule whole: O (the lowest index of each bonded group) inside
    /// the box, and every bonded atom at the copy nearest its partner, so the drawn bond vector equals
    /// the core's minimum-image bond vector. The test water box has molecules that cross the boundary.
    /// </summary>
    [Test]
    public static void ViewerLayout_DrawsMoleculesWhole()
    {
        var (ms, _, top, _) = WriteWaterRun();
        var t = Viewer.TrajectoryFile.Load(ms);
        var layout = new Viewer.DisplayLayout();
        var p = new double[3 * t.AtomCount];
        int crossing = 0;
        for (int k = 0; k + 1 < t.Frames.Count; k++)
        {
            var a = t.Frames[k];
            var box = new Box(new Vec3(a.BoxX, a.BoxY, a.BoxZ));
            foreach (double alpha in new[] { 0.0, 0.5 })
            {
                layout.Positions(a, t.Frames[k + 1], alpha, t.AtomCount, p);
                for (int m = 0; m < top.MoleculeCount; m++)
                {
                    int o = top.FirstAtomOf(m);
                    var po = new Vec3(p[3 * o], p[3 * o + 1], p[3 * o + 2]);
                    Assert.True(po.X >= 0 && po.X < box.Lengths.X && po.Y >= 0 && po.Y < box.Lengths.Y && po.Z >= 0 && po.Z < box.Lengths.Z,
                        $"frame {k}, molecule {m}: O drawn outside the box");
                    for (int h = o + 1; h <= o + 2; h++)
                    {
                        var ph = new Vec3(p[3 * h], p[3 * h + 1], p[3 * h + 2]);
                        if (alpha == 0.0)
                        {
                            var wo = new Vec3(a.Positions[3 * o], a.Positions[3 * o + 1], a.Positions[3 * o + 2]);
                            var wh = new Vec3(a.Positions[3 * h], a.Positions[3 * h + 1], a.Positions[3 * h + 2]);
                            Vec3 expected = box.MinimumImage(wh - wo);
                            Assert.Close(0.0, (ph - po - expected).Length, 1e-12, 0, $"frame {k}, atom {h}: drawn O–H vector");
                            if ((wh - wo).Length > 2.0) crossing++;
                        }
                        // Between frames the blend is a straight line, which can shorten a fast O–H
                        // stretch (frames here are 20 fs apart, the stretch period about 9 fs). A split
                        // molecule would show a bond of about half the 6.2 Å box or more.
                        double r = (ph - po).Length;
                        Assert.True(r > 0.5 && r < 1.5, $"frame {k}, alpha {alpha}, atom {h}: drawn O–H length {r:F3} Å");
                    }
                }
            }
        }
        Assert.True(crossing > 0, "test setup: some molecule should cross the box boundary in its wrapped positions");
    }

    [Test]
    public static void ViewerLayout_WrapsUnbondedAtomsLikeTheCore()
    {
        var (ms, _, _) = WriteRun();
        var t = Viewer.TrajectoryFile.Load(ms);
        var layout = new Viewer.DisplayLayout();
        var p = new double[3 * t.AtomCount];
        var f = t.Frames[1]; // no bonds in this frame
        Assert.Equal(0, f.Bonds.Length, "test setup: frame without bonds");
        layout.Positions(f, f, 0.0, t.AtomCount, p);
        for (int i = 0; i < 3 * t.AtomCount; i++) Assert.Close(f.Positions[i], p[i], 1e-12, 0, $"coordinate {i}");
    }
}
