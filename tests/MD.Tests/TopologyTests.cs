using System;
using MD.Core;

namespace MD.Tests;

public static class TopologyTests
{
    // Placeholder three-site molecule with water's shape. These are NOT SPC/Fw parameters, which
    // must come from the paper (WATER_PLAN.md, W4). Only counts, indices and sums are tested here.
    private static readonly AtomType Ox = new("OX", "O", mass: 16.0, charge: -0.8, ljSigma: 3.0, ljEpsilon: 0.6);
    private static readonly AtomType Hx = new("HX", "H", mass: 1.0, charge: 0.4);

    private static MoleculeTemplate WaterShaped() => new(
        "X3",
        new[] { Ox, Hx, Hx },
        new[] { new BondTerm(0, 1), new BondTerm(0, 2) },
        new[] { new AngleTerm(1, 0, 2) });

    [Test]
    public static void WaterShapedMolecules_HaveExpectedCountsAndIndices()
    {
        var template = WaterShaped();
        var top = new Topology();
        const int n = 5;
        for (int m = 0; m < n; m++) Assert.Equal(m, top.AddMolecule(template), "molecule index");

        Assert.Equal(3 * n, top.AtomCount, "atoms");
        Assert.Equal(n, top.MoleculeCount, "molecules");
        Assert.Equal(2 * n, top.Bonds.Count, "bonds");
        Assert.Equal(n, top.Angles.Count, "angles");

        for (int m = 0; m < n; m++)
        {
            int o = top.FirstAtomOf(m);
            Assert.Equal(3 * m, o, "first atom");
            Assert.Equal(new BondTerm(o, o + 1), top.Bonds[2 * m], "first bond");
            Assert.Equal(new BondTerm(o, o + 2), top.Bonds[2 * m + 1], "second bond");
            Assert.Equal(new AngleTerm(o + 1, o, o + 2), top.Angles[m], "angle (vertex on O)");
            Assert.Equal("O", top.TypeOf(o).Element, "element of first atom");
            Assert.Equal("H", top.TypeOf(o + 2).Element, "element of third atom");
            for (int k = 0; k < 3; k++) Assert.Equal(m, top.MoleculeOf(o + k), "molecule of atom");
        }

        Assert.Close(0.0, template.TotalCharge(), 1e-15, 0, "molecule charge");
        Assert.Close(0.0, top.TotalCharge(), 1e-14, 0, "system charge");
        var q = top.Charges();
        var mass = top.Masses();
        var el = top.Elements();
        Assert.Equal(3 * n, q.Length, "charges length");
        Assert.Close(-0.8, q[3], 0, 0, "charge of O in molecule 1");
        Assert.Close(1.0, mass[4], 0, 0, "mass of H in molecule 1");
        Assert.Equal("H", el[5], "element of H in molecule 1");
        Assert.True(top.TypeOf(0).HasLj && !top.TypeOf(1).HasLj, "only O has an LJ site");
    }

    [Test]
    public static void Exclusions_AreExactlyTheIntramolecularPairs()
    {
        var top = new Topology();
        for (int m = 0; m < 4; m++) top.AddMolecule(WaterShaped());
        for (int i = 0; i < top.AtomCount; i++)
        {
            for (int j = 0; j < top.AtomCount; j++)
            {
                if (i == j) continue;
                bool sameMolecule = i / 3 == j / 3;
                Assert.Equal(sameMolecule, top.IsExcluded(i, j), $"exclusion of ({i}, {j})");
            }
        }
        // For a three-site molecule, intramolecular pairs are exactly the 1-2 (bond) and 1-3 (angle ends) pairs.
        int excludedPairs = 0;
        for (int i = 0; i < top.AtomCount; i++)
            for (int j = i + 1; j < top.AtomCount; j++)
                if (top.IsExcluded(i, j)) excludedPairs++;
        Assert.Equal(top.Bonds.Count + top.Angles.Count, excludedPairs, "excluded pairs = bonds + angle end pairs");
    }

    [Test]
    public static void ParticleSystem_BuildsFromTopology()
    {
        var top = new Topology();
        top.AddMolecule(WaterShaped());
        top.AddMolecule(WaterShaped());
        var positions = new Vec3[top.AtomCount];
        for (int i = 0; i < positions.Length; i++) positions[i] = new Vec3(i, 1, 1);
        var s = new ParticleSystem(new Box(new Vec3(20, 20, 20)), positions, top.Masses(), top.Elements());
        Assert.Equal(6, s.Count, "atom count");
        Assert.Close(2 * 18.0, s.TotalMass(), 0, 1e-15, "total mass");
    }

    [Test]
    public static void InvalidDefinitions_AreRejected()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => new MoleculeTemplate("bad", new[] { Ox, Hx }, new[] { new BondTerm(0, 2) }, Array.Empty<AngleTerm>()), "bond index out of range");
        Assert.Throws<ArgumentException>(() => new MoleculeTemplate("bad", new[] { Ox, Hx }, new[] { new BondTerm(1, 1) }, Array.Empty<AngleTerm>()), "self bond");
        Assert.Throws<ArgumentException>(() => new MoleculeTemplate("bad", new[] { Ox, Hx, Hx }, Array.Empty<BondTerm>(), new[] { new AngleTerm(1, 0, 1) }), "angle repeats an atom");
        Assert.Throws<ArgumentException>(() => new MoleculeTemplate("bad", Array.Empty<AtomType>(), Array.Empty<BondTerm>(), Array.Empty<AngleTerm>()), "empty molecule");
        Assert.Throws<ArgumentOutOfRangeException>(() => new AtomType("Z", "Z", mass: 0, charge: 0), "zero mass");
        Assert.Throws<ArgumentException>(() => new AtomType("Z", "Z", mass: 1, charge: 0, ljSigma: 3.0), "sigma without epsilon");
    }
}
