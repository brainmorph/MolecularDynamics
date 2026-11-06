using System;

namespace MD.Core;

/// <summary>
/// SPC/Fw flexible water (Wu, Tepper &amp; Voth, J. Chem. Phys. 124, 024503 (2006)), decision D5.
///
/// The primary paper could not be read (paywalled). Values were taken from sources that cite it, and
/// are recorded with those sources in WATER_PLAN.md Section 7:
/// - Bonded terms, U = ½k_b(r − r0)² and U = ½k_a(θ − θ0)²: k_b = 1059.162 kcal/(mol·Å²),
///   r0 = 1.012 Å, k_a = 75.90 kcal/(mol·rad²), θ0 = 113.24° — Faulkner, Qin, Maggs &amp; Krauth,
///   arXiv:1804.05795, Section V.A (its unit labels for k_b and k_a are swapped; dimensions fix them).
/// - Charges q_O = −0.82 e, q_H = +0.41 e — same source; equal to SPC.
/// - O–O LJ: the SPC values (SPC/Fw keeps SPC's LJ and charges, Raabe &amp; Sadus 2007).
///   Sources differ in the last digit: ε = 0.1553 (LAMMPS SPC page), 0.155 (arXiv, 4ε = 0.62),
///   0.650 kJ/mol ≈ 0.1554 kcal/mol; σ = 3.165 to 3.166 Å. Used here: ε = 0.1554 kcal/mol,
///   σ = 3.1655 Å. **The exact digits in the paper's table are still unverified.**
/// - Masses: O 15.9994, H 1.008 g/mol (LAMMPS SPC page; standard atomic weights).
///
/// Energies are converted with the thermochemical calorie, 1 kcal = 4.184 kJ (exact).
/// No hydrogen LJ site. All intramolecular non-bonded pairs are excluded (Topology rule).
/// </summary>
public static class SpcFw
{
    public const double KilojoulesPerKilocalorie = 4.184;

    public const double BondForceConstantKcal = 1059.162;
    public const double BondLengthAngstrom = 1.012;
    public const double AngleForceConstantKcal = 75.90;
    public const double AngleDegrees = 113.24;
    public const double ChargeO = -0.82;
    public const double ChargeH = 0.41;
    public const double LjEpsilonKcal = 0.1554;
    public const double LjSigmaAngstrom = 3.1655;
    public const double MassO = 15.9994;
    public const double MassH = 1.008;

    public static double MolarMass => MassO + 2 * MassH;

    public static readonly AtomType Oxygen = new("OW", "O", MassO, ChargeO,
        ljSigma: LjSigmaAngstrom, ljEpsilon: LjEpsilonKcal * KilojoulesPerKilocalorie);
    public static readonly AtomType Hydrogen = new("HW", "H", MassH, ChargeH);

    public static readonly HarmonicBond Bond = new(BondForceConstantKcal * KilojoulesPerKilocalorie, BondLengthAngstrom);
    public static readonly HarmonicAngle Angle = new(AngleForceConstantKcal * KilojoulesPerKilocalorie, AngleDegrees * Math.PI / 180.0);

    /// <summary>Atom order O, H, H; bonds O–H1, O–H2; angle H1–O–H2.</summary>
    public static readonly MoleculeTemplate Molecule = new(
        "SPC/Fw",
        new[] { Oxygen, Hydrogen, Hydrogen },
        new[] { new BondTerm(0, 1, Bond), new BondTerm(0, 2, Bond) },
        new[] { new AngleTerm(1, 0, 2, Angle) });

    /// <summary>Equilibrium geometry (r0, θ0) with O at <paramref name="o"/>, oriented by the orthonormal pair (u, w):
    /// the H atoms lie in the u–w plane, symmetric about u.</summary>
    public static Vec3[] EquilibriumGeometry(Vec3 o, Vec3 u, Vec3 w)
    {
        double half = 0.5 * Angle.Angle;
        double r = Bond.Length;
        return new[]
        {
            o,
            o + r * (Math.Cos(half) * u + Math.Sin(half) * w),
            o + r * (Math.Cos(half) * u - Math.Sin(half) * w),
        };
    }

    /// <summary>
    /// The full SPC/Fw force field for a topology: bonded terms, O–O LJ (truncated and energy-shifted
    /// at <paramref name="ljCutoff"/>), and Ewald Coulomb with the given real-space cutoff and accuracy.
    /// </summary>
    public static CompositeForceField ForceField(Topology topology, double ljCutoff, double ewaldRealCutoff, double ewaldAccuracy) =>
        new(new BondedForces(topology),
            new MolecularLennardJones(topology, ljCutoff),
            new EwaldForceField(Ewald.ForAccuracy(ewaldRealCutoff, ewaldAccuracy, RealUnits.CoulombConstant), topology));
}

/// <summary>Builds a periodic box of water molecules at a target density.</summary>
public static class WaterBox
{
    /// <summary>
    /// n³ molecules with O atoms on a simple cubic grid, random orientations, equilibrium geometry.
    /// The cubic box side gives the requested mass density (g/cm³).
    /// </summary>
    public static (Box Box, Vec3[] Positions, Topology Topology) Cubic(int perSide, double densityGramsPerCm3, int seed)
    {
        if (perSide < 1) throw new ArgumentOutOfRangeException(nameof(perSide));
        if (!(densityGramsPerCm3 > 0)) throw new ArgumentOutOfRangeException(nameof(densityGramsPerCm3));
        int n = perSide * perSide * perSide;
        double massDensity = RealUnits.GramsPerCubicCentimetreToDensity(densityGramsPerCm3); // (g/mol)/Å³
        double side = Math.Cbrt(n * SpcFw.MolarMass / massDensity);
        var box = new Box(new Vec3(side, side, side));
        double spacing = side / perSide;
        var rng = new Random(seed);
        var topology = new Topology();
        var positions = new Vec3[3 * n];
        int m = 0;
        for (int ix = 0; ix < perSide; ix++)
            for (int iy = 0; iy < perSide; iy++)
                for (int iz = 0; iz < perSide; iz++)
                {
                    topology.AddMolecule(SpcFw.Molecule);
                    var o = spacing * new Vec3(ix + 0.5, iy + 0.5, iz + 0.5);
                    Vec3 u = RandomUnit(rng);
                    Vec3 w;
                    do w = Vec3.Cross(u, RandomUnit(rng)); while (w.Length < 0.1);
                    w /= w.Length;
                    var g = SpcFw.EquilibriumGeometry(o, u, w);
                    for (int k = 0; k < 3; k++) positions[3 * m + k] = box.Wrap(g[k], Int3.Zero).Position;
                    m++;
                }
        return (box, positions, topology);
    }

    private static Vec3 RandomUnit(Random rng)
    {
        while (true)
        {
            var v = new Vec3(2 * rng.NextDouble() - 1, 2 * rng.NextDouble() - 1, 2 * rng.NextDouble() - 1);
            double l = v.Length;
            if (l > 0.1 && l <= 1) return v / l;
        }
    }
}
