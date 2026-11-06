using System;
using System.Collections.Generic;

namespace MD.Core;

/// <summary>
/// An atom type for molecular simulations, in real units (<see cref="RealUnits"/>).
/// A type with <see cref="LjEpsilon"/> = 0 has no LJ interaction (charge only, like H in SPC/Fw).
/// </summary>
public sealed record AtomType
{
    public AtomType(string name, string element, double mass, double charge, double ljSigma = 0, double ljEpsilon = 0)
    {
        if (string.IsNullOrWhiteSpace(name)) throw new ArgumentException("Atom type needs a name.", nameof(name));
        if (string.IsNullOrWhiteSpace(element)) throw new ArgumentException("Atom type needs an element.", nameof(element));
        if (!(mass > 0) || double.IsInfinity(mass)) throw new ArgumentOutOfRangeException(nameof(mass));
        if (!double.IsFinite(charge)) throw new ArgumentOutOfRangeException(nameof(charge));
        if (!(ljSigma >= 0) || double.IsInfinity(ljSigma)) throw new ArgumentOutOfRangeException(nameof(ljSigma));
        if (!(ljEpsilon >= 0) || double.IsInfinity(ljEpsilon)) throw new ArgumentOutOfRangeException(nameof(ljEpsilon));
        if ((ljEpsilon > 0) != (ljSigma > 0)) throw new ArgumentException("LJ sigma and epsilon must both be set or both be zero.");
        Name = name;
        Element = element;
        Mass = mass;
        Charge = charge;
        LjSigma = ljSigma;
        LjEpsilon = ljEpsilon;
    }

    public string Name { get; }
    public string Element { get; }
    /// <summary>g/mol.</summary>
    public double Mass { get; }
    /// <summary>e.</summary>
    public double Charge { get; }
    /// <summary>Å; 0 when the type has no LJ site.</summary>
    public double LjSigma { get; }
    /// <summary>kJ/mol; 0 when the type has no LJ site.</summary>
    public double LjEpsilon { get; }
    public bool HasLj => LjEpsilon > 0;
}

/// <summary>
/// Harmonic bond U = ½·K·(r − R0)², with K in kJ/(mol·Å²) and R0 in Å. Note the ½: a parameter
/// quoted without it (U = K(r − R0)²) must be doubled.
/// </summary>
public sealed record HarmonicBond
{
    public HarmonicBond(double forceConstant, double length)
    {
        if (!(forceConstant > 0) || double.IsInfinity(forceConstant)) throw new ArgumentOutOfRangeException(nameof(forceConstant));
        if (!(length > 0) || double.IsInfinity(length)) throw new ArgumentOutOfRangeException(nameof(length));
        ForceConstant = forceConstant;
        Length = length;
    }

    public double ForceConstant { get; }
    public double Length { get; }
}

/// <summary>
/// Harmonic angle U = ½·K·(θ − θ0)², with K in kJ/(mol·rad²) and θ0 in radians. Same ½ convention
/// as <see cref="HarmonicBond"/>.
/// </summary>
public sealed record HarmonicAngle
{
    public HarmonicAngle(double forceConstant, double angle)
    {
        if (!(forceConstant > 0) || double.IsInfinity(forceConstant)) throw new ArgumentOutOfRangeException(nameof(forceConstant));
        if (!(angle > 0 && angle < Math.PI)) throw new ArgumentOutOfRangeException(nameof(angle), "Equilibrium angle must be strictly between 0 and π.");
        ForceConstant = forceConstant;
        Angle = angle;
    }

    public double ForceConstant { get; }
    public double Angle { get; }
}

/// <summary>
/// Bond between atoms <see cref="I"/> and <see cref="J"/> (indices local to a template, or global in a
/// topology). <see cref="Parameters"/> may be null for topology-only use; <see cref="BondedForces"/>
/// refuses a bond without parameters.
/// </summary>
public readonly record struct BondTerm(int I, int J, HarmonicBond? Parameters = null);

/// <summary>Angle I–J–K with <see cref="J"/> at the vertex. Parameters as for <see cref="BondTerm"/>.</summary>
public readonly record struct AngleTerm(int I, int J, int K, HarmonicAngle? Parameters = null);

/// <summary>
/// A molecule definition: its atom types, bonds and angles, with atom indices local to the molecule.
/// Bonded parameters are attached to each <see cref="BondTerm"/> and <see cref="AngleTerm"/>.
/// </summary>
public sealed class MoleculeTemplate
{
    public MoleculeTemplate(string name, IReadOnlyList<AtomType> atoms, IReadOnlyList<BondTerm> bonds, IReadOnlyList<AngleTerm> angles)
    {
        if (string.IsNullOrWhiteSpace(name)) throw new ArgumentException("Molecule needs a name.", nameof(name));
        ArgumentNullException.ThrowIfNull(atoms);
        ArgumentNullException.ThrowIfNull(bonds);
        ArgumentNullException.ThrowIfNull(angles);
        if (atoms.Count == 0) throw new ArgumentException("Molecule needs at least one atom.", nameof(atoms));
        int n = atoms.Count;
        foreach (var b in bonds)
        {
            CheckIndex(b.I, n);
            CheckIndex(b.J, n);
            if (b.I == b.J) throw new ArgumentException($"BondTerm {b} joins an atom to itself.");
        }
        foreach (var a in angles)
        {
            CheckIndex(a.I, n);
            CheckIndex(a.J, n);
            CheckIndex(a.K, n);
            if (a.I == a.J || a.J == a.K || a.I == a.K) throw new ArgumentException($"AngleTerm {a} repeats an atom.");
        }
        Name = name;
        Atoms = new List<AtomType>(atoms).AsReadOnly();
        Bonds = new List<BondTerm>(bonds).AsReadOnly();
        Angles = new List<AngleTerm>(angles).AsReadOnly();
    }

    public string Name { get; }
    public IReadOnlyList<AtomType> Atoms { get; }
    public IReadOnlyList<BondTerm> Bonds { get; }
    public IReadOnlyList<AngleTerm> Angles { get; }

    public double TotalCharge()
    {
        double q = 0;
        foreach (var a in Atoms) q += a.Charge;
        return q;
    }

    private static void CheckIndex(int i, int n)
    {
        if (i < 0 || i >= n) throw new ArgumentOutOfRangeException(nameof(i), $"Atom index {i} is outside 0..{n - 1}.");
    }
}

/// <summary>
/// Per-atom types and charges, and per-molecule bonds and angles, for a whole system.
/// Atom order is molecule by molecule, in the order molecules were added.
///
/// Exclusion rule: two atoms in the same molecule have no LJ or Coulomb interaction; they interact
/// only through bonded terms. For a three-site water this is the same as excluding 1-2 and 1-3
/// pairs. It must be checked against the SPC/Fw paper before W4 (WATER_PLAN.md, W1).
/// </summary>
public sealed class Topology
{
    private readonly List<AtomType> _types = new();
    private readonly List<int> _moleculeOf = new();
    private readonly List<int> _moleculeStart = new();
    private readonly List<MoleculeTemplate> _moleculeTemplate = new();
    private readonly List<BondTerm> _bonds = new();
    private readonly List<AngleTerm> _angles = new();

    public int AtomCount => _types.Count;
    public int MoleculeCount => _moleculeStart.Count;
    public IReadOnlyList<BondTerm> Bonds => _bonds;
    public IReadOnlyList<AngleTerm> Angles => _angles;

    public AtomType TypeOf(int atom) => _types[atom];
    public int MoleculeOf(int atom) => _moleculeOf[atom];
    public int FirstAtomOf(int molecule) => _moleculeStart[molecule];
    public MoleculeTemplate TemplateOf(int molecule) => _moleculeTemplate[molecule];

    /// <summary>Appends one molecule and returns its index.</summary>
    public int AddMolecule(MoleculeTemplate template)
    {
        ArgumentNullException.ThrowIfNull(template);
        int molecule = _moleculeStart.Count;
        int offset = _types.Count;
        _moleculeStart.Add(offset);
        _moleculeTemplate.Add(template);
        foreach (var t in template.Atoms)
        {
            _types.Add(t);
            _moleculeOf.Add(molecule);
        }
        foreach (var b in template.Bonds) _bonds.Add(b with { I = b.I + offset, J = b.J + offset });
        foreach (var a in template.Angles) _angles.Add(a with { I = a.I + offset, J = a.J + offset, K = a.K + offset });
        return molecule;
    }

    /// <summary>True when atoms i and j (i ≠ j) have no non-bonded (LJ or Coulomb) interaction.</summary>
    public bool IsExcluded(int i, int j) => _moleculeOf[i] == _moleculeOf[j];

    public double[] Masses() => Map(t => t.Mass);
    public double[] Charges() => Map(t => t.Charge);

    public string[] Elements()
    {
        var e = new string[_types.Count];
        for (int i = 0; i < e.Length; i++) e[i] = _types[i].Element;
        return e;
    }

    public double TotalCharge()
    {
        double q = 0;
        foreach (var t in _types) q += t.Charge;
        return q;
    }

    /// <summary>Every excluded pair (i &lt; j) once: the atom pairs within each molecule.</summary>
    public IEnumerable<(int I, int J)> ExcludedPairs()
    {
        for (int m = 0; m < _moleculeStart.Count; m++)
        {
            int first = _moleculeStart[m];
            int end = first + _moleculeTemplate[m].Atoms.Count;
            for (int i = first; i < end; i++)
                for (int j = i + 1; j < end; j++)
                    yield return (i, j);
        }
    }

    private double[] Map(Func<AtomType, double> f)
    {
        var v = new double[_types.Count];
        for (int i = 0; i < v.Length; i++) v[i] = f(_types[i]);
        return v;
    }
}
