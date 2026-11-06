using System;

namespace MD.Core;

/// <summary>A source of forces: writes forces for the given positions and returns energy and virial.</summary>
public interface IForceField
{
    /// <summary>Computes forces into <paramref name="forces"/> (overwritten) and returns energy and virial.</summary>
    ForceResult Compute(ReadOnlySpan<Vec3> positions, Box box, Span<Vec3> forces);
}

/// <summary>Sum of several force fields (for example LJ + Coulomb + bonded terms).</summary>
public sealed class CompositeForceField : IForceField
{
    private readonly IForceField[] _terms;
    private Vec3[] _scratch = Array.Empty<Vec3>();

    public CompositeForceField(params IForceField[] terms)
    {
        ArgumentNullException.ThrowIfNull(terms);
        if (terms.Length == 0) throw new ArgumentException("At least one force field is needed.", nameof(terms));
        foreach (var t in terms) ArgumentNullException.ThrowIfNull(t, nameof(terms));
        _terms = (IForceField[])terms.Clone();
    }

    public ForceResult Compute(ReadOnlySpan<Vec3> positions, Box box, Span<Vec3> forces)
    {
        if (_scratch.Length != forces.Length) _scratch = new Vec3[forces.Length];
        ForceResult total = _terms[0].Compute(positions, box, forces);
        double energy = total.PotentialEnergy;
        double virial = total.Virial;
        for (int t = 1; t < _terms.Length; t++)
        {
            ForceResult r = _terms[t].Compute(positions, box, _scratch);
            energy += r.PotentialEnergy;
            virial += r.Virial;
            for (int i = 0; i < forces.Length; i++) forces[i] += _scratch[i];
        }
        return new ForceResult(energy, virial);
    }
}
