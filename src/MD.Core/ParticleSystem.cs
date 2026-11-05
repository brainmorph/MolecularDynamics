using System;

namespace MD.Core;

/// <summary>
/// State of all particles. Positions are kept wrapped into the box; <see cref="Images"/> records
/// boundary crossings so unwrapped positions are exact. Everything is in reduced LJ units.
/// </summary>
public sealed class ParticleSystem
{
    public ParticleSystem(Box box, Vec3[] positions, double[] masses, string[] elements)
    {
        Box = box ?? throw new ArgumentNullException(nameof(box));
        ArgumentNullException.ThrowIfNull(positions);
        ArgumentNullException.ThrowIfNull(masses);
        ArgumentNullException.ThrowIfNull(elements);
        int n = positions.Length;
        if (masses.Length != n || elements.Length != n)
        {
            throw new ArgumentException("positions, masses and elements must have the same length.");
        }
        foreach (double m in masses)
        {
            if (!(m > 0) || double.IsInfinity(m)) throw new ArgumentException("Masses must be positive and finite.");
        }

        Positions = new Vec3[n];
        Images = new Int3[n];
        for (int i = 0; i < n; i++)
        {
            (Positions[i], Images[i]) = box.Wrap(positions[i], Int3.Zero);
        }
        Velocities = new Vec3[n];
        Forces = new Vec3[n];
        Masses = (double[])masses.Clone();
        Elements = (string[])elements.Clone();
    }

    public Box Box { get; }
    public int Count => Positions.Length;
    public Vec3[] Positions { get; }
    public Int3[] Images { get; }
    public Vec3[] Velocities { get; }
    public Vec3[] Forces { get; }
    public double[] Masses { get; }
    public string[] Elements { get; }

    public Vec3 UnwrappedPosition(int i) => Box.Unwrap(Positions[i], Images[i]);

    public Vec3 TotalMomentum()
    {
        Vec3 p = Vec3.Zero;
        for (int i = 0; i < Count; i++) p += Masses[i] * Velocities[i];
        return p;
    }

    public double TotalMass()
    {
        double m = 0;
        foreach (double mi in Masses) m += mi;
        return m;
    }
}
