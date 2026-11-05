using System;

namespace MD.Core;

/// <summary>
/// Orthorhombic periodic box with corners at (0,0,0) and (Lx,Ly,Lz).
/// Atoms that leave one side re-enter from the opposite side; image counts record how many
/// times each atom has crossed, so unwrapped positions can be recovered exactly.
/// </summary>
public sealed class Box
{
    public Box(Vec3 lengths)
    {
        if (!(lengths.X > 0 && lengths.Y > 0 && lengths.Z > 0) ||
            double.IsInfinity(lengths.X) || double.IsInfinity(lengths.Y) || double.IsInfinity(lengths.Z))
        {
            throw new ArgumentException("Box lengths must be positive and finite.", nameof(lengths));
        }
        Lengths = lengths;
    }

    public Vec3 Lengths { get; }
    public double Volume => Lengths.X * Lengths.Y * Lengths.Z;

    /// <summary>Shortest of the three box lengths.</summary>
    public double MinLength => Math.Min(Lengths.X, Math.Min(Lengths.Y, Lengths.Z));

    /// <summary>Applies the minimum-image convention to a separation vector.</summary>
    public Vec3 MinimumImage(Vec3 d) => new(
        d.X - Lengths.X * Math.Round(d.X / Lengths.X),
        d.Y - Lengths.Y * Math.Round(d.Y / Lengths.Y),
        d.Z - Lengths.Z * Math.Round(d.Z / Lengths.Z));

    /// <summary>Wraps a position into [0, L) on each axis and updates its image counts.</summary>
    public (Vec3 Position, Int3 Image) Wrap(Vec3 position, Int3 image)
    {
        var (x, ix) = Wrap1(position.X, Lengths.X, image.X);
        var (y, iy) = Wrap1(position.Y, Lengths.Y, image.Y);
        var (z, iz) = Wrap1(position.Z, Lengths.Z, image.Z);
        return (new Vec3(x, y, z), new Int3(ix, iy, iz));
    }

    /// <summary>Unwrapped position: wrapped position plus image count times box length.</summary>
    public Vec3 Unwrap(Vec3 position, Int3 image) => new(
        position.X + image.X * Lengths.X,
        position.Y + image.Y * Lengths.Y,
        position.Z + image.Z * Lengths.Z);

    private static (double, int) Wrap1(double x, double length, int image)
    {
        if (x >= 0 && x < length)
        {
            return (x, image);
        }
        double n = Math.Floor(x / length);
        x -= n * length;
        int shift = checked((int)n);
        // Rounding can leave x == length (e.g. x = -1e-17 maps to length); fold it back to 0.
        if (x >= length)
        {
            x -= length;
            shift++;
        }
        if (x < 0)
        {
            x = 0;
        }
        return (x, checked(image + shift));
    }
}
