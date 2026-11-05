using System;

namespace MD.Core;

/// <summary>Double-precision 3-vector. The BCL has no double-precision Vector3, so this is hand-written.</summary>
public readonly struct Vec3 : IEquatable<Vec3>
{
    public readonly double X;
    public readonly double Y;
    public readonly double Z;

    public Vec3(double x, double y, double z)
    {
        X = x;
        Y = y;
        Z = z;
    }

    public static Vec3 Zero => new(0, 0, 0);

    public double LengthSquared => X * X + Y * Y + Z * Z;
    public double Length => Math.Sqrt(LengthSquared);

    public double this[int axis] => axis switch
    {
        0 => X,
        1 => Y,
        2 => Z,
        _ => throw new ArgumentOutOfRangeException(nameof(axis)),
    };

    public static double Dot(Vec3 a, Vec3 b) => a.X * b.X + a.Y * b.Y + a.Z * b.Z;

    public static Vec3 operator +(Vec3 a, Vec3 b) => new(a.X + b.X, a.Y + b.Y, a.Z + b.Z);
    public static Vec3 operator -(Vec3 a, Vec3 b) => new(a.X - b.X, a.Y - b.Y, a.Z - b.Z);
    public static Vec3 operator -(Vec3 a) => new(-a.X, -a.Y, -a.Z);
    public static Vec3 operator *(double s, Vec3 a) => new(s * a.X, s * a.Y, s * a.Z);
    public static Vec3 operator *(Vec3 a, double s) => new(s * a.X, s * a.Y, s * a.Z);
    public static Vec3 operator /(Vec3 a, double s) => new(a.X / s, a.Y / s, a.Z / s);

    public bool Equals(Vec3 other) => X == other.X && Y == other.Y && Z == other.Z;
    public override bool Equals(object? obj) => obj is Vec3 v && Equals(v);
    public override int GetHashCode() => HashCode.Combine(X, Y, Z);
    public static bool operator ==(Vec3 a, Vec3 b) => a.Equals(b);
    public static bool operator !=(Vec3 a, Vec3 b) => !a.Equals(b);

    public override string ToString() => $"({X:G17}, {Y:G17}, {Z:G17})";
}

/// <summary>Integer 3-vector, used for periodic image counts.</summary>
public readonly record struct Int3(int X, int Y, int Z)
{
    public static Int3 Zero => new(0, 0, 0);
}
