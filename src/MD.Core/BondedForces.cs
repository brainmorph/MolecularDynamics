using System;

namespace MD.Core;

/// <summary>
/// Harmonic bond and angle forces from a <see cref="Topology"/> (WATER_PLAN.md D3, W3).
///
///   Bond  U = ½K(r − r0)²            F_i = −K(r − r0)·d/r,  d = r_i − r_j,  F_j = −F_i
///   Angle U = ½K(θ − θ0)²,  θ = ∠(a, b),  a = r_i − r_j,  b = r_k − r_j  (j at the vertex)
///         F_i = K(θ − θ0)/(|a| sin θ)·(b̂ − cos θ·â),  F_k likewise with a and b swapped,
///         F_j = −(F_i + F_k)
///
/// Separation vectors use the minimum image, so a molecule must be smaller than half the box.
/// The virial is Σ d·F over each term; for an angle it is zero (its forces are perpendicular to
/// the arms). An angle at exactly 0 or π has no defined gradient and is reported as an error.
/// </summary>
public sealed class BondedForces : IForceField
{
    public BondedForces(Topology topology)
    {
        Topology = topology ?? throw new ArgumentNullException(nameof(topology));
        foreach (var b in topology.Bonds)
        {
            if (b.Parameters is null) throw new ArgumentException($"Bond {b.I}–{b.J} has no parameters.", nameof(topology));
        }
        foreach (var a in topology.Angles)
        {
            if (a.Parameters is null) throw new ArgumentException($"Angle {a.I}–{a.J}–{a.K} has no parameters.", nameof(topology));
        }
    }

    public Topology Topology { get; }

    public ForceResult Compute(ReadOnlySpan<Vec3> positions, Box box, Span<Vec3> forces)
    {
        if (positions.Length != Topology.AtomCount) throw new ArgumentException("Topology and positions have different atom counts.");
        if (forces.Length != positions.Length) throw new ArgumentException("Force and position arrays must have the same length.");
        forces.Clear();
        double energy = 0;
        double virial = 0;

        foreach (var bond in Topology.Bonds)
        {
            HarmonicBond p = bond.Parameters!;
            Vec3 d = box.MinimumImage(positions[bond.I] - positions[bond.J]);
            double r = d.Length;
            double dr = r - p.Length;
            energy += 0.5 * p.ForceConstant * dr * dr;
            Vec3 fi = (-p.ForceConstant * dr / r) * d;
            forces[bond.I] += fi;
            forces[bond.J] -= fi;
            virial += Vec3.Dot(d, fi);
        }

        foreach (var angle in Topology.Angles)
        {
            HarmonicAngle p = angle.Parameters!;
            Vec3 a = box.MinimumImage(positions[angle.I] - positions[angle.J]);
            Vec3 b = box.MinimumImage(positions[angle.K] - positions[angle.J]);
            double la = a.Length;
            double lb = b.Length;
            Vec3 ua = a / la;
            Vec3 ub = b / lb;
            double cos = Math.Clamp(Vec3.Dot(ua, ub), -1.0, 1.0);
            double sin = Vec3.Cross(ua, ub).Length;
            if (!(sin > 0))
            {
                throw new InvalidOperationException($"Angle {angle.I}–{angle.J}–{angle.K} is linear; its gradient is undefined.");
            }
            double theta = Math.Atan2(sin, cos);
            double dTheta = theta - p.Angle;
            energy += 0.5 * p.ForceConstant * dTheta * dTheta;
            double c = p.ForceConstant * dTheta / sin;
            Vec3 fi = (c / la) * (ub - cos * ua);
            Vec3 fk = (c / lb) * (ua - cos * ub);
            forces[angle.I] += fi;
            forces[angle.K] += fk;
            forces[angle.J] -= fi + fk;
            virial += Vec3.Dot(a, fi) + Vec3.Dot(b, fk);
        }

        return new ForceResult(energy, virial);
    }
}
