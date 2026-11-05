using System;

namespace MD.Core;

/// <summary>
/// Velocity Verlet integrator (NVE). Time-reversible and symplectic.
///   v(t+½dt) = v(t) + ½dt·F(t)/m
///   x(t+dt)  = x(t) + dt·v(t+½dt)          (then wrapped into the box)
///   F(t+dt)  = F(x(t+dt))
///   v(t+dt)  = v(t+½dt) + ½dt·F(t+dt)/m
/// </summary>
public sealed class VelocityVerlet
{
    public VelocityVerlet(ForceCalculator forces, double timestep)
    {
        Forces = forces ?? throw new ArgumentNullException(nameof(forces));
        if (!(timestep > 0) || double.IsInfinity(timestep)) throw new ArgumentOutOfRangeException(nameof(timestep));
        Timestep = timestep;
    }

    public ForceCalculator Forces { get; }
    public double Timestep { get; }

    /// <summary>Computes forces for the current positions. Call once before the first <see cref="Step"/>.</summary>
    public ForceResult ComputeForces(ParticleSystem s) => Forces.Compute(s.Positions, s.Box, s.Forces);

    /// <summary>Advances one step. Assumes <see cref="ParticleSystem.Forces"/> hold forces for the current positions.</summary>
    public ForceResult Step(ParticleSystem s)
    {
        double dt = Timestep;
        double halfDt = 0.5 * dt;
        for (int i = 0; i < s.Count; i++)
        {
            Vec3 v = s.Velocities[i] + (halfDt / s.Masses[i]) * s.Forces[i];
            s.Velocities[i] = v;
            (s.Positions[i], s.Images[i]) = s.Box.Wrap(s.Positions[i] + dt * v, s.Images[i]);
        }

        ForceResult result = ComputeForces(s);

        for (int i = 0; i < s.Count; i++)
        {
            s.Velocities[i] += (halfDt / s.Masses[i]) * s.Forces[i];
        }
        return result;
    }
}
