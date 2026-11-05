using System;

namespace MD.Core;

/// <summary>Instantaneous thermodynamic quantities, in reduced LJ units (k_B = 1).</summary>
public static class Thermo
{
    public static double KineticEnergy(ParticleSystem s)
    {
        double ke = 0;
        for (int i = 0; i < s.Count; i++) ke += 0.5 * s.Masses[i] * s.Velocities[i].LengthSquared;
        return ke;
    }

    /// <summary>
    /// Degrees of freedom: 3N − 3, because total momentum is removed at the start and conserved after.
    /// </summary>
    public static int DegreesOfFreedom(int atomCount) => Math.Max(3 * atomCount - 3, 0);

    /// <summary>T = 2·KE / (3N − 3). Zero for a single atom (no free degrees of freedom).</summary>
    public static double Temperature(double kineticEnergy, int atomCount)
    {
        int dof = DegreesOfFreedom(atomCount);
        return dof == 0 ? 0 : 2.0 * kineticEnergy / dof;
    }

    /// <summary>
    /// Virial pressure P = (N·T + W/3) / V, with T from <see cref="Temperature"/> and W = Σ r_ij·F_ij.
    /// No long-range tail correction is applied (the potential is truncated and shifted).
    /// </summary>
    public static double Pressure(double kineticEnergy, double virial, int atomCount, double volume)
    {
        double t = Temperature(kineticEnergy, atomCount);
        return (atomCount * t + virial / 3.0) / volume;
    }
}
