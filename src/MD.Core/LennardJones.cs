using System;

namespace MD.Core;

/// <summary>
/// Lennard-Jones 12-6 pair potential, truncated at a cutoff and shifted so the energy is zero there:
///   U(r) = 4ε[(σ/r)^12 − (σ/r)^6] − U_LJ(r_c)   for r &lt; r_c, and 0 otherwise.
/// The force is the unshifted LJ force (the shift is a constant, so it adds no force).
/// Note: the force itself still jumps at r_c; only the energy is continuous.
/// </summary>
public sealed class LennardJones
{
    public LennardJones(double sigma = 1.0, double epsilon = 1.0, double cutoff = 2.5)
    {
        if (!(sigma > 0)) throw new ArgumentOutOfRangeException(nameof(sigma));
        if (!(epsilon > 0)) throw new ArgumentOutOfRangeException(nameof(epsilon));
        if (!(cutoff > 0) || double.IsInfinity(cutoff)) throw new ArgumentOutOfRangeException(nameof(cutoff));
        Sigma = sigma;
        Epsilon = epsilon;
        Cutoff = cutoff;
        CutoffSquared = cutoff * cutoff;
        EnergyShift = UnshiftedEnergy(cutoff * cutoff);
    }

    public double Sigma { get; }
    public double Epsilon { get; }
    public double Cutoff { get; }
    public double CutoffSquared { get; }

    /// <summary>U_LJ(r_c), subtracted from every pair energy inside the cutoff.</summary>
    public double EnergyShift { get; }

    /// <summary>Plain LJ energy with no cutoff or shift, for a squared distance r².</summary>
    public double UnshiftedEnergy(double r2)
    {
        double s2 = Sigma * Sigma / r2;
        double s6 = s2 * s2 * s2;
        return 4.0 * Epsilon * (s6 * s6 - s6);
    }

    /// <summary>
    /// Evaluates the truncated, shifted potential at squared distance r².
    /// <paramref name="forceOverR"/> is f such that the force on atom i from atom j is f·(r_i − r_j).
    /// Returns false (with zero outputs) when r ≥ r_c.
    /// </summary>
    public bool Evaluate(double r2, out double energy, out double forceOverR)
    {
        if (r2 >= CutoffSquared)
        {
            energy = 0;
            forceOverR = 0;
            return false;
        }
        double s2 = Sigma * Sigma / r2;
        double s6 = s2 * s2 * s2;
        double s12 = s6 * s6;
        energy = 4.0 * Epsilon * (s12 - s6) - EnergyShift;
        forceOverR = 24.0 * Epsilon * (2.0 * s12 - s6) / r2;
        return true;
    }
}
