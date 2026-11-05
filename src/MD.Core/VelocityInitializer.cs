using System;

namespace MD.Core;

/// <summary>Assigns random starting velocities for a target temperature.</summary>
public static class VelocityInitializer
{
    /// <summary>
    /// Draws each velocity component from a Gaussian with variance T/m (Maxwell–Boltzmann), removes
    /// the centre-of-mass velocity so the box does not drift, then rescales so the instantaneous
    /// temperature equals <paramref name="temperature"/> exactly. Deterministic for a given seed.
    /// </summary>
    public static void Assign(ParticleSystem s, double temperature, int seed)
    {
        if (!(temperature >= 0) || double.IsInfinity(temperature)) throw new ArgumentOutOfRangeException(nameof(temperature));

        var rng = new Random(seed);
        for (int i = 0; i < s.Count; i++)
        {
            double sd = Math.Sqrt(temperature / s.Masses[i]);
            s.Velocities[i] = sd * new Vec3(Gaussian(rng), Gaussian(rng), Gaussian(rng));
        }
        RemoveNetMomentum(s);

        double ke = Thermo.KineticEnergy(s);
        double current = Thermo.Temperature(ke, s.Count);
        double scale = current > 0 ? Math.Sqrt(temperature / current) : 0;
        for (int i = 0; i < s.Count; i++) s.Velocities[i] = scale * s.Velocities[i];
    }

    public static void RemoveNetMomentum(ParticleSystem s)
    {
        Vec3 vcm = s.TotalMomentum() / s.TotalMass();
        for (int i = 0; i < s.Count; i++) s.Velocities[i] -= vcm;
    }

    /// <summary>Standard normal sample by the Box–Muller transform.</summary>
    private static double Gaussian(Random rng)
    {
        double u1 = 1.0 - rng.NextDouble(); // in (0, 1], avoids log(0)
        double u2 = rng.NextDouble();
        return Math.Sqrt(-2.0 * Math.Log(u1)) * Math.Cos(2.0 * Math.PI * u2);
    }
}
