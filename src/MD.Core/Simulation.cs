using System;

namespace MD.Core;

/// <summary>Snapshot of energies and state variables at one step, all in reduced LJ units.</summary>
public readonly record struct ThermoState(
    long Step,
    double Time,
    double KineticEnergy,
    double PotentialEnergy,
    double TotalEnergy,
    double Temperature,
    double Pressure);

/// <summary>NVE simulation: a particle system advanced by velocity Verlet.</summary>
public sealed class Simulation
{
    private ForceResult _lastForces;

    public Simulation(ParticleSystem system, VelocityVerlet integrator)
    {
        System = system ?? throw new ArgumentNullException(nameof(system));
        Integrator = integrator ?? throw new ArgumentNullException(nameof(integrator));
        _lastForces = Integrator.ComputeForces(System);
        Initial = State;
    }

    public ParticleSystem System { get; }
    public VelocityVerlet Integrator { get; }
    public long StepCount { get; private set; }
    public double Time => StepCount * Integrator.Timestep;

    /// <summary>State at construction (step 0). Its total energy is the reference for drift.</summary>
    public ThermoState Initial { get; }

    public ThermoState State
    {
        get
        {
            double ke = Thermo.KineticEnergy(System);
            double pe = _lastForces.PotentialEnergy;
            return new ThermoState(
                StepCount,
                Time,
                ke,
                pe,
                ke + pe,
                Thermo.Temperature(ke, System.Count),
                Thermo.Pressure(ke, _lastForces.Virial, System.Count, System.Box.Volume));
        }
    }

    public void Step()
    {
        _lastForces = Integrator.Step(System);
        StepCount++;
    }

    /// <summary>Reverses all velocities (used by the time-reversibility test).</summary>
    public void ReverseVelocities()
    {
        for (int i = 0; i < System.Count; i++) System.Velocities[i] = -System.Velocities[i];
    }
}
