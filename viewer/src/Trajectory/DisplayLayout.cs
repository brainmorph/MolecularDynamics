// Display placement of atoms. It uses no Godot types and no MD.Core, so the core test suite can
// compile this file and check it against the core's output.
//
// Display arithmetic only: it moves atoms by whole box lengths to choose which periodic copy to
// draw, and interpolates between two stored frames. It never computes forces, energies or bonds;
// bonds come from the file.
using System;
using System.Collections.Generic;

namespace MDViewer.Trajectory;

public sealed class DisplayLayout
{
    private uint[]? _bondsFor;
    private int _atomsFor = -1;
    private int[] _adjacencyStart = Array.Empty<int>();
    private int[] _adjacency = Array.Empty<int>();
    private bool[] _placed = Array.Empty<bool>();
    private double[] _unwrapped = Array.Empty<double>();
    private readonly Stack<int> _stack = new();

    /// <summary>
    /// Fills <paramref name="output"/> (x0 y0 z0 x1 ...) with display positions between frames
    /// <paramref name="a"/> and <paramref name="b"/> (<paramref name="alpha"/> from 0 to 1).
    ///
    /// Unwrapped positions are interpolated, so an atom crossing a boundary does not sweep across the
    /// box. Atoms joined by bonds (frame a's bond list) are drawn as one piece: the lowest-index atom
    /// of each bonded group is wrapped into the box, and each other atom is drawn at the periodic copy
    /// nearest the atom it is bonded to. An atom with no bonds is simply wrapped into the box.
    /// </summary>
    public void Positions(TrajectoryFrame a, TrajectoryFrame b, double alpha, int atomCount, double[] output)
    {
        if (output.Length < 3 * atomCount) throw new ArgumentException("Output array is too short.", nameof(output));
        BuildAdjacency(a.Bonds, atomCount);
        if (_unwrapped.Length < 3 * atomCount) _unwrapped = new double[3 * atomCount];
        double[] box = { a.BoxX, a.BoxY, a.BoxZ };

        for (int i = 0; i < atomCount; i++)
            for (int k = 0; k < 3; k++)
            {
                double x = a.Unwrapped(i, k);
                _unwrapped[3 * i + k] = x + (b.Unwrapped(i, k) - x) * alpha;
            }

        Array.Clear(_placed, 0, atomCount);
        for (int root = 0; root < atomCount; root++)
        {
            if (_placed[root]) continue;
            for (int k = 0; k < 3; k++) output[3 * root + k] = Wrap(_unwrapped[3 * root + k], box[k]);
            _placed[root] = true;
            _stack.Push(root);
            while (_stack.Count > 0)
            {
                int i = _stack.Pop();
                for (int e = _adjacencyStart[i]; e < _adjacencyStart[i + 1]; e++)
                {
                    int j = _adjacency[e];
                    if (_placed[j]) continue;
                    for (int k = 0; k < 3; k++)
                    {
                        double d = _unwrapped[3 * j + k] - _unwrapped[3 * i + k];
                        output[3 * j + k] = output[3 * i + k] + d - box[k] * Math.Round(d / box[k]);
                    }
                    _placed[j] = true;
                    _stack.Push(j);
                }
            }
        }
    }

    /// <summary>Wraps x into [0, length).</summary>
    public static double Wrap(double x, double length)
    {
        double w = x - length * Math.Floor(x / length);
        return w >= length ? 0 : w;
    }

    private void BuildAdjacency(uint[] bonds, int atomCount)
    {
        if (ReferenceEquals(bonds, _bondsFor) && atomCount == _atomsFor) return;
        // Frames usually repeat the same bonds; reuse the adjacency when the content is unchanged.
        if (_bondsFor != null && atomCount == _atomsFor && bonds.AsSpan().SequenceEqual(_bondsFor))
        {
            _bondsFor = bonds;
            return;
        }
        _bondsFor = bonds;
        _atomsFor = atomCount;
        _placed = new bool[atomCount];
        var degree = new int[atomCount + 1];
        for (int e = 0; e < bonds.Length; e++) degree[bonds[e]]++;
        _adjacencyStart = new int[atomCount + 1];
        for (int i = 0; i < atomCount; i++) _adjacencyStart[i + 1] = _adjacencyStart[i] + degree[i];
        _adjacency = new int[bonds.Length];
        var fill = (int[])_adjacencyStart.Clone();
        for (int e = 0; e + 1 < bonds.Length; e += 2)
        {
            int i = (int)bonds[e], j = (int)bonds[e + 1];
            _adjacency[fill[i]++] = j;
            _adjacency[fill[j]++] = i;
        }
    }
}
