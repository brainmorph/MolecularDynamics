using System;

namespace MD.Core;

/// <summary>Builds starting configurations on a face-centred cubic lattice.</summary>
public static class Lattice
{
    public const int AtomsPerFccCell = 4;

    /// <summary>
    /// FCC lattice of nx × ny × nz conventional cells (4 atoms each). The lattice constant is
    /// a = (4/ρ)^(1/3), where ρ is <paramref name="density"/>.
    /// By default the box is exactly (nx·a, ny·a, nz·a) and the density is the true number density.
    /// If <paramref name="box"/> is given, the box is that size instead and the lattice is centred in it,
    /// so the atoms form a compact block with empty space around them. In that case
    /// <paramref name="density"/> only sets the lattice spacing, not the system's density.
    /// </summary>
    public static (Box Box, Vec3[] Positions) Fcc(int nx, int ny, int nz, double density, Vec3? box = null)
    {
        if (nx < 1 || ny < 1 || nz < 1) throw new ArgumentOutOfRangeException(nameof(nx), "Cell counts must be at least 1.");
        if (!(density > 0) || double.IsInfinity(density)) throw new ArgumentOutOfRangeException(nameof(density));

        double a = Math.Cbrt(AtomsPerFccCell / density);
        var extent = new Vec3(nx * a, ny * a, nz * a);
        Vec3 lengths = box ?? extent;
        if (lengths.X < extent.X || lengths.Y < extent.Y || lengths.Z < extent.Z)
            throw new ArgumentException($"Box {lengths} is smaller than the lattice {extent}.", nameof(box));
        Vec3 offset = (lengths - extent) / 2;
        var basis = new[]
        {
            new Vec3(0.0, 0.0, 0.0),
            new Vec3(0.5, 0.5, 0.0),
            new Vec3(0.5, 0.0, 0.5),
            new Vec3(0.0, 0.5, 0.5),
        };

        var positions = new Vec3[checked(AtomsPerFccCell * nx * ny * nz)];
        int k = 0;
        for (int ix = 0; ix < nx; ix++)
        for (int iy = 0; iy < ny; iy++)
        for (int iz = 0; iz < nz; iz++)
        foreach (Vec3 b in basis)
        {
            positions[k++] = offset + a * (new Vec3(ix, iy, iz) + b);
        }
        return (new Box(lengths), positions);
    }

    /// <summary>
    /// <paramref name="count"/> identical FCC blocks (nx × ny × nz cells each) in a straight row along x,
    /// separated by <paramref name="gap"/> σ between the cell boundaries of neighbouring blocks (the outermost
    /// atom planes sit a/2 inside those boundaries). The row is centred
    /// in <paramref name="box"/>, and each block is centred on y and z. The box must be large enough
    /// to hold the whole row.
    /// </summary>
    public static (Box Box, Vec3[] Positions) Chain(int nx, int ny, int nz, double density, int count, double gap, Vec3 box)
    {
        if (nx < 1 || ny < 1 || nz < 1) throw new ArgumentOutOfRangeException(nameof(nx), "Cell counts must be at least 1.");
        if (count < 1) throw new ArgumentOutOfRangeException(nameof(count), "Cluster count must be at least 1.");
        if (!(gap >= 0) || double.IsInfinity(gap)) throw new ArgumentOutOfRangeException(nameof(gap), "Gap must be zero or positive.");
        if (!(density > 0) || double.IsInfinity(density)) throw new ArgumentOutOfRangeException(nameof(density));

        double a = Math.Cbrt(AtomsPerFccCell / density);
        var block = new Vec3(nx * a, ny * a, nz * a);
        double rowLength = count * block.X + (count - 1) * gap;
        if (rowLength > box.X || block.Y > box.Y || block.Z > box.Z)
            throw new ArgumentException($"Box {box} cannot hold {count} blocks of {block} with gap {gap} (row length {rowLength}).", nameof(box));

        var basis = new[]
        {
            new Vec3(0.0, 0.0, 0.0),
            new Vec3(0.5, 0.5, 0.0),
            new Vec3(0.5, 0.0, 0.5),
            new Vec3(0.0, 0.5, 0.5),
        };
        var positions = new Vec3[checked(AtomsPerFccCell * nx * ny * nz * count)];
        double startX = (box.X - rowLength) / 2;
        var yz = new Vec3(0, (box.Y - block.Y) / 2, (box.Z - block.Z) / 2);
        int k = 0;
        for (int c = 0; c < count; c++)
        {
            var origin = new Vec3(startX + c * (block.X + gap), yz.Y, yz.Z);
            for (int ix = 0; ix < nx; ix++)
            for (int iy = 0; iy < ny; iy++)
            for (int iz = 0; iz < nz; iz++)
            foreach (Vec3 b in basis)
            {
                positions[k++] = origin + a * (new Vec3(ix, iy, iz) + b);
            }
        }
        return (new Box(box), positions);
    }
}
