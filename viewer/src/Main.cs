using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using Godot;
using MDViewer.Trajectory;

namespace MDViewer;

/// <summary>
/// Plays back a trajectory file. Display only: every energy shown comes from the file, written by
/// the core, and bonds are drawn from the file's bond list. The viewer never computes forces,
/// energies or bonds.
///
/// One Godot unit is one length unit of the file: σ for "lj-reduced" runs, Å for "real" runs.
///
/// Usage: godot --path viewer -- path/to/run.traj   (default: ../out/argon-8.traj)
/// </summary>
public partial class Main : Node3D
{
    /// <summary>Trajectory frames shown per second of playback.</summary>
    private const double PlaybackFramesPerSecond = 30.0;

    /// <summary>Display radius of an atom sphere in "lj-reduced" files, in σ. Smaller than the
    /// contact distance (~0.5σ) so the structure inside the box stays visible.</summary>
    private const float ReducedAtomRadius = 0.3f;

    /// <summary>Bond stick radius, in the file's length unit (σ or Å).</summary>
    private const float BondRadius = 0.09f;

    /// <summary>Ball-and-stick sphere radii for "real" files, in Å (display sizes, not physical radii).</summary>
    private static float RealAtomRadius(string element) => element switch
    {
        "H" => 0.22f,
        "O" => 0.36f,
        _ => 0.32f,
    };

    private const string RealUnits = "real";

    private TrajectoryFile? _trajectory;
    private string _loadedPath = "";
    private double _playhead; // fractional frame index

    private FlyCamera _camera = null!;
    private MultiMesh _atoms = null!;
    private MultiMesh _bonds = null!;
    private readonly DisplayLayout _layout = new();
    private double[] _layoutPositions = Array.Empty<double>();
    private float[] _radii = Array.Empty<float>();
    private Color[] _colors = Array.Empty<Color>();
    private MeshInstance3D _boxOutline = null!;
    private ImmediateMesh _boxMesh = null!;
    private Vector3 _drawnBox = Vector3.Zero;
    private Vector3[] _displayPositions = Array.Empty<Vector3>();

    private PanelContainer _diagnostics = null!;
    private Label _diagnosticsText = null!;
    private Control _menu = null!;
    private Label _message = null!;
    private bool _menuOpen;

    public override void _Ready()
    {
        BuildWorld();
        BuildUi();
        LoadTrajectory(ResolveTrajectoryPath());
        // Render at 60 frames per second. This changes only how often the screen updates, not playback speed.
        Engine.MaxFps = 60;
        // The mouse starts free. It is captured only by an explicit click in the window (see _Input).
    }

    // ---------------------------------------------------------------- input

    public override void _Input(InputEvent @event)
    {
        if (@event is InputEventKey { Pressed: true, Echo: false } key)
        {
            if (key.PhysicalKeycode == Key.Escape)
            {
                SetMenuOpen(!_menuOpen);
                GetViewport().SetInputAsHandled();
            }
            else if (key.PhysicalKeycode == Key.F1)
            {
                // Diagnostics only toggle visibility; they never change camera, playback or data.
                _diagnostics.Visible = !_diagnostics.Visible;
                GetViewport().SetInputAsHandled();
            }
        }
        else if (@event is InputEventMouseButton { Pressed: true } && !_menuOpen &&
                 Input.MouseMode != Input.MouseModeEnum.Captured)
        {
            // Recapture after the window lost focus.
            Input.MouseMode = Input.MouseModeEnum.Captured;
        }
    }

    private void SetMenuOpen(bool open)
    {
        _menuOpen = open;
        _menu.Visible = open;
        _camera.Active = !open;
        Input.MouseMode = open ? Input.MouseModeEnum.Visible : Input.MouseModeEnum.Captured;
    }

    private void Restart()
    {
        _playhead = 0;
        SetMenuOpen(false);
    }

    // ---------------------------------------------------------------- frame update

    public override void _Process(double delta)
    {

        if (_trajectory is { Frames.Count: > 0 } t)
        {
            // Playback pauses while the menu is open. It stops on the last frame.
            if (!_menuOpen)
            {
                _playhead = Math.Min(_playhead + delta * PlaybackFramesPerSecond, t.Frames.Count - 1);
            }
            ShowPlayhead(t);
        }
        if (_diagnostics.Visible) UpdateDiagnostics(delta);
    }

    private void ShowPlayhead(TrajectoryFile t)
    {
        int index = (int)Math.Floor(_playhead);
        TrajectoryFrame a = t.Frames[index];
        TrajectoryFrame b = t.Frames[Math.Min(index + 1, t.Frames.Count - 1)];
        double alpha = _playhead - index;

        var box = new Vector3((float)a.BoxX, (float)a.BoxY, (float)a.BoxZ);
        if (box != _drawnBox) DrawBox(box);

        // Display-only placement: interpolate unwrapped positions, wrap into the box, and keep
        // bonded atoms together (see DisplayLayout).
        _layout.Positions(a, b, alpha, t.AtomCount, _layoutPositions);
        for (int i = 0; i < t.AtomCount; i++)
        {
            var p = new Vector3((float)_layoutPositions[3 * i], (float)_layoutPositions[3 * i + 1], (float)_layoutPositions[3 * i + 2]);
            _displayPositions[i] = p;
            _atoms.SetInstanceTransform(i, new Transform3D(Basis.Identity.Scaled(Vector3.One * _radii[i]), p));
        }

        // Each bond is two half-sticks, each in the colour of its own atom.
        int bondCount = a.Bonds.Length / 2;
        if (_bonds.InstanceCount != 2 * bondCount) _bonds.InstanceCount = 2 * bondCount;
        for (int k = 0; k < bondCount; k++)
        {
            int i = (int)a.Bonds[2 * k], j = (int)a.Bonds[2 * k + 1];
            Vector3 mid = (_displayPositions[i] + _displayPositions[j]) / 2;
            _bonds.SetInstanceTransform(2 * k, Stick(_displayPositions[i], mid));
            _bonds.SetInstanceTransform(2 * k + 1, Stick(mid, _displayPositions[j]));
            _bonds.SetInstanceColor(2 * k, _colors[i]);
            _bonds.SetInstanceColor(2 * k + 1, _colors[j]);
        }
    }

    /// <summary>Transform taking the unit cylinder (radius 1, height 1, along y, centred) to a stick from <paramref name="from"/> to <paramref name="to"/>.</summary>
    private static Transform3D Stick(Vector3 from, Vector3 to)
    {
        Vector3 d = to - from;
        float length = d.Length();
        if (length < 1e-6f) return new Transform3D(Basis.Identity.Scaled(Vector3.Zero), from);
        Vector3 axis = d / length;
        Vector3 helper = Math.Abs(axis.X) < 0.9f ? Vector3.Right : Vector3.Up;
        Vector3 u = axis.Cross(helper).Normalized();
        Vector3 w = axis.Cross(u);
        return new Transform3D(new Basis(u * BondRadius, d, w * BondRadius), from + d / 2);
    }

    // ---------------------------------------------------------------- diagnostics (F1)

    private void UpdateDiagnostics(double delta)
    {
        var ci = CultureInfo.InvariantCulture;
        var text = new System.Text.StringBuilder();
        text.AppendLine(string.Format(ci, "FPS {0:F0}   frame time {1:F2} ms", Engine.GetFramesPerSecond(), delta * 1000.0));

        if (_trajectory is not { Frames.Count: > 0 } t)
        {
            text.AppendLine("No trajectory loaded.");
            _diagnosticsText.Text = text.ToString();
            return;
        }

        text.AppendLine(string.Format(ci, "Atoms on screen {0} / {1}   bonds {2}", CountAtomsInView(), t.AtomCount, t.Frames[0].Bonds.Length / 2));
        text.AppendLine();

        // Energies are shown for the current trajectory frame, exactly as the core wrote them.
        int index = (int)Math.Floor(_playhead);
        TrajectoryFrame f = t.Frames[index];
        text.AppendLine(string.Format(ci, "Frame {0} / {1}   step {2}", index, t.Frames.Count - 1, f.Step));
        bool real = t.Units == RealUnits;
        string energyUnit = real ? "kJ/mol" : "ε";
        if (real)
        {
            text.AppendLine(string.Format(ci, "Time  {0:F4} ps", f.Time * t.TauFemtoseconds / 1000.0));
        }
        else
        {
            text.AppendLine(string.Format(ci, "Time  {0:F3} τ  ({1:F3} ps argon)", f.Time, f.Time * t.TauFemtoseconds / 1000.0));
        }
        text.AppendLine();
        text.AppendLine($"Energy (from trajectory, {(real ? "kJ/mol" : "units of ε")}):");
        text.AppendLine(string.Format(ci, "  Kinetic    {0,14:F6}", f.KineticEnergy));
        text.AppendLine(string.Format(ci, "  Potential  {0,14:F6}", f.PotentialEnergy));
        text.AppendLine(string.Format(ci, "  Total      {0,14:F6}", f.TotalEnergy));
        text.AppendLine(string.Format(ci, "  Start E0   {0,14:F6}", t.InitialTotalEnergy));
        text.AppendLine();
        // E(t) − E0 is display arithmetic on two values the core computed.
        double drift = f.TotalEnergy - t.InitialTotalEnergy;
        string percent = t.InitialTotalEnergy != 0
            ? string.Format(ci, "{0:+0.0000;-0.0000} %", 100.0 * drift / Math.Abs(t.InitialTotalEnergy))
            : "n/a";
        text.AppendLine("Drift E − E0 (display arithmetic on core values):");
        text.AppendLine(string.Format(ci, "  {0:+0.000E+0;-0.000E+0} {1}   {2} of |E0|", drift, energyUnit, percent));

        if (t.Warning != null)
        {
            text.AppendLine();
            text.AppendLine("Warning: " + t.Warning);
        }
        _diagnosticsText.Text = text.ToString();
    }

    /// <summary>Atoms whose display sphere is inside the camera frustum (display information only).</summary>
    private int CountAtomsInView()
    {
        var planes = _camera.GetFrustum();
        int count = 0;
        for (int i = 0; i < _displayPositions.Length; i++)
        {
            Vector3 p = _displayPositions[i];
            bool inside = true;
            foreach (Plane plane in planes)
            {
                // Frustum plane normals point outwards.
                if (plane.DistanceTo(p) > _radii[i])
                {
                    inside = false;
                    break;
                }
            }
            if (inside) count++;
        }
        return count;
    }

    // ---------------------------------------------------------------- loading

    private static string ResolveTrajectoryPath()
    {
        foreach (string arg in OS.GetCmdlineUserArgs())
        {
            if (arg.StartsWith("--", StringComparison.Ordinal)) continue;
            if (Path.IsPathRooted(arg)) return arg;
            // Godot switches the working directory to the project folder, so resolve relative
            // paths against the shell's directory (PWD), where the user typed them.
            string? launchDir = System.Environment.GetEnvironmentVariable("PWD");
            return Path.GetFullPath(arg, string.IsNullOrEmpty(launchDir) ? System.Environment.CurrentDirectory : launchDir);
        }
        string projectDir = ProjectSettings.GlobalizePath("res://");
        return Path.GetFullPath(Path.Combine(projectDir, "..", "out", "argon-8.traj"));
    }

    private void LoadTrajectory(string path)
    {
        _loadedPath = path;
        try
        {
            _trajectory = TrajectoryFile.Load(path);
        }
        catch (Exception e) when (e is IOException or InvalidDataException or UnauthorizedAccessException)
        {
            _trajectory = null;
            ShowMessage($"Could not load trajectory:\n{path}\n\n{e.Message}\n\nRun: godot --path viewer -- <file.traj>");
            return;
        }

        TrajectoryFile t = _trajectory;
        if (t.Frames.Count == 0)
        {
            ShowMessage($"Trajectory has no frames:\n{path}" + (t.Warning != null ? "\n\n" + t.Warning : ""));
            return;
        }
        if (t.Warning != null) ShowMessage(t.Warning, seconds: 6);

        bool real = t.Units == RealUnits;
        _atoms.InstanceCount = t.AtomCount;
        _displayPositions = new Vector3[t.AtomCount];
        _layoutPositions = new double[3 * t.AtomCount];
        _radii = new float[t.AtomCount];
        _colors = new Color[t.AtomCount];
        for (int i = 0; i < t.AtomCount; i++)
        {
            _radii[i] = real ? RealAtomRadius(t.Elements[i]) : ReducedAtomRadius;
            _colors[i] = ElementColor(t.Elements[i]);
            _atoms.SetInstanceColor(i, _colors[i]);
        }

        TrajectoryFrame first = t.Frames[0];
        var box = new Vector3((float)first.BoxX, (float)first.BoxY, (float)first.BoxZ);
        Vector3 centre = box / 2;
        float size = Math.Max(box.X, Math.Max(box.Y, box.Z));
        _camera.LookFrom(centre + new Vector3(0.6f * size, 0.5f * size, 1.6f * size), centre);
        // Argon files keep the original 30 σ/s. Real-unit boxes are a few tens of Å, so move at
        // about half a box side per second there.
        _camera.Speed = real ? Math.Max(0.5f * size, 2f) : 30f;
        _playhead = 0;
        ShowPlayhead(t);
        GetWindow().Title = $"MD Viewer – {Path.GetFileName(path)} ({t.AtomCount} atoms, {t.Frames.Count} frames)";
    }

    /// <summary>Display colours by element label. In "lj-reduced" files labels do not change the physics (all atoms are argon-like LJ).</summary>
    private static Color ElementColor(string element) => element switch
    {
        "Ar" => new Color(0.35f, 0.75f, 0.95f),
        "H" => new Color(0.95f, 0.95f, 0.95f),
        "O" => new Color(0.95f, 0.25f, 0.2f),
        "C" => new Color(0.3f, 0.3f, 0.3f),
        "N" => new Color(0.25f, 0.4f, 0.95f),
        _ => new Color(0.8f, 0.5f, 0.9f),
    };

    // ---------------------------------------------------------------- scene construction

    private void BuildWorld()
    {
        var env = new Godot.Environment
        {
            BackgroundMode = Godot.Environment.BGMode.Color,
            BackgroundColor = new Color(0.06f, 0.07f, 0.09f),
            AmbientLightSource = Godot.Environment.AmbientSource.Color,
            AmbientLightColor = new Color(1, 1, 1),
            AmbientLightEnergy = 0.35f,
        };
        AddChild(new WorldEnvironment { Environment = env });

        var sun = new DirectionalLight3D { LightEnergy = 1.0f };
        sun.RotationDegrees = new Vector3(-50, -30, 0);
        AddChild(sun);

        // 60° vertical (≈ 87° horizontal at 16:10): close to a natural view on a desktop monitor.
        // Godot's default of 75° looks stretched at the edges.
        _camera = new FlyCamera
        {
            Fov = 60f,
            KeepAspect = Camera3D.KeepAspectEnum.Height,
            Near = 0.02f,
            Far = 500f,
            Current = true,
        };
        AddChild(_camera);

        // Instanced spheres: one draw call for all atoms (rendering approach is still Open in the plan).
        // Unit sphere; each instance is scaled to its atom's display radius.
        var sphere = new SphereMesh { Radius = 1f, Height = 2f, RadialSegments = 24, Rings = 12 };
        sphere.Material = new StandardMaterial3D { VertexColorUseAsAlbedo = true, Roughness = 0.5f };
        _atoms = new MultiMesh
        {
            TransformFormat = MultiMesh.TransformFormatEnum.Transform3D,
            UseColors = true,
            Mesh = sphere,
        };
        AddChild(new MultiMeshInstance3D { Multimesh = _atoms });

        // Instanced bond sticks: unit cylinder, scaled and oriented per half-bond.
        var cylinder = new CylinderMesh { TopRadius = 1f, BottomRadius = 1f, Height = 1f, RadialSegments = 12, Rings = 1, CapTop = false, CapBottom = false };
        cylinder.Material = new StandardMaterial3D { VertexColorUseAsAlbedo = true, Roughness = 0.5f };
        _bonds = new MultiMesh
        {
            TransformFormat = MultiMesh.TransformFormatEnum.Transform3D,
            UseColors = true,
            Mesh = cylinder,
        };
        AddChild(new MultiMeshInstance3D { Multimesh = _bonds });

        _boxMesh = new ImmediateMesh();
        _boxOutline = new MeshInstance3D
        {
            Mesh = _boxMesh,
            MaterialOverride = new StandardMaterial3D
            {
                ShadingMode = BaseMaterial3D.ShadingModeEnum.Unshaded,
                AlbedoColor = new Color(0.9f, 0.85f, 0.4f),
            },
        };
        AddChild(_boxOutline);
    }

    private void DrawBox(Vector3 l)
    {
        _drawnBox = l;
        _boxMesh.ClearSurfaces();
        _boxMesh.SurfaceBegin(Mesh.PrimitiveType.Lines);
        var c = new List<Vector3>();
        for (int i = 0; i < 8; i++) c.Add(new Vector3((i & 1) * l.X, ((i >> 1) & 1) * l.Y, ((i >> 2) & 1) * l.Z));
        for (int i = 0; i < 8; i++)
        {
            for (int bit = 0; bit < 3; bit++)
            {
                int j = i | (1 << bit);
                if (j == i) continue;
                _boxMesh.SurfaceAddVertex(c[i]);
                _boxMesh.SurfaceAddVertex(c[j]);
            }
        }
        _boxMesh.SurfaceEnd();
    }

    private void BuildUi()
    {
        var layer = new CanvasLayer();
        AddChild(layer);

        // F1 diagnostics panel, hidden by default.
        _diagnostics = new PanelContainer { Visible = false, Position = new Vector2(12, 12) };
        _diagnosticsText = new Label();
        _diagnosticsText.AddThemeFontSizeOverride("font_size", 15);
        var margin = new MarginContainer();
        foreach (string side in new[] { "left", "right", "top", "bottom" }) margin.AddThemeConstantOverride("margin_" + side, 10);
        margin.AddChild(_diagnosticsText);
        _diagnostics.AddChild(margin);
        layer.AddChild(_diagnostics);

        // Esc menu: a dimmed full-screen overlay with a centred panel.
        _menu = new Control { Visible = false };
        _menu.SetAnchorsPreset(Control.LayoutPreset.FullRect);
        var shade = new ColorRect { Color = new Color(0, 0, 0, 0.45f) };
        shade.SetAnchorsPreset(Control.LayoutPreset.FullRect);
        _menu.AddChild(shade);
        var centre = new CenterContainer();
        centre.SetAnchorsPreset(Control.LayoutPreset.FullRect);
        _menu.AddChild(centre);

        var panel = new PanelContainer();
        panel.AddThemeStyleboxOverride("panel", new StyleBoxFlat
        {
            BgColor = new Color(0.13f, 0.14f, 0.17f),
            CornerRadiusTopLeft = 6, CornerRadiusTopRight = 6, CornerRadiusBottomLeft = 6, CornerRadiusBottomRight = 6,
        });
        var pad = new MarginContainer();
        foreach (string side in new[] { "left", "right", "top", "bottom" }) pad.AddThemeConstantOverride("margin_" + side, 20);
        var box = new VBoxContainer();
        box.AddThemeConstantOverride("separation", 10);
        var title = new Label { Text = "Menu", HorizontalAlignment = HorizontalAlignment.Center };
        title.AddThemeFontSizeOverride("font_size", 22);
        var restart = new Button { Text = "Restart", CustomMinimumSize = new Vector2(220, 40) };
        restart.Pressed += Restart;
        var quit = new Button { Text = "Quit", CustomMinimumSize = new Vector2(220, 40) };
        quit.Pressed += () => GetTree().Quit();
        var hint = new Label { Text = "Esc to close", HorizontalAlignment = HorizontalAlignment.Center, Modulate = new Color(1, 1, 1, 0.6f) };
        box.AddChild(title);
        box.AddChild(restart);
        box.AddChild(quit);
        box.AddChild(hint);
        pad.AddChild(box);
        panel.AddChild(pad);
        centre.AddChild(panel);
        layer.AddChild(_menu);

        // Centre message for load errors and warnings.
        _message = new Label
        {
            Visible = false,
            HorizontalAlignment = HorizontalAlignment.Center,
            VerticalAlignment = VerticalAlignment.Center,
            AutowrapMode = TextServer.AutowrapMode.WordSmart,
        };
        _message.SetAnchorsPreset(Control.LayoutPreset.FullRect);
        _message.AddThemeFontSizeOverride("font_size", 18);
        layer.AddChild(_message);
    }

    private async void ShowMessage(string text, double seconds = 0)
    {
        _message.Text = text;
        _message.Visible = true;
        if (seconds <= 0) return;
        await ToSignal(GetTree().CreateTimer(seconds), SceneTreeTimer.SignalName.Timeout);
        if (_message.Text == text) _message.Visible = false;
    }
}
