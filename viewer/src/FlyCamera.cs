using Godot;

namespace MDViewer;

/// <summary>
/// Free-fly camera. W/A/S/D move relative to the view direction, Q/E move down/up,
/// Shift moves faster, the mouse looks around while the cursor is captured.
/// </summary>
public partial class FlyCamera : Camera3D
{
    /// <summary>When false (menu open) the camera ignores input and does not move.</summary>
    public bool Active { get; set; } = true;

    public float Speed { get; set; } = 3.0f;          // σ per second
    public float FastMultiplier { get; set; } = 4.0f;
    public float MouseSensitivity { get; set; } = 0.0025f; // radians per pixel

    private float _yaw;
    private float _pitch;
    private const float MaxPitch = 1.55f; // just under 90°, avoids flipping over

    public void LookFrom(Vector3 position, Vector3 target)
    {
        Position = position;
        Vector3 dir = (target - position).Normalized();
        _yaw = Mathf.Atan2(-dir.X, -dir.Z);
        _pitch = Mathf.Clamp(Mathf.Asin(dir.Y), -MaxPitch, MaxPitch);
        ApplyRotation();
    }

    public override void _UnhandledInput(InputEvent @event)
    {
        if (!Active || Input.MouseMode != Input.MouseModeEnum.Captured) return;
        if (@event is InputEventMouseMotion motion)
        {
            _yaw -= motion.Relative.X * MouseSensitivity;
            _pitch = Mathf.Clamp(_pitch - motion.Relative.Y * MouseSensitivity, -MaxPitch, MaxPitch);
            ApplyRotation();
        }
    }

    public override void _Process(double delta)
    {
        if (!Active) return;

        Basis basis = GlobalTransform.Basis;
        Vector3 forward = -basis.Z;
        Vector3 right = basis.X;
        Vector3 dir = Vector3.Zero;
        if (Input.IsPhysicalKeyPressed(Key.W)) dir += forward;
        if (Input.IsPhysicalKeyPressed(Key.S)) dir -= forward;
        if (Input.IsPhysicalKeyPressed(Key.D)) dir += right;
        if (Input.IsPhysicalKeyPressed(Key.A)) dir -= right;
        if (Input.IsPhysicalKeyPressed(Key.E)) dir += Vector3.Up;
        if (Input.IsPhysicalKeyPressed(Key.Q)) dir -= Vector3.Up;
        if (dir == Vector3.Zero) return;

        float speed = Speed * (Input.IsPhysicalKeyPressed(Key.Shift) ? FastMultiplier : 1.0f);
        Position += dir.Normalized() * speed * (float)delta;
    }

    private void ApplyRotation() => Rotation = new Vector3(_pitch, _yaw, 0);
}
