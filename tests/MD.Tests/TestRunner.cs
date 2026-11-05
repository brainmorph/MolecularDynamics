using System;
using System.Diagnostics;
using System.Linq;
using System.Reflection;

namespace MD.Tests;

/// <summary>Marks a public static parameterless method as a test.</summary>
[AttributeUsage(AttributeTargets.Method)]
public sealed class TestAttribute : Attribute { }

public sealed class AssertionException : Exception
{
    public AssertionException(string message) : base(message) { }
}

public static class Assert
{
    public static void True(bool condition, string message)
    {
        if (!condition) throw new AssertionException(message);
    }

    /// <summary>Passes when |actual − expected| ≤ absTol + relTol·|expected|.</summary>
    public static void Close(double expected, double actual, double absTol, double relTol = 0, string what = "value")
    {
        double tol = absTol + relTol * Math.Abs(expected);
        if (!(Math.Abs(actual - expected) <= tol))
        {
            throw new AssertionException(
                $"{what}: expected {expected:G17}, got {actual:G17} (diff {actual - expected:G3}, tol {tol:G3})");
        }
    }

    public static void Equal<T>(T expected, T actual, string what = "value")
    {
        if (!Equals(expected, actual)) throw new AssertionException($"{what}: expected {expected}, got {actual}");
    }

    public static void Throws<TException>(Action action, string what) where TException : Exception
    {
        try
        {
            action();
        }
        catch (TException)
        {
            return;
        }
        catch (Exception e)
        {
            throw new AssertionException($"{what}: expected {typeof(TException).Name}, got {e.GetType().Name}: {e.Message}");
        }
        throw new AssertionException($"{what}: expected {typeof(TException).Name}, nothing was thrown");
    }
}

/// <summary>
/// Minimal test runner (no third-party test framework, per the project rules).
/// Usage: dotnet run --project tests/MD.Tests [-- filter]
/// Runs every [Test] method whose "Class.Method" name contains the filter. Exit code 1 on any failure.
/// </summary>
public static class TestRunner
{
    public static int Main(string[] args)
    {
        string filter = args.Length > 0 ? args[0] : "";
        var tests = typeof(TestRunner).Assembly.GetTypes()
            .SelectMany(t => t.GetMethods(BindingFlags.Public | BindingFlags.Static))
            .Where(m => m.GetCustomAttribute<TestAttribute>() != null)
            .Select(m => (Name: $"{m.DeclaringType!.Name}.{m.Name}", Method: m))
            .Where(t => t.Name.Contains(filter, StringComparison.OrdinalIgnoreCase))
            .OrderBy(t => t.Name, StringComparer.Ordinal)
            .ToList();

        int failed = 0;
        var total = Stopwatch.StartNew();
        foreach (var (name, method) in tests)
        {
            var sw = Stopwatch.StartNew();
            try
            {
                method.Invoke(null, null);
                Console.WriteLine($"PASS  {name} ({sw.ElapsedMilliseconds} ms)");
            }
            catch (TargetInvocationException e) when (e.InnerException != null)
            {
                failed++;
                Console.WriteLine($"FAIL  {name} ({sw.ElapsedMilliseconds} ms)");
                Console.WriteLine($"      {e.InnerException.GetType().Name}: {e.InnerException.Message}");
                if (e.InnerException is not AssertionException)
                {
                    Console.WriteLine(e.InnerException.StackTrace);
                }
            }
        }
        Console.WriteLine();
        Console.WriteLine($"{tests.Count - failed}/{tests.Count} passed in {total.Elapsed.TotalSeconds:F1} s");
        return failed == 0 && tests.Count > 0 ? 0 : 1;
    }
}
