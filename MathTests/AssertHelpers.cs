using MathNet.Numerics;

namespace MathTests
{
    internal static class AssertHelpers
    {
        public static void AlmostEqualRelative(double expected, double actual, int decimalPlaces)
        {
            if (double.IsNaN(expected) && double.IsNaN(actual))
            {
                return;
            }

            if (!expected.AlmostEqualRelative(actual, decimalPlaces))
            {
                Assert.Fail("Not equal within {0} places. Expected:{1}; Actual:{2}", decimalPlaces, expected, actual);
            }
        }
    }
}
