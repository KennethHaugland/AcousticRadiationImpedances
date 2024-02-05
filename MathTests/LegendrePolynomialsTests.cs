using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MathTests
{
    [TestFixture, Category("Functions")]
    public class LegendrePolynomialsTests
    {
        [TestCase(1, 0.5, 0.5)]
        [TestCase(5, 0.5, 0.08984375)]
        [TestCase(10, 0.5, -0.188228607177734375)]
        public void LegendrePolynomialsExact(double n, double x, double expected)
        {
            AssertHelpers.AlmostEqualRelative(expected, Mathematics.SpecialFunctions.LegendrePolynomials(n,x).Last(), 10);
        }
    }
}
