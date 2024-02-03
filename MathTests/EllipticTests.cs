using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MathTests
{
    [TestFixture, Category("Functions")]
    public  class EllipticTests
    {
        [TestCase(0.0, 1.570796326794897, 15)]
        [TestCase(0.5, 1.854074677301372, 15)]
        [TestCase(0.9, 2.578092113348173, 15)]
        public void EllipticKExact(double n, double expected, int decimalPlaces)
        {
            AssertHelpers.AlmostEqualRelative(expected, Mathematics.SpecialFunctions.EllipticK(n), decimalPlaces);
        }

        [TestCase(0.0, 1.570796326794897, 15)]
        [TestCase(0.5, 1.350643881047675, 15)]
        [TestCase(0.9, 1.104774732704073, 15)]
        public void EllipticEExact(double n, double expected, int decimalPlaces)
        {
            AssertHelpers.AlmostEqualRelative(expected, Mathematics.SpecialFunctions.EllipticE(n), decimalPlaces);
        }
    }
}
