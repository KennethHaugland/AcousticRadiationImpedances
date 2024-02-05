using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace Acoustics.RadiationImpedance
{
    public static partial class Objects
    {

        public static Complex Spherical(double k, double a, double Mode, Complex Z_0)
        {
            Complex i = new Complex(0, 1);
            if (k == 0) return Complex.Zero;

            return -i * Z_0 * Mathematics.SpecialFunctions.SphericalHankel02(Mode, k*a)/ Mathematics.SpecialFunctions.SphericalHankel02Derivative(Mode,k*a);
        }
    }
}
