using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace Mathematics
{
    public static partial class SpecialFunctions
    {
        public static Complex AsymptoticHankel01(double n, double z)
        {
            Complex i = new Complex(0, 1);
            return Complex.Sqrt(2 / (Math.PI * z)) * Complex.Exp(i * (z - n * Math.PI / 2 - Math.PI / 4));
        }


        public static Complex SphericalHankel01(double n, Complex z)
        {
            Complex i = new Complex(0, 1);
            return MathNet.Numerics.SpecialFunctions.SphericalBesselJ(n, z) + i * MathNet.Numerics.SpecialFunctions.SphericalBesselY(n, z);
        }


        public static Complex SphericalHankel01Derivative(double n, Complex z)
        {
            return sphericalDerivative(n, z, (n, z) => { return SphericalHankel01(n, z); });
        }


        public static Complex SphericalHankel02Derivative(double n, Complex z)
        {
            return sphericalDerivative(n, z, (n, z) => { return SphericalHankel02(n, z); });
        }

        public static Complex SphericalHankel02(double n, Complex z)
        {
            Complex i = new Complex(0, -1);
            return MathNet.Numerics.SpecialFunctions.SphericalBesselJ(n, z) + i * MathNet.Numerics.SpecialFunctions.SphericalBesselY(n, z);
        }

        private static Complex sphericalDerivative(double n, Complex z, Func<double, Complex, Complex> f)
        {
            return -f(n + 1, z) + (n / z) * f(n, z);
        }
    }
    
}
