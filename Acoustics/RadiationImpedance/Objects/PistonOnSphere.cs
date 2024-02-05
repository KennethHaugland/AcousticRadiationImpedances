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

        public static Complex PistonOnSphere(double k, double theta, double radius, double Z_0 = 1)
        {
            if (k == 0) return Complex.Zero;
            
            // Needs to increase as the angle gets shallower
            int order = 20;

            Complex i = new Complex(0, 1);
            Complex Z_n(double n)
            {
                return -i * Z_0 *
                Mathematics.SpecialFunctions.SphericalHankel02(n, k * radius) /
                Mathematics.SpecialFunctions.SphericalHankel02Derivative(n, k * radius);                
            }

            double cosV_0 = Math.Cos(theta);

            double[] L_0 = Mathematics.SpecialFunctions.LegendrePolynomials(order+1, cosV_0);

            var L = L_0.ToList();
            L.Insert(0, 1);
            Complex Sum = Complex.Zero;
            for (double j = 0; j < L.Count-2; j++)
            {
                var P = L[(int)j] - L[(int)j + 2];
                Sum += Z_n(j) / (2d * j + 1d) * P * P;            
            }
            
            return Z_0 * 0.5 /(1d - cosV_0) *Sum;
        }
    }
}
