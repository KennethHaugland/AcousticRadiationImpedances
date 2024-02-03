using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace Acoustics.RadiationImpedance
{
    public static partial class FieldExcited
    {
        /// <summary>
        /// 
        /// </summary>
        /// <param name="k0">Wave number</param>
        /// <param name="theta"></param>
        /// <param name="phi">azimut angle</param>
        /// <param name="a">Total height of strip (has to be less than 1 meter) </param>
        /// <param name="Z_0">Impedance of medium</param>
        /// <returns></returns>
        public static Complex WideStrip(double k0, double theta, double phi, double a, double Z_0 = 1)
        {
            // This s a more stable calculation
            return Rectangular(k0, theta, phi, a, a * 10, Z_0);

            if (a > 1.6)
                throw new ArgumentOutOfRangeException("The width of the wide strip is too high, numerical integration will not work");

            if (k0 == 0)
                return Complex.Zero;

            double sinth = Math.Sin(theta);
            double costh = Math.Cos(theta);

            double sinphi = Math.Sin(phi);

            double root = Math.Sqrt(1 - sinphi * sinphi * sinth * sinth);
            double s = Math.Sqrt(1 - costh * costh / (1 - sinphi * sinphi * sinth * sinth));
            double b = k0 * a * root;

            // integration resulution
            double dx = b / 1000;

            double range = b / dx;

            Complex C = Complex.Zero;

            // Integration
            for (double i = 0; i < range; i++)
            {
                double x = dx * i;
                var T1 = ((1 - x / b) * Complex.Cos(x * s));
                var T2 = MathNet.Numerics.SpecialFunctions.HankelH2(0, x);
                if (Complex.IsNaN(T2))
                    T2 = 0;

                C += (T1) * (T2);
            }

            C *= dx;
            return Z_0 * C / root;
        }


    }
}
