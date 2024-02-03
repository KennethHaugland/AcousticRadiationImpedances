using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace Acoustics.RadiationImpedance    
{
    public static partial class Piston
    {

        /// <summary>
        /// Calculates radiation impedance for a Rectangular piston placed in an infinite baffel.
        /// </summary>
        /// <param name="k">Wavenumber 2*pi*f/c_0 of surrounding propagation medium</param>
        /// <param name="a">Long side of rectangle</param>
        /// <param name="b">Short side of rectangle</param>
        /// <param name="Z_0">Specific impedance of surrounding propagation medium. If not specified it returns normalized impedance</param>
        /// <returns></returns>
        public static Complex RetangularBaffel(double k, double a, double b, double Z_0 = 1)
        {

            return FieldExcited.Rectangular(k, 0, 0, a, b, Z_0);


            if (k == 0)
                return Complex.Zero;

            // 'a' should be the long side and 'b' the sort side
            double beta = Math.Min(a, b) / Math.Max(a, b);
            double beta2 = beta * beta;
            double ka = k * Math.Max(a, b);
            double kb = ka * beta;

            // Set upper iteration count, setting it to
            // zero will do it automatically
            double m = 0;

            // Real part
            double a1 = I_a(ka, Math.Sqrt(1 + beta2), m);
            double a2 = I_a(kb, Math.Sqrt(1 + 1 / beta2), m);

            double sqrBeta = Math.Sqrt(1 + beta2);

            double zs1 = 1 - 2 / (Math.PI * ka * kb) *
                (1 + Math.Cos(ka * sqrBeta)
                + ka * sqrBeta * Math.Sin(ka * sqrBeta)
                - Math.Cos(ka) - Math.Cos(kb))
                + 2 / Math.PI * (a1 + a2);


            // Imaginary part            
            double b1 = I_b(ka, Math.Sqrt(1 + beta2), m);
            double b2 = I_b(kb, Math.Sqrt(1 + 1 / beta2), m);

            double zs2 = 2 / (Math.PI * ka * kb) *
                        (ka + kb + Math.Sin(ka * Math.Sqrt(1 + beta2)) -
                         ka * Math.Sqrt(1 + beta2) * Math.Cos(ka * Math.Sqrt(1 + beta2)) -
                         Math.Sin(ka) - Math.Sin(kb))
                         - 2 / Math.PI * (b1 + b2);

            return Z_0 * new Complex(zs1, zs2);
        }

        /// <summary>
        /// Intermediate integral int_1^B { sqrt(1-1/x^2) * cos(x*A) dx}
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B">Upper integration limit</param>
        /// <param name="m"></param>
        /// <returns></returns>
        private static double I_a(double A, double B, double m = 0)
        {
            // Upper iteration limit
            if (m == 0)
                m = Math.Ceiling(10 * A) + 4;

            // Precalcualted values
            double B2 = B * B;
            double sqrB = Math.Sqrt(B2 - 1);
            double sqrB3 = sqrB * (B2 - 1);

            // Iteration values
            double C_n = 1d;
            double I_n = sqrB - Math.Acos(1 / B);
            double I_n1;

            // Result
            double ReturnValue = I_n;

            // Estimate error for each iteration step 
            double LastValue;
            double tol = 1e-10;

            // Using Wynn's epsilon algoritm
            //List<double> result = new List<double>() { I_n };
            for (double n = 1; n < m; n++)
            {
                C_n = -C_n * A * A / (2 * n * (2 * n - 1));
                I_n1 = I_n;
                I_n = Math.Pow(B, (2 * n - 2)) / (2 * n + 1) * sqrB3 + (2 * n - 2) / (2 * n + 1) * I_n1;
                LastValue = ReturnValue;
                ReturnValue += C_n * I_n;

                if (Math.Abs(LastValue - ReturnValue) < tol)
                    return ReturnValue;

            }
            return ReturnValue;
        }

        /// <summary>
        /// Intermediate integral int_1^B { sqrt(1-1/x^2) * sin(x*A) dx}
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B">Upper integration limit</param>
        /// <param name="m"></param>
        /// <returns></returns>
        private static double I_b(double A, double B, double m = 0)
        {
            // Upper iteration limit
            if (m == 0)
                m = Math.Ceiling(10 * A) + 4;

            double B2 = B * B;
            double sqrB = Math.Sqrt(B2 - 1);
            double sqrB3 = sqrB * (B2 - 1);

            double C_n = A;

            double I_n = 0.5 * B * sqrB - 0.5 * Math.Log(B + sqrB);
            double ReturnValue = A * I_n;
            double I_n1;
            double LastValue = 0;

            double tol = 1e-10;
            for (double n = 1; n < m; n++)
            {
                C_n = -C_n * A * A / (2 * n * (2 * n + 1));
                I_n1 = I_n;
                I_n = Math.Pow(B, 2d * n - 1d) / (2 * n + 2) * sqrB3 + (2 * n - 1) / (2 * n + 2) * I_n1;
                LastValue = ReturnValue;
                ReturnValue += C_n * I_n;

                if (Math.Abs(LastValue - ReturnValue) < tol)
                    return ReturnValue;

            }

            return ReturnValue;
        }

    }

}
