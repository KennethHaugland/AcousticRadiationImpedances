using Mathematics;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace Acoustics.RadiationImpedance
{
    public static partial class Pistons
    {

        private static double EllipticSumZ_r(double k0a2, double beta, double gn)
        {
            List<double> Summation = new List<double>();
            double C_n = 0.5d * k0a2;
            double beta2 = beta * beta;
            double I_n1 = 1d / beta;
            double I_n = 1d;
            double I_n2;

            double SumReal = 0.5 * k0a2;

            Summation.Add(SumReal);

            for (double n = 2; n < gn; n++)
            {
                C_n = -C_n * k0a2 / (n * (1 + n));
                I_n2 = I_n1;
                I_n1 = I_n;
                I_n = ((2d * n - 3d) / (2d * n - 2d)) * (1 + beta2) * I_n1 -
                        ((2d * n - 4d) / (2d * n - 2d)) * beta2 * I_n2;
                SumReal += C_n * I_n;
                Summation.Add(SumReal);
            }

            return SpecialFunctions.EpsilonAlgorithm(Summation);
        }

        private static double EllipticSumZ_i(double k0a2, double beta, double gn)
        {

            double beta2 = beta * beta;
            double C_n = -16d / 45d * k0a2;

            double I_n1 = SpecialFunctions.EllipticK(1 - beta2);
            double I_n = SpecialFunctions.EllipticE(1 - beta2);
            double I_n2;

            double SumImag = 4d / 3d * I_n1 - 16d / 45d * k0a2 * I_n;

            List<double> Summation = new List<double>() {
                4d / 3d * I_n1,
                SumImag };


            for (double n = 2; n < gn; n++)
            {
                C_n = -4 * C_n * k0a2 / ((1 + 2 * n) * (3 + 2 * n));
                I_n2 = I_n1;
                I_n1 = I_n;
                I_n = (2d * n - 2d) / (2d * n - 1d) * (1 + beta2) * I_n1 -
                        (2d * n - 3d) / (2d * n - 1d) * beta2 * I_n2;
                SumImag += C_n * I_n;
                Summation.Add(SumImag);

            }

            return SpecialFunctions.EpsilonAlgorithm(Summation);
        }

        /// <summary>
        /// Calculates radiation impedance for a Elliptic piston placed in an infinite baffel.
        /// </summary>
        /// <param name="k">Wavenumber 2*pi*f/c_0 of surrounding propagation medium</param>
        /// <param name="a">Long radius of ellipse</param>
        /// <param name="b">Short radius of ellipse</param>
        /// <param name="Z_0">Specific impedance of surrounding propagation medium. If not specified it returns normalized impedance</param>
        /// <returns></returns>
        public static Complex EllipticBaffle(double k, double a, double b, double Z_0 = 1)
        {
            if (k == 0)
                return Complex.Zero;

            double ka = k * a;
            double ka2 = ka * ka;
            double beta = Math.Min(a, b) / Math.Max(a, b);

            double Limit = Math.Ceiling(4 * ka) + 4d;
            double MaxLimit = 80;

            if (Limit < MaxLimit)
            {
                //Taylor series
                double Z_r = beta * (EllipticSumZ_r(ka2, beta, Limit));
                double Z_i = 4 / Math.PI / Math.PI * ka * beta * EllipticSumZ_i(ka2, beta, Limit);
                return Z_0 * new Complex(Z_r, Z_i);
            }
            else
            {
                // Asymptotic expansion
                Complex LongSideBaffel = CircularBaffle(k, a);
                Complex ShortSideBaffel = CircularBaffle(k, b);
                return Z_0 * (LongSideBaffel + ShortSideBaffel) / 2;
            }
        }

    }
}
