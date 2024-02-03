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
        /// <param name="k_0"></param>
        /// <param name="theta"></param>
        /// <param name="psi"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="Z_0"></param>
        /// <returns></returns>
        public static Complex Rectangular(double k_0, double theta, double psi, double a, double b, double Z_0 = 1)
        {
            if (k_0 == 0)
                return Complex.Zero;

            // defining some constant terms
            double k0a = k_0 * a;
            double k0b = k_0 * b;
            double mu_x = Math.Sin(theta) * Math.Cos(psi);
            double mu_y = Math.Sin(theta) * Math.Sin(psi);
            double arctan = Math.Atan(b / a);

            Complex U = k0a * k0b;

            // Defining some functions that is integrated over the angle
            Complex V(double phi)
            {
                return -(k0a * Math.Sin(phi) + k0b * Math.Cos(phi));
            };

            Complex W(double phi)
            {
                return 0.5 * Math.Sin(2 * phi);
            };

            Complex alpha(double phi)
            {
                return mu_x * Math.Cos(phi);
            };

            Complex beta(double phi)
            {
                return mu_y * Math.Sin(phi);
            };

            // Analytical solution to the intermediate integral
            Complex I_R(Complex R, Complex C_n, Complex U, Complex V, Complex W)
            {
                Complex C_n2 = C_n * C_n;
                Complex C_n3 = C_n2 * C_n;
                Complex i = new Complex(0, 1);
                Complex expR = Complex.Exp(i * C_n * R);

                Complex result =
                    (expR * (-U * i * C_n2
                              + V * (C_n - i * C_n2 * R)
                              + W * (-i * C_n2 * R * R + 2 * C_n * R + 2 * i)
                                )
                        + (U * i * C_n2 - V * C_n - W * 2 * i)
                    ) / (C_n3 * 4);

                return result;

            }

            // Default start value
            Complex result = Complex.Zero;

            // Integration resolution
            double deltaPhi = 0.01;


            // Integral - part one
            for (double i = 0; i < arctan; i += deltaPhi)
            {
                for (int n = 0; n < 4; n++)
                {
                    // Generate sequence 1, -1, 1, -1
                    double mod1 = (1 & n) == 0 ? 1 : -1;

                    // Generate sequence 1, 1, -1, -1
                    double mod2 = (2 & n) == 0 ? 1 : -1;

                    Complex C_n = mod1 * alpha(i) + mod2 * beta(i) - 1;
                    result += I_R(k0a / Math.Cos(i), C_n, U, V(i), W(i));
                }
            }

            // Integral - part two
            for (double j = arctan; j < Math.PI / 2; j += deltaPhi)
            {
                for (int n = 0; n < 4; n++)
                {
                    // Generate sequence 1, -1, 1, -1
                    double mod1 = (1 & n) == 0 ? 1 : -1;

                    // Generate sequence 1, 1, -1, -1
                    double mod2 = (2 & n) == 0 ? 1 : -1;

                    Complex C_n = mod1 * alpha(j) + mod2 * beta(j) - 1;
                    result += I_R(k0b / Math.Sin(j), C_n, U, V(j), W(j));
                }
            }

            return Z_0 * result * deltaPhi * new Complex(0, 2) / (Math.PI * U);
        }
    }
    
}
