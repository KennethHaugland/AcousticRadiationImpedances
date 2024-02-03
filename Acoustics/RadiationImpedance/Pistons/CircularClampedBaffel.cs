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
        public static Complex CircularClampedBaffle(double k, double a, double Z_0 = 1.0)
        {
            if (k == 0.0)
                return Complex.Zero;
            double ka = k * a;

            double[] PistonBoundaryConditions = new double[3] { 1.0, 0.0, 0.0 };
            double[] ClampedBoundaryConditions = new double[3] { 1.0, -2.0, 1.0 };
            double[] SimplySupportedBoundaryCondition = new double[3]{1.0,-1.2453,0.2453};

            double[] a_n = ClampedBoundaryConditions;
            Complex ModalSum = Complex.Zero;
            Complex i = new Complex(0.0, 1.0);
            Complex[] S_nm = Pistons.R(k);
            Complex[] T_nm = Pistons.X(k);
            for (double n = 0.0; n < 3.0; ++n)
            {
                for (double m = 0.0; m < 3.0; ++m)
                {
                    double max = Math.Max(n, m);
                    double min = Math.Min(n, m);
                    Complex PartialResult = Complex.Zero;
                    if (max == 0.0)
                        PartialResult = S_nm[0] + i * T_nm[0];
                    else if (max == 1.0)
                        PartialResult = min != 1.0 ? S_nm[1] + i * T_nm[1] : S_nm[2] + i * T_nm[2];
                    else if (max == 2.0)
                        PartialResult = min != 0.0 ? (min != 1.0 ? S_nm[5] + i * T_nm[5] : S_nm[4] + i * T_nm[4]) : S_nm[3] + i * T_nm[3];
                    ModalSum += a_n[(int)n] * a_n[(int)m] * PartialResult;
                }
            }
            double VelocityCorrection = 0.0;
            for (double index = 0.0; index < (double)a_n.Length; ++index)
                VelocityCorrection += a_n[(int)index] / (2.0 * index + 2.0);
            double V = VelocityCorrection * VelocityCorrection * 2.0;
            return Z_0 * ModalSum / V;

        }
        private static Complex[] R(double b)
        {
            double b2 = b * b;
            double b3 = b2 * b;
            double b4 = b2 * b2;
            double b5 = b * b4;
            double b6 = b2 * b4;
            double b7 = b3 * b4;
            double b8 = b4 * b4;
            double b9 = b5 * b4;
            double J0 = MathNet.Numerics.SpecialFunctions.BesselJ(0.0, 2.0 * b);
            double J1 = MathNet.Numerics.SpecialFunctions.BesselJ(1.0, 2.0 * b);
            return new Complex[6]
            {
        (Complex) (0.5 - J1 / (2.0 * b)),
        (Complex) (0.25 - 1.0 / (2.0 * b2) - J1 / (2.0 * b3) * (b2 - 3.0) - J0 / b2),
        (Complex) (1.0 / 6.0 - 1.0 / (2.0 * b2) - J1 / (2.0 * b5) * (b4 - 10.0 * b2 + 10.0) - J0 / b4 * (2.0 * b2 - 5.0)),
        (Complex) (1.0 / 6.0 - 1.0 / b2 + 6.0 / b4 - J1 / (2.0 * b5) * (b4 - 14.0 * b2 + 40.0) - J0 / b4 * (2.0 * b2 - 14.0)),
        (Complex) (0.125 - 5.0 / (6.0 * b2) + 3.0 / b4 - J1 / (2.0 * b7) * (b6 - 25.0 * b4 + 140.0 * b2 - 140.0) - J0 / b6 * (3.0 * b4 - 32.0 * b2 + 70.0)),
        (Complex) (0.1 - 1.0 / b2 + 4.0 / b4 - J1 / (2.0 * b9) * (b8 - 44.0 * b6 + 504.0 * b4 - 2016.0 * b2 + 2016.0) - J0 / b8 * (4.0 * b6 - 80.0 * b4 + 504.0 * b2 - 1008.0))
            };
        }

        private static Complex[] X(double b)
        {
            Complex complex1 = Mathematics.SpecialFunctions.StruveH(0.0, (Complex)(2.0 * b));
            Complex complex2 = Mathematics.SpecialFunctions.StruveH(1.0, (Complex)(2.0 * b));
            double num1 = b * b;
            double num2 = num1 * b;
            double num3 = num1 * num1;
            double num4 = num3 * b;
            double num5 = num4 * b;
            double num6 = num5 * b;
            double num7 = num6 * b;
            double num8 = num4 * num3;
            double pi = Math.PI;
            return new Complex[6]
            {
        complex2 / (2.0 * b),
        complex2 / (2.0 * num2) * (num1 - 3.0) + complex1 / num1,
        complex2 / (2.0 * num4) * (num3 - 10.0 * num1 + 10.0) + complex1 / num3 * (2.0 * num1 - 5.0) + 20.0 / (3.0 * pi * num2),
        complex2 / (2.0 * num4) * (num3 - 14.0 * num1 + 40.0) + complex1 / num3 * (2.0 * num1 - 14.0) + 8.0 / (3.0 * pi * num2),
        complex2 / (2.0 * num6) * (num5 - 25.0 * num3 + 140.0 * num1 - 140.0) + complex1 / num5 * (3.0 * num3 - 32.0 * num1 + 70.0) + 16.0 / (pi * num2) - 280.0 / (3.0 * pi * num4),
        complex2 / (2.0 * num8) * (num7 - 44.0 * num5 + 504.0 * num3 - 2016.0 * num1 + 2016.0) + complex1 / num7 * (4.0 * num5 - 80.0 * num3 + 504.0 * num1 - 1008.0) + 1.0 / (5.0 * pi * num6) * (160.0 * num3 - 2016.0 * num1 + 6720.0)
            };
        }
    }
}
