using MathNet.Numerics;
using MathNet.Numerics.Financial;
using System;
using System.Collections.Generic;
using System.Diagnostics.Eventing.Reader;
using System.Linq;
using System.Numerics;
using System.Security.Cryptography.X509Certificates;

namespace Acoustics
{
    public static class SpecialFunctions
    {
        /// <summary>
        /// 
        /// </summary>
        /// <param name="n"></param>
        /// <param name="x"></param>
        /// <returns></returns>
        public static double[] LegendrePolynomials(double n, double x)
        {
            if (n == 0)
                return new double[] { 1 };

            if (n == 1)
                return new double[] { 1, x };

            double[] P = new double[(int)n];
            P[0] = 1;
            P[1] = x;

            for (int i = 2; i < n; i++)          
                P[i] = ((x * (2 * i + 1) * x * P[i - 1] - i * P[i - 2]) / (i + 1));
            
            return P;
        }

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
            { return sphericalDerivative(n, z, (n, z) => { return SphericalHankel01(n, z); });
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

        private static Complex sphericalDerivative(double n, Complex z, Func<double, Complex,  Complex> f)
        {
            return -f(n + 1, z) + (n / z) * f(n, z);
        }

        /// <summary>
        /// Return the Complete Elliptic integral of the 1st kind
        /// </summary>
        /// <param name="k">K(k^2) absolute value has to be below 1</param>
        /// <returns></returns>
        /// <remarks>Abramowitz and Stegun p.591, formula 17.3.11</remarks>
        public static double EllipticK(double k)
        {
            return EllipticK(90, k);
        }

        /// <summary>
        /// Return the Complete Elliptic integral of the 2nd kind
        /// </summary>
        /// <param name="k">E(k^2) absolute value has to be below 1</param>
        /// <returns></returns>
        /// <remarks>Abramowitz and Stegun p.591, formula 17.3.12</remarks>
        public static double EllipticE(double k)
        {
            return EllipticE(90, k);
        }

        /// <summary>
        /// Returns the imcomplete elliptic integral of the first kind 
        /// </summary>
        /// <param name="angle">In degrees, valid value range is from 0 to 90</param>
        /// <param name="k">This function thakes k^2 as the parameter</param>
        /// <returns></returns>
        /// <remarks></remarks>
        public static double EllipticK(double angle, double k)
        {
            if (k < 0 || k > 1)
                throw new ArgumentException("Value has to be within range: 0 - 1");

            if (k == 1) 
                return double.PositiveInfinity;

            double ang = Math.PI * angle / 180;
            double Sin = Math.Sin(ang);
            double Sin2 = Sin * Sin;
            double Cos = Math.Cos(ang);
            double Cos2 = Cos * Cos;
            return Sin * RF(Cos2, 1 - k * Sin2, 1);   
        }

        /// <summary>
        /// Returns the imcomplete elliptic integral of the second kind 
        /// </summary>
        /// <param name="angle">In degrees, valid value range is from 0 to 90</param>
        /// <param name="k">This function thakes k^2 as the parameter</param>
        /// <returns></returns>
        /// <remarks></remarks>
        public static double EllipticE(double angle, double k)
        {
            if (k < 0 || k > 1)
                throw new ArgumentException("Value has to be within range: 0 - 1");

            if (k == 1)
                return 1;

            double ang = Math.PI * angle / 180;
            double Cos2 = Math.Cos(ang)*Math.Cos(ang);
            double Sin = Math.Sin(ang);
            double Sin2 = Sin*Sin;
            double Sin3 = Sin2 * Sin;

            return Sin * RF(Cos2, 1 - k * Sin2, 1) + -1d / 3d * k * Sin3 * RD(Cos2, 1 - k * Sin2, 1);
        }


        /// <summary>
        /// Computes the R_F from Carlson symmetric form
        /// </summary>
        /// <param name="X"></param>
        /// <param name="Y"></param>
        /// <param name="Z"></param>
        /// <returns></returns>
        /// <remarks>http://en.wikipedia.org/wiki/Carlson_symmetric_form#Series_Expansion</remarks>
        private static double RF(double X, double Y, double Z)
        {
            double result;
            double A;
            double lamda;
            double dx;
            double dy;
            double dz;
            double MinError = 1E-10;

            do
            {
                lamda = Math.Sqrt(X * Y) + Math.Sqrt(Y * Z) + Math.Sqrt(Z * X);

                X = (X + lamda) / 4;
                Y = (Y + lamda) / 4;
                Z = (Z + lamda) / 4;

                A = (X + Y + Z) / 3;

                dx = (1 - X / A);
                dy = (1 - Y / A);
                dz = (1 - Z / A);

            } while (Math.Max(Math.Max(Math.Abs(dx), Math.Abs(dy)), Math.Abs(dz)) > MinError);

            double E2 = dx * dy + dy * dz + dz * dx;
            double E3 = dy * dx * dz;

            //http://dlmf.nist.gov/19.36#E1
            result = 1 - (1 / 10) * E2 + (1 / 14) * E3 + (1 / 24) * Math.Pow(E2, 2) - (3 / 44) * E2 * E3 - (5 / 208) * Math.Pow(E2, 3) + (3 / 104) * Math.Pow(E3, 2) + (1 / 16) * Math.Pow(E2, 2) * E3;

            result *= (1 / Math.Sqrt(A));
            return result;

        }

        /// <summary>
        /// Computes the R_D from Carlson symmetric form
        /// </summary>
        /// <param name="X"></param>
        /// <param name="Y"></param>
        /// <param name="Z"></param>
        /// <returns>Construced from R_J(x,y,z,z) which is equal to R_D(x,y,z)</returns>
        /// <remarks>http://en.wikipedia.org/wiki/Carlson_symmetric_form#Series_Expansion</remarks>
        private static double RD(double X, double Y, double Z)
        {
            double sum = 0;
            double fac = 0;
            double lamda = 0;
            double dx = 0;
            double dy = 0;
            double dz = 0;
            double A = 0;
            double MinError = 0;
            MinError = 1E-10;

            sum = 0;
            fac = 1;

            do
            {
                lamda = Math.Sqrt(X * Y) + Math.Sqrt(Y * Z) + Math.Sqrt(Z * X);
                sum += fac / (Math.Sqrt(Z) * (Z + lamda));

                fac /= 4;

                X = (X + lamda) / 4;
                Y = (Y + lamda) / 4;
                Z = (Z + lamda) / 4;

                A = (X + Y + 3 * Z) / 5;

                dx = (1 - X / A);
                dy = (1 - Y / A);
                dz = (1 - Z / A);

            } while (Math.Max(Math.Max(Math.Abs(dx), Math.Abs(dy)), Math.Abs(dz)) > MinError);

            double E2 = 0;
            double E3 = 0;
            double E4 = 0;
            double E5 = 0;
            double result = 0;
            E2 = dx * dy + dy * dz + 3 * Math.Pow(dz, 2) + 2 * dz * dx + dx * dz + 2 * dy * dz;
            E3 = Math.Pow(dz, 3) + dx * Math.Pow(dz, 2) + 3 * dx * dy * dz + 2 * dy * Math.Pow(dz, 2) + dy * Math.Pow(dz, 2) + 2 * dx * Math.Pow(dz, 2);
            E4 = dy * Math.Pow(dz, 3) + dx * Math.Pow(dz, 3) + dx * dy * Math.Pow(dz, 2) + 2 * dx * dy * Math.Pow(dz, 2);
            E5 = dx * dy * Math.Pow(dz, 3);

            //http://dlmf.nist.gov/19.36#E2
            result = (1 - (3 / 14) * E2 + (1 / 6) * E3 + (9 / 88) * Math.Pow(E2, 2) - (3 / 22) * E4 - (9 / 52) * E2 * E3 + (3 / 26) * E5 - (1 / 16) * Math.Pow(E2, 3) + (3 / 40) * Math.Pow(E3, 2) + (3 / 20) * E2 * E4 + (45 / 272) * Math.Pow(E2, 2) * E3 - (9 / 68) * (E3 * E4 + E2 * E5));

            result = 3.0 * sum + fac * result / (A * Math.Sqrt(A));
            return result;

        }

        /// <summary>
        /// Struve function
        /// </summary>
        /// <param name="n">Order</param>
        /// <param name="z">Argument</param>
        /// <returns>H_n(z)</returns>
        public static Complex StruveH(double n, Complex z)
        {
            if (Complex.Abs(z) < 20)
                // Taylor series
                return StruveHTaylor(n, z);
            else
            {
                // Asymptotic expansion
                return StruveHnNeumannYn(n, z) + MathNet.Numerics.SpecialFunctions.BesselY(n, z);
            }
        }

        /// <summary>
        /// Struve Neumann function that returns H_n(z) - Y_n(z) by asymptotic series
        /// </summary>
        /// <param name="n">Order</param>
        /// <param name="z">Argument</param>
        /// <returns> H_n(z) - Y_n(z). Uses steepest decent evaluation </returns>
        public static Complex StruveHnNeumannYn(double n, Complex z)
        {
            Complex Result = Complex.Zero;

            // This term is actually multiplied with the summation term. I just want to take the logaritm and 
            // try to cancel out all high order terms 
            Complex ConstantTerm = (n - 1) * Complex.Log(0.5 * z)
                            - 0.5 * Complex.Log(Math.PI)
                            - MathNet.Numerics.SpecialFunctions.GammaLn(n + 0.5);

            // The coefficients from H_n(z) => n/Z
            Complex[] Ck = C_k(n / z);

            for (int k = 0; k < Ck.Length; k++)
            {
                // Trying to avoid massive devisions
                Complex item = ConstantTerm + Complex.Log(Ck[k])
                    + MathNet.Numerics.SpecialFunctions.GammaLn(k + 1)
                    - k * Complex.Log(z);

                // Remove logarithm and add the asymptotic series
                Result += Complex.Exp(item);
            }

            return Result;
        }

        /// <summary>
        /// Steepest decent evaluation values 
        /// </summary>
        /// <param name="q">The value n/Z is used for the asymtotic expansion</param>
        /// <returns>All first 10 expansions</returns>
        private static Complex[] C_k(Complex q)
        {
            // I find the formulas easier to read by precalculating the powers
            Complex q2 = q * q;
            Complex q3 = q2 * q;
            Complex q4 = q3 * q;
            Complex q5 = q4 * q;
            Complex q6 = q5 * q;
            Complex q7 = q6 * q;
            Complex q8 = q7 * q;
            Complex q9 = q8 * q;
            Complex q10 = q9 * q;

            return new Complex[]{
                1d,
                2d * q,
                6d * q2 - 0.5d,
                20d * q3 - 4d * q,
                70d * q4 - 45d / 2d * q2 + 3d / 8d,
                252d * q5 - 112d * q3 + 23d / 4d * q,
                924d * q6 - 525d * q4 + 301d / 6d * q2 - 5d / 16d,
                3432d * q7 - 2376d * q5 + 345d * q3 - 22d / 3d * q,
                12870d * q8 - 21021d / 2d * q6 + 16665d / 8d * q4 - 1425d / 16d * q2 + 35d / 128d,
                48620d * q9 - 45760 * q7 + 139139d / 12d * q5 - 1595d / 2d * q3 + 563d / 64d * q,
                184756d * q10 - 196911d * q8 +61061d * q6 - 287289d / 48d * q4 + 133529d / 960d * q2 - 63d / 256d
            };

        }

        /// <summary>
        /// Struve Neumann function that returns H_n(z) - Y_n(z) by asymptotic series
        /// </summary>
        /// <param name="n">Order</param>
        /// <param name="z">Argument</param>
        /// <returns>Obsolete! H_n(z) - Y_n(z). Will only give O(n^-3) accuracy for H_0 and O(n^-6) for H_1</returns>
        private static Complex StruveHnNeumannYn_ver1(double n, Complex z)
        {
            // Avoid negative gamma function
            double k = n + 2;

            Complex sum = 0;
            Complex Z2 = Complex.Log(z / 2);
            for (double m = 0; m < k; m++)
            {
                Complex One = (n - 1 - 2 * m) * Z2;
                Complex Two = -Complex.Log(Math.PI);

                double PossibleNegative = n + 0.5 - m;
                Complex LogValue;

                // Extend to negativ factorial with gamma(z)*gamma(1-z) = pi/sin(pi*z)
                if (PossibleNegative > 0)
                    LogValue = -MathNet.Numerics.SpecialFunctions.GammaLn(n + 0.5 - m);
                else
                    LogValue = Complex.Log(Math.PI)
                    - Complex.Log(Complex.Sin(Math.PI * PossibleNegative))
                    - MathNet.Numerics.SpecialFunctions.GammaLn(1 - PossibleNegative);

                Complex Three = MathNet.Numerics.SpecialFunctions.GammaLn(m + 0.5);
                sum += Complex.Exp(One + Two + LogValue + Three);
            }
            return sum;

        }

        /// <summary>
        /// A straight forward implementation of the Struve function valid for Abs(z) < 16
        /// </summary>
        /// <param name="n">Order</param>
        /// <param name="z">Argument</param>
        /// <returns>Struve function - H_n(z) by Taylor series expansion</returns>
        private static Complex StruveHTaylor(double n, Complex z)
        {            

            // Termwise result
            Complex TermResult = Complex.Zero;
            // For epsilon algorthm
            List<Complex> SummationTerms = new List<Complex>();

            // If zero just return that
            if (z == Complex.Zero) { return Complex.Zero; }


            // Cap the number of iterations for the loop
            double MaxIteration = 30d;

            // Precalculate a value that does not change
            Complex Log_z2 = Complex.Log(z * 0.5d);
            Complex iPi = new Complex(0, Math.PI);

            // Estimated error
            Complex error;

            // Accepted tolerance               
            double tol = 1e-12;

            // Standard Taylor seris implementation except for taking the logaritmic values instead
            for (double m = 0; m < MaxIteration; m += 1d)
            {
                Complex LogarithmicTermEvaluation =
                    // This is just i*pi*m since Log(-1) = i*pi
                    m * iPi
                    // Use precalcualted value
                    + (2 * m + n + 1) * Log_z2
                    // Natural logarithm of the gamma function
                    - MathNet.Numerics.SpecialFunctions.GammaLn(m + 1.5d)
                    - MathNet.Numerics.SpecialFunctions.GammaLn(m + n + 1.5d);

                // The exponential will remove the logarithm
                TermResult += Complex.Exp(LogarithmicTermEvaluation);

                // Termwise results
                SummationTerms.Add(TermResult);

                // Using the Epilon algorithm generally seems to shave
                // off one or two iteration for the same precition
                // Should be a good algorithm to use
                // since its an alternating Taylor series
                Complex Summation = EpsilonAlgorithm(SummationTerms.ToArray());

                // Must have at least two values to use Wynns epsilon algorithm
                if (m > 0)
                {
                    // Assume that the Epsilon algortim improves the summation by at least order of 1.
                    // So the error is estimated by simple substraction
                    error = SummationTerms.Last() - Summation;

                    // Compare magnitude of error with accepted tolerance
                    if (tol > Complex.Abs(error))
                    {
                        //Debug.WriteLine("Number of iterations: " + m.ToString() + " and error " + Complex.Abs(error).ToString("N10") + " with Wynn: ");
                        return Summation;
                    }
                }
            }

            // Tolerance not reached within maximum iteration
            return EpsilonAlgorithm(SummationTerms.ToArray());
        }

        /// <summary>
        /// Peter Wynns epsilon algorithm for calculating accelerated convergence
        /// </summary>
        /// <param name="S_n">The partial sums</param>
        /// <returns>The best accelerated sum it finds</returns>
        public static Complex EpsilonAlgorithm(Complex[] S_n, bool Logaritmic = false)
        {

            int m = S_n.Length;

            Complex[,] r = new Complex[m + 1, m + 1];

            // Fill in the partial sums in the 1 column
            for (int i = 0; i < m; i++)
                r[i, 1] = S_n[i];

            // Epsilon algorithm
            for (int column = 2; column <= m; column++)
            {
                for (int row = 0; row < m - column + 1; row++)
                {
                    //Check for divisions of zero (other checks should be done here)
                    Complex divisor = (r[row + 1, column - 1] - r[row, column - 1]);

                    // Epsilon
                    Complex numerator = 1;

                    if (Logaritmic)
                        numerator = column + 1;

                    if (divisor != 0)
                        r[row, column] = r[row + 1, column - 2] + numerator / divisor;
                    else
                        r[row, column] = 0;
                }
            }

            // Clean up, only interested in the odd number columns
            int items = (int)System.Math.Floor((double)((m + 1) / 2));
            Complex[,] e = new Complex[m, items];

            for (int row = 0; row < m; row++)
            {
                int index = 0;
                for (int column = 1; column < m + 1; column = column + 2)
                {
                    if (row + index == m)
                        break;

                    //e[row + index, index] = r[row, column];
                    e[row, index] = r[row, column];
                    index += 1;
                }
            }
            return e[0, e.GetLength(1) - 1];
        }

        /// <summary>
        /// Peter Wynns epsilon algorithm for calculating accelerated convergence
        /// </summary>
        /// <param name="S_n">The partial sums</param>
        /// <returns>The best accelerated sum it finds</returns>
        public static double EpsilonAlgorithm(List<double> S_n, bool Logaritmic = false)
        {
            return EpsilonAlgorithm(S_n.ToArray() , Logaritmic);
        }
            /// <summary>
            /// Peter Wynns epsilon algorithm for calculating accelerated convergence
            /// </summary>
            /// <param name="S_n">The partial sums</param>
            /// <returns>The best accelerated sum it finds</returns>
            public static double EpsilonAlgorithm(double[] S_n, bool Logaritmic = false)
        {

            int m = S_n.Length;

            double[,] r = new double[m + 1, m + 1];

            // Fill in the partial sums in the 1 column
            for (int i = 0; i < m; i++)
                r[i, 1] = S_n[i];

            // Epsilon algorithm
            for (int column = 2; column <= m; column++)
            {
                for (int row = 0; row < m - column + 1; row++)
                {
                    //Check for divisions of zero (other checks should be done here)
                    double divisor = (r[row + 1, column - 1] - r[row, column - 1]);

                    // Epsilon
                    double numerator = 1;

                    if (Logaritmic)
                        numerator = column + 1;

                    if (divisor != 0)
                        r[row, column] = r[row + 1, column - 2] + numerator / divisor;
                    else
                        r[row, column] = 0;
                }
            }

            // Clean up, only interested in the odd number columns
            int items = (int)System.Math.Floor((double)((m + 1) / 2));
            double[,] e = new double[m, items];

            for (int row = 0; row < m; row++)
            {
                int index = 0;
                for (int column = 1; column < m + 1; column = column + 2)
                {
                    if (row + index == m)
                        break;

                    //e[row + index, index] = r[row, column];
                    e[row, index] = r[row, column];
                    index += 1;
                }
            }
            return e[0, e.GetLength(1) - 1];
        }

        /// <summary>
        /// Euler transform that transforms the alternating series
        /// a_0 into a faster convergence with no negative coefficients
        /// </summary>
        /// <param name="a_0">The alternating power series</param>
        /// <returns></returns>
        public static double[] EulerTransformation(double[] a_0)
        {
            // Each series item
            List<double> a_k = new List<double>();

            // finite difference of each item
            double delta_a_0 = 0;

            for (int k = 0; k < a_0.Length; k++)
            {
                delta_a_0 = 0;
                for (int m = 0; m <= k; m++)
                {
                    double choose_k_over_m = Math.Exp(MathNet.Numerics.SpecialFunctions.GammaLn(k + 1) - MathNet.Numerics.SpecialFunctions.GammaLn(m + 1) - MathNet.Numerics.SpecialFunctions.GammaLn(k - m + 1));
                    //double choose_k_over_m = (double)GetBinCoeff(k, m);
                    delta_a_0 += Math.Pow(-1d, (double)m) *
                                 choose_k_over_m * Math.Abs(a_0[(int)(m)]);
                }

                a_k.Add(Math.Pow(1d / 2d, (double)(k + 1)) * delta_a_0);
            }
            return a_k.ToArray();
        }
    }
}
