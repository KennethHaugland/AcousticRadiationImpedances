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

    }
}

