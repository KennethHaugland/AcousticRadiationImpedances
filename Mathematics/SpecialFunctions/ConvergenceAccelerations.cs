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
            return EpsilonAlgorithm(S_n.ToArray(), Logaritmic);
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
