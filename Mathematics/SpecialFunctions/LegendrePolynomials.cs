using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Mathematics
{
    public static partial class SpecialFunctions
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

            double[] P = new double[(int)n+1];
            P[0] = 1;
            P[1] = x;

            for (int i = 2; i < P.Length; i++)
                P[i] = ((x * (2 * (double)i - 1) * P[i - 1] - ((double)i -1) * P[i - 2]) / (double)i);

            return P;
        }
    }
}
