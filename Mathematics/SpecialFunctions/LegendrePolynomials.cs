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

            double[] P = new double[(int)n];
            P[0] = 1;
            P[1] = x;

            for (int i = 2; i < n; i++)
                P[i] = ((x * (2 * i + 1) * x * P[i - 1] - i * P[i - 2]) / (i + 1));

            return P;
        }
    }
}
