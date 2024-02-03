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
        public static Complex CircularSupportedBaffle(double k, double a, double Z_0 = 1.0)
        {
            
            if (k == 0.0)
                return Complex.Zero;
            
            // BC
            double[] a_n = new double[3] { 1.0, -2.0, 1.0 }; 

            // Velocity
            double sum = 0.0;
            for (double index = 0.0; index < (double)a_n.Length; ++index)
                sum += a_n[(int)index] / (2.0 * index + 2.0);
            double V = sum * sum * 2.0;

            Complex Sum = Complex.Zero;
            Complex i = new Complex(0.0, 1.0);
            Complex[] S_nm = Pistons.R(k);
            Complex[] T_nm = Pistons.X(k);
            for (double n = 0.0; n < 3.0; ++n)
            {
                for (double m = 0.0; m < 3.0; ++m)
                {
                    // Modal sum 
                    double max = Math.Max(n, m);
                    double min = Math.Min(n, m);
                    Complex PartialSum = Complex.Zero;
                    if (max == 0.0)
                        PartialSum = S_nm[0] + i * T_nm[0];
                    else if (max == 1.0)
                        PartialSum = min != 1.0 ? S_nm[1] + i * T_nm[1] : S_nm[2] + i * T_nm[2];
                    else if (max == 2.0)
                        PartialSum = min != 0.0 ? (min != 1.0 ? S_nm[5] + i * T_nm[5] : S_nm[4] + i * T_nm[4]) : S_nm[3] + i * T_nm[3];

                    Sum += a_n[(int)n] * a_n[(int)m] * PartialSum;
                }
            }
            return Z_0 * Sum / V;
        }
    }
    
}
