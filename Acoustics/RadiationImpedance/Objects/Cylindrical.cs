using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace Acoustics.RadiationImpedance
{
    public static partial class Objects
    {

        public static Complex Cylindrical(double k, double r, double n, double m, double Z_0 = 1)
        {

            if (k==0) return Complex.Zero;

            double Krm = m!= 0 ? Math.Sqrt(k * k - k * k / m / m) : k;
            Complex i = new Complex(0, 1);

            Complex derivativeHankelH2 = 0.5 * (
                    MathNet.Numerics.SpecialFunctions.HankelH2(n - 1, Krm * r) 
                -   MathNet.Numerics.SpecialFunctions.HankelH2(n + 1, Krm * r));

            return -i * k / Krm * Z_0 * MathNet.Numerics.SpecialFunctions.HankelH2(n, Krm * r)/ derivativeHankelH2;
        }
    }
    
}
