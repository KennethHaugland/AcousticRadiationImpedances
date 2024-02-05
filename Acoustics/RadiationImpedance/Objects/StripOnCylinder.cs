using MathNet.Numerics.Integration;
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

        public static Complex StripOnCylinder(double k0, double r, double v0, double n, double Z0 = 1)
        { 
            if (k0==0) return Complex.Zero;

            Complex kr = Complex.Sqrt(k0 * k0 - k0 * k0 / r / r* n * n);
            Complex i = new Complex(0, 1);
            Complex sum = Complex.Zero;

            int gn = 20;

            double dt = 1;

            Complex derivativeHankelH2(double nn, Complex kra)
            {
                return 0.5 * (
                        MathNet.Numerics.SpecialFunctions.HankelH2(nn - 1, kra)
                        - MathNet.Numerics.SpecialFunctions.HankelH2(nn + 1, kra));
            }

            for (int m = 0 ; m < gn+1; m++) 
            {
                double sin = m==0?1:Math.Sin(m * v0) / (m * v0);
                sum += dt *
                    MathNet.Numerics.SpecialFunctions.HankelH2(m, kr * r) / derivativeHankelH2(m, kr * r) * sin * sin;                   

                if (dt == 1)
                    dt = 2;
            }

            return Z0 * -i * v0 *k0/ kr / Math.PI * sum; 
        
        }
    }
    
}
