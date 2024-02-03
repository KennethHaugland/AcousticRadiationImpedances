using Mathematics;
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
        /// <summary>
        /// Calculates radiation impedance for a Circular piston placed in an infinite baffel.
        /// </summary>
        /// <param name="k">Wavenumber 2*pi*f/c_0 of surrounding propagation medium</param>
        /// <param name="a">radius of circle</param>
        /// <param name="Z_0">Specific impedance of surrounding propagation medium. If not specified it returns normalized impedance</param>
        /// <returns></returns>
        public static Complex CircularBaffle(double k, double a, double Z_0 = 1)
        {
            if (k == 0)
                return Complex.Zero;

            double ka = k * a;
            Complex i = new Complex(0, 1);
            return Z_0 * (1 - MathNet.Numerics.SpecialFunctions.BesselJ(1, 2 * ka) / ka + i * SpecialFunctions.StruveH(1, 2 * ka) / ka);
        }
    }    
}
