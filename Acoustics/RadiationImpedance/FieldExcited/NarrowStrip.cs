using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace Acoustics.RadiationImpedance
{
    public static partial class FieldExcited
    {

        /// <summary>
        /// Returns the radiation impedance of an acoustical narrow strip
        /// </summary>
        /// <param name="k0">Wave number</param>
        /// <param name="theta">Normal incident wave at 0 degree and maximum is 90 degree</param>
        /// <param name="phi">Azimut angle</param>
        /// <param name="a">Strip height</param>
        /// <param name="Z_0">Specific impedance of sourronding medium. Returns normalzed impedance if not set</param>
        /// <returns></returns>
        public static Complex NarrowStrip(double k0, double theta, double phi, double a, double Z_0 = 1)
        {
            double sin_th2 = Math.Sin(theta) * Math.Sin(theta);
            double cos_ph2 = Math.Cos(phi) * Math.Cos(phi);

            double sqr = Math.Sqrt(1 - sin_th2 * cos_ph2);

            if (k0 == 0)
                return new Complex(0, 0);

            double k0a = k0 * a;
            double ka = k0a * sqr;

            Complex i = new Complex(0, 1);

            Complex H0 = MathNet.Numerics.SpecialFunctions.HankelH2(0, ka);
            Complex H1 = MathNet.Numerics.SpecialFunctions.HankelH2(1, ka);


            Complex S0 = Mathematics.SpecialFunctions.StruveH(0d, ka);
            Complex S1 = Mathematics.SpecialFunctions.StruveH(1d, ka);

            Complex ImgTerm = (H1 * S0 - H0 * S1);
            Complex NormalizedRadiationImpedance = k0a * (H0 - H1 / ka + 2 * i / (Math.PI * ka * ka) + Math.PI / 2 * (ImgTerm));

            return Z_0 * NormalizedRadiationImpedance;
        }
    }
}
