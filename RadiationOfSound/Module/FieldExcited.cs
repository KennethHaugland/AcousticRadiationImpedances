using OxyPlot.Series;
using System;
using System.Linq.Expressions;
using System.Numerics;
using System.Runtime.CompilerServices;
using System.Security.Cryptography.Xml;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Markup;
using System.Diagnostics;

namespace Acoustics.Radiation
{

    public static class FieldExcited
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


            Complex S0 = SpecialFunctions.StruveH(0d, ka);
            Complex S1 = SpecialFunctions.StruveH(1d, ka);
            
            Complex ImgTerm = (H1 * S0 - H0 * S1);
            Complex NormalizedRadiationImpedance = k0a * (H0 - H1 / ka + 2 * i / (Math.PI * ka * ka) + Math.PI / 2 * (ImgTerm));

            return Z_0 *NormalizedRadiationImpedance;
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="k0">Wave number</param>
        /// <param name="theta"></param>
        /// <param name="phi">azimut angle</param>
        /// <param name="a">Total height of strip (has to be less than 1 meter) </param>
        /// <param name="Z_0">Impedance of medium</param>
        /// <returns></returns>
        public static Complex WideStrip(double k0, double theta, double phi, double a, double Z_0 = 1)
        {
            // This s a more stable calculation
            return Rectangular(k0, theta, phi, a, a * 10, Z_0);

            if (a > 1.6)
                throw new ArgumentOutOfRangeException("The width of the wide strip is too high, numerical integration will not work");

            if (k0 == 0)
                return Complex.Zero;

            double sinth = Math.Sin(theta);
            double costh = Math.Cos(theta);

            double sinphi = Math.Sin(phi);

            double root = Math.Sqrt(1 - sinphi * sinphi * sinth * sinth);
            double s = Math.Sqrt(1 - costh * costh / (1 - sinphi * sinphi * sinth * sinth));
            double b = k0 * a * root;

            // integration resulution
            double dx = b / 1000;

            double range = b / dx;

            Complex C = Complex.Zero;

            // Integration
            for (double i = 0; i < range; i++)
            {
                double x = dx * i;
                var T1 = ((1 - x / b) * Complex.Cos(x * s));
                var T2 = MathNet.Numerics.SpecialFunctions.HankelH2(0, x);
                if (Complex.IsNaN(T2))
                    T2 = 0;

                C += (T1) * (T2);
            }

            C *= dx;
            return Z_0 * C / root;
        }

        private static Complex setting(double b, double t)
        {
            double b2 = b * b;
            double b3 = b2 * b;
            double b4 = b3 * b;
            double b5 = b4 * b;
            double b6 = b5 * b;
            double b7 = b6 * b;

            Complex i = new Complex(0, 1);
            return
Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * (-1.697916666666667e+2)
- Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 9.418489583333333e+2
+ Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 1.697916666666667e+2 * i
+ Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 9.418489583333333e+2 * i
- b2 * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * 5.272569444444444e+1
- b2 * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 2.117508680555556e+2
+ b3 * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * 1.003559027777778e+1
+ b3 * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 4.296128472222222e+1
- b4 * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * 1.126458333333333
- b4 * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 5.214520833333333
+ b5 * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * 6.929861111111111e-2
+ b5 * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 3.503013888888889e-1
- b6 * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * 1.811111111111111e-3
- b6 * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 1.004638888888889e-2
+ b2 * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 5.272569444444444e+1 * i
+ b2 * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 2.117508680555556e+2 * i
- b3 * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 1.003559027777778e+1 * i
- b3 * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 4.296128472222222e+1 * i
+ b4 * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 1.126458333333333 * i
+ b4 * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 5.214520833333333 * i
- b5 * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 6.929861111111111e-2 * i
- b5 * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 3.503013888888889e-1 * i
+ b6 * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 1.811111111111111e-3 * i
+ b6 * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 1.004638888888889e-2 * i
- Math.Sin(Math.Sin(t) * 5.0) * Math.Sin(t) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * 1.630859375e+3
+ Math.Sin(Math.Sin(t) * 5.0) * Math.Sin(t) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 4.244791666666667e+2
+ Math.Sin(Math.Sin(t) * 5.0) * Math.Sin(t) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 1.630859375e+3 * i
- Math.Sin(Math.Sin(t) * 5.0) * Math.Sin(t) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 4.244791666666667e+2 * i
+ Math.Cos(t * 2.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * 1.302083333333333e+2
+ Math.Cos(t * 2.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 9.244791666666667e+2
- Math.Cos(t * 4.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 8.138020833333333e+1
- Math.Cos(t * 2.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 1.302083333333333e+2 * i
- Math.Cos(t * 2.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 9.244791666666667e+2 * i
+ Math.Cos(t * 4.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 8.138020833333333e+1 * i
+ Math.Sin(t * 3.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * 3.662109375e+2
- Math.Sin(t * 3.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 6.510416666666667e+1
- Math.Sin(t * 5.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * 8.138020833333333
- Math.Sin(t * 3.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 3.662109375e+2 * i
+ Math.Sin(t * 3.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 6.510416666666667e+1 * i
+ Math.Sin(t * 5.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 8.138020833333333 * i
+ b * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * 1.632994791666667e+2
+ b * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 6.254401041666667e+2
- b * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 1.632994791666667e+2 * i
                - b * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 6.254401041666667e+2 * i
                - b * Math.Cos(t * 2.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * 1.3359375e+2
                - b * Math.Cos(t * 2.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 6.236979166666667e+2
                + b * Math.Cos(t * 4.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * (6.25e+2 / 3.84e+2)
                + b * Math.Cos(t * 4.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 5.696614583333333e+1
                + b * Math.Cos(t * 2.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 1.3359375e+2 * i
                + b * Math.Cos(t * 2.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 6.236979166666667e+2 * i
                - b * Math.Cos(t * 4.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 1.627604166666667 * i
                - b * Math.Cos(t * 4.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 5.696614583333333e+1 * i
                - b * Math.Sin(t * 3.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * 2.537434895833333e+2
                + b * Math.Sin(t * 3.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 5.859375e+1
                + b * Math.Sin(t * 5.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * 5.696614583333333
                + b * Math.Sin(t * 3.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 2.537434895833333e+2 * i
                - b * Math.Sin(t * 3.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 5.859375e+1 * i
                - b * Math.Sin(t * 5.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 5.696614583333333 * i
                - b2 * Math.Sin(Math.Sin(t) * 5.0) * Math.Sin(t) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * 3.701102430555556e+2
                + b2 * Math.Sin(Math.Sin(t) * 5.0) * Math.Sin(t) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 1.193229166666667e+2
                + b3 * Math.Sin(Math.Sin(t) * 5.0) * Math.Sin(t) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * 7.490885416666667e+1
                - b3 * Math.Sin(Math.Sin(t) * 5.0) * Math.Sin(t) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 2.307986111111111e+1
                - b4 * Math.Sin(Math.Sin(t) * 5.0) * Math.Sin(t) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * 9.068229166666667
                + b4 * Math.Sin(Math.Sin(t) * 5.0) * Math.Sin(t) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 2.650416666666667
                + b5 * Math.Sin(Math.Sin(t) * 5.0) * Math.Sin(t) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * 6.077430555555556e-1
                - b5 * Math.Sin(Math.Sin(t) * 5.0) * Math.Sin(t) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 1.6775e-1
                - b6 * Math.Sin(Math.Sin(t) * 5.0) * Math.Sin(t) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * 1.739583333333333e-2
                + b6 * Math.Sin(Math.Sin(t) * 5.0) * Math.Sin(t) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 4.527777777777778e-3
                + b2 * Math.Sin(Math.Sin(t) * 5.0) * Math.Sin(t) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 3.701102430555556e+2 * i
                - b2 * Math.Sin(Math.Sin(t) * 5.0) * Math.Sin(t) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 1.193229166666667e+2 * i
                - b3 * Math.Sin(Math.Sin(t) * 5.0) * Math.Sin(t) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 7.490885416666667e+1 * i
                + b3 * Math.Sin(Math.Sin(t) * 5.0) * Math.Sin(t) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 2.307986111111111e+1 * i
                + b4 * Math.Sin(Math.Sin(t) * 5.0) * Math.Sin(t) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 9.068229166666667 * i
                - b4 * Math.Sin(Math.Sin(t) * 5.0) * Math.Sin(t) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 2.650416666666667 * i
                - b5 * Math.Sin(Math.Sin(t) * 5.0) * Math.Sin(t) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 6.077430555555556e-1 * i
                + b5 * Math.Sin(Math.Sin(t) * 5.0) * Math.Sin(t) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 1.6775e-1 * i
                + b6 * Math.Sin(Math.Sin(t) * 5.0) * Math.Sin(t) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 1.739583333333333e-2 * i
                - b6 * Math.Sin(Math.Sin(t) * 5.0) * Math.Sin(t) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 4.527777777777778e-3 * i
                + b2 * Math.Cos(t * 2.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * 4.211805555555556e+1
                + b2 * Math.Cos(t * 2.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 2.105034722222222e+2
                - b3 * Math.Cos(t * 2.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * 7.878472222222222
                - b2 * Math.Cos(t * 4.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * (1.25e+2 / 2.88e+2)
                - b3 * Math.Cos(t * 2.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 4.251736111111111e+1
                + b4 * Math.Cos(t * 2.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * (7.0 / 8.0)
                - b2 * Math.Cos(t * 4.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 1.898871527777778e+1
                + b3 * Math.Cos(t * 4.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * (2.5e+1 / 3.84e+2)
                + b4 * Math.Cos(t * 2.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 5.139583333333333
                - b5 * Math.Cos(t * 2.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * 5.347222222222222e-2
                + b3 * Math.Cos(t * 4.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 3.797743055555556
                - (b4 * Math.Cos(t * 4.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0)) / 1.92e+2
                - b5 * Math.Cos(t * 2.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 3.443055555555556e-1
                + (b6 * Math.Cos(t * 2.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0)) / 7.2e+2
                - b4 * Math.Cos(t * 4.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * (1.75e+2 / 3.84e+2)
                + (b5 * Math.Cos(t * 4.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0)) / 5.76e+3
                + b6 * Math.Cos(t * 2.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 9.861111111111111e-3
                + b5 * Math.Cos(t * 4.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 3.038194444444444e-2
                - (b6 * Math.Cos(t * 4.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0)) / 1.152e+3
                - b2 * Math.Cos(t * 2.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 4.211805555555556e+1 * i
                - b2 * Math.Cos(t * 2.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 2.105034722222222e+2 * i
                + b3 * Math.Cos(t * 2.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 7.878472222222222 * i
                + b2 * Math.Cos(t * 4.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 4.340277777777778e-1 * i
                + b3 * Math.Cos(t * 2.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 4.251736111111111e+1 * i
                - b4 * Math.Cos(t * 2.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 8.75e-1 * i
                + b2 * Math.Cos(t * 4.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 1.898871527777778e+1 * i
                - b3 * Math.Cos(t * 4.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 6.510416666666667e-2 * i
                - b4 * Math.Cos(t * 2.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 5.139583333333333 * i
                + b5 * Math.Cos(t * 2.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 5.347222222222222e-2 * i
                - b3 * Math.Cos(t * 4.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 3.797743055555556 * i
                + b4 * Math.Cos(t * 4.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 5.208333333333333e-3 * i
                + b5 * Math.Cos(t * 2.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 3.443055555555556e-1 * i
                - b6 * Math.Cos(t * 2.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 1.388888888888889e-3 * i
                + b4 * Math.Cos(t * 4.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 4.557291666666667e-1 * i
                - b5 * Math.Cos(t * 4.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 1.736111111111111e-4 * i
                - b6 * Math.Cos(t * 2.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 9.861111111111111e-3 * i
                - b5 * Math.Cos(t * 4.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 3.038194444444444e-2 * i
                + b6 * Math.Cos(t * 4.0) * Math.Cos(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 8.680555555555556e-4 * i
                + b2 * Math.Sin(t * 3.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * 8.492838541666667e+1
                - b2 * Math.Sin(t * 3.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 1.866319444444444e+1
                - b3 * Math.Sin(t * 3.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * 1.703776041666667e+1
                - b2 * Math.Sin(t * 5.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * 1.898871527777778
                + b3 * Math.Sin(t * 3.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 3.559027777777778
                + b4 * Math.Sin(t * 3.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * 2.048697916666667
                + b3 * Math.Sin(t * 5.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * 3.797743055555556e-1
                - b4 * Math.Sin(t * 3.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * (1.3e+1 / 3.2e+1)
                - b5 * Math.Sin(t * 3.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * (3.5e+1 / 2.56e+2)
                - b4 * Math.Sin(t * 5.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * (3.5e+1 / 7.68e+2)
                + b5 * Math.Sin(t * 3.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 2.569444444444444e-2
                + (b6 * Math.Sin(t * 3.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0)) / 2.56e+2
                + b5 * Math.Sin(t * 5.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * 3.038194444444444e-3
                - (b6 * Math.Sin(t * 3.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0)) / 1.44e+3
                - (b6 * Math.Sin(t * 5.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0)) / 1.152e+4
                - b2 * Math.Sin(t * 3.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 8.492838541666667e+1 * i
                + b2 * Math.Sin(t * 3.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 1.866319444444444e+1 * i
                + b3 * Math.Sin(t * 3.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 1.703776041666667e+1 * i
                + b2 * Math.Sin(t * 5.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 1.898871527777778 * i
                - b3 * Math.Sin(t * 3.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 3.559027777777778 * i
                - b4 * Math.Sin(t * 3.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 2.048697916666667 * i
                - b3 * Math.Sin(t * 5.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 3.797743055555556e-1 * i
                + b4 * Math.Sin(t * 3.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 4.0625e-1 * i
                + b5 * Math.Sin(t * 3.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 1.3671875e-1 * i
                + b4 * Math.Sin(t * 5.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 4.557291666666667e-2 * i
                - b5 * Math.Sin(t * 3.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 2.569444444444444e-2 * i
                - b6 * Math.Sin(t * 3.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 3.90625e-3 * i
                - b5 * Math.Sin(t * 5.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 3.038194444444444e-3 * i
                + b6 * Math.Sin(t * 3.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 6.944444444444444e-4 * i
                + b6 * Math.Sin(t * 5.0) * Math.Sin(Math.Sin(t) * 5.0) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 8.680555555555556e-5 * i
                + b * Math.Sin(Math.Sin(t) * 5.0) * Math.Sin(t) * MathNet.Numerics.SpecialFunctions.BesselJ(0, 5.0) * 1.094622395833333e+3
                - b * Math.Sin(Math.Sin(t) * 5.0) * Math.Sin(t) * MathNet.Numerics.SpecialFunctions.BesselJ(1, 5.0) * 3.653645833333333e+2
                - b * Math.Sin(Math.Sin(t) * 5.0) * Math.Sin(t) * MathNet.Numerics.SpecialFunctions.BesselY(0, 5.0) * 1.094622395833333e+3 * i
                + b * Math.Sin(Math.Sin(t) * 5.0) * Math.Sin(t) * MathNet.Numerics.SpecialFunctions.BesselY(1, 5.0) * 3.653645833333333e+2 * i;

        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="k_0"></param>
        /// <param name="theta"></param>
        /// <param name="psi"></param>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="Z_0"></param>
        /// <returns></returns>
        public static Complex Rectangular(double k_0, double theta, double psi, double a, double b, double Z_0 = 1)
        {
            if (k_0 == 0)
                return Complex.Zero;

            // defining some constant terms
            double k0a = k_0 * a;
            double k0b = k_0 * b;
            double mu_x = Math.Sin(theta) * Math.Cos(psi);
            double mu_y = Math.Sin(theta) * Math.Sin(psi);
            double arctan = Math.Atan(b / a);

            Complex U = k0a * k0b;

            // Defining some functions that is integrated over the angle
            Complex V(double phi)
            {
                return -(k0a * Math.Sin(phi) + k0b * Math.Cos(phi));
            };

            Complex W(double phi)
            {
                return 0.5 * Math.Sin(2 * phi);
            };

            Complex alpha(double phi)
            {
                return mu_x * Math.Cos(phi);
            };

            Complex beta(double phi)
            {
                return mu_y * Math.Sin(phi);
            };

            // Analytical solution to the intermediate integral
            Complex I_R(Complex R, Complex C_n, Complex U, Complex V, Complex W)
            {
                Complex C_n2 = C_n * C_n;
                Complex C_n3 = C_n2 * C_n;
                Complex i = new Complex(0, 1);
                Complex expR = Complex.Exp(i * C_n * R);

                Complex result =
                    (expR * ( - U * i * C_n2
                              + V * (C_n - i * C_n2 * R)
                              + W * (-i * C_n2 * R * R + 2 * C_n * R + 2 * i)
                                )
                        + (U * i * C_n2 - V * C_n - W * 2 * i)
                    ) / (C_n3 * 4);

                return result;

            }

            // Default start value
            Complex result = Complex.Zero;

            // Integration resolution
            double deltaPhi = 0.01;


            // Integral - part one
            for (double i = 0; i < arctan; i += deltaPhi)
            {
                for (int n = 0; n < 4; n++)
                {
                    // Generate sequence 1, -1, 1, -1
                    double mod1 = (1 & n) == 0 ? 1 : -1;

                    // Generate sequence 1, 1, -1, -1
                    double mod2 = (2 & n) == 0 ? 1 : -1;

                    Complex C_n = mod1 * alpha(i) + mod2 * beta(i) - 1;
                    result += I_R(k0a / Math.Cos(i), C_n, U, V(i), W(i));
                }
            }

            // Integral - part two
            for (double j = arctan; j < Math.PI / 2; j += deltaPhi)
            {
                for (int n = 0; n < 4; n++)
                {
                    // Generate sequence 1, -1, 1, -1
                    double mod1 = (1 & n) == 0 ? 1 : -1;

                    // Generate sequence 1, 1, -1, -1
                    double mod2 = (2 & n) == 0 ? 1 : -1;

                    Complex C_n = mod1 * alpha(j) + mod2 * beta(j) - 1;
                    result += I_R(k0b / Math.Sin(j), C_n, U, V(j), W(j));
                }
            }

            return Z_0 * result * deltaPhi * new Complex(0, 2) / (Math.PI * U);
        }
    }
}
