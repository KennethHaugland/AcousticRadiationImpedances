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
        /// Return the Complete Elliptic integral of the 1st kind
        /// </summary>
        /// <param name="k">K(k^2) absolute value has to be below 1</param>
        /// <returns></returns>
        /// <remarks>Abramowitz and Stegun p.591, formula 17.3.11</remarks>
        public static double EllipticK(double k)
        {
            return EllipticK(90, k);
        }

        /// <summary>
        /// Return the Complete Elliptic integral of the 2nd kind
        /// </summary>
        /// <param name="k">E(k^2) absolute value has to be below 1</param>
        /// <returns></returns>
        /// <remarks>Abramowitz and Stegun p.591, formula 17.3.12</remarks>
        public static double EllipticE(double k)
        {
            return EllipticE(90, k);
        }

        /// <summary>
        /// Returns the imcomplete elliptic integral of the first kind 
        /// </summary>
        /// <param name="angle">In degrees, valid value range is from 0 to 90</param>
        /// <param name="k">This function thakes k^2 as the parameter</param>
        /// <returns></returns>
        /// <remarks></remarks>
        public static double EllipticK(double angle, double k)
        {
            if (k < 0 || k > 1)
                throw new ArgumentException("Value has to be within range: 0 - 1");

            if (k == 1)
                return double.PositiveInfinity;

            double ang = Math.PI * angle / 180;
            double Sin = Math.Sin(ang);
            double Sin2 = Sin * Sin;
            double Cos = Math.Cos(ang);
            double Cos2 = Cos * Cos;
            return Sin * RF(Cos2, 1 - k * Sin2, 1);
        }

        /// <summary>
        /// Returns the imcomplete elliptic integral of the second kind 
        /// </summary>
        /// <param name="angle">In degrees, valid value range is from 0 to 90</param>
        /// <param name="k">This function thakes k^2 as the parameter</param>
        /// <returns></returns>
        /// <remarks></remarks>
        public static double EllipticE(double angle, double k)
        {
            if (k < 0 || k > 1)
                throw new ArgumentException("Value has to be within range: 0 - 1");

            if (k == 1)
                return 1;

            double ang = Math.PI * angle / 180;
            double Cos2 = Math.Cos(ang) * Math.Cos(ang);
            double Sin = Math.Sin(ang);
            double Sin2 = Sin * Sin;
            double Sin3 = Sin2 * Sin;

            return Sin * RF(Cos2, 1 - k * Sin2, 1) + -1d / 3d * k * Sin3 * RD(Cos2, 1 - k * Sin2, 1);
        }


        /// <summary>
        /// Computes the R_F from Carlson symmetric form
        /// </summary>
        /// <param name="X"></param>
        /// <param name="Y"></param>
        /// <param name="Z"></param>
        /// <returns></returns>
        /// <remarks>http://en.wikipedia.org/wiki/Carlson_symmetric_form#Series_Expansion</remarks>
        private static double RF(double X, double Y, double Z)
        {
            double result;
            double A;
            double lamda;
            double dx;
            double dy;
            double dz;
            double MinError = 1E-10;

            do
            {
                lamda = Math.Sqrt(X * Y) + Math.Sqrt(Y * Z) + Math.Sqrt(Z * X);

                X = (X + lamda) / 4;
                Y = (Y + lamda) / 4;
                Z = (Z + lamda) / 4;

                A = (X + Y + Z) / 3;

                dx = (1 - X / A);
                dy = (1 - Y / A);
                dz = (1 - Z / A);

            } while (Math.Max(Math.Max(Math.Abs(dx), Math.Abs(dy)), Math.Abs(dz)) > MinError);

            double E2 = dx * dy + dy * dz + dz * dx;
            double E3 = dy * dx * dz;

            //http://dlmf.nist.gov/19.36#E1
            result = 1 - (1 / 10) * E2 + (1 / 14) * E3 + (1 / 24) * Math.Pow(E2, 2) - (3 / 44) * E2 * E3 - (5 / 208) * Math.Pow(E2, 3) + (3 / 104) * Math.Pow(E3, 2) + (1 / 16) * Math.Pow(E2, 2) * E3;

            result *= (1 / Math.Sqrt(A));
            return result;

        }

        /// <summary>
        /// Computes the R_D from Carlson symmetric form
        /// </summary>
        /// <param name="X"></param>
        /// <param name="Y"></param>
        /// <param name="Z"></param>
        /// <returns>Construced from R_J(x,y,z,z) which is equal to R_D(x,y,z)</returns>
        /// <remarks>http://en.wikipedia.org/wiki/Carlson_symmetric_form#Series_Expansion</remarks>
        private static double RD(double X, double Y, double Z)
        {
            double sum = 0;
            double fac = 0;
            double lamda = 0;
            double dx = 0;
            double dy = 0;
            double dz = 0;
            double A = 0;
            double MinError = 0;
            MinError = 1E-10;

            sum = 0;
            fac = 1;

            do
            {
                lamda = Math.Sqrt(X * Y) + Math.Sqrt(Y * Z) + Math.Sqrt(Z * X);
                sum += fac / (Math.Sqrt(Z) * (Z + lamda));

                fac /= 4;

                X = (X + lamda) / 4;
                Y = (Y + lamda) / 4;
                Z = (Z + lamda) / 4;

                A = (X + Y + 3 * Z) / 5;

                dx = (1 - X / A);
                dy = (1 - Y / A);
                dz = (1 - Z / A);

            } while (Math.Max(Math.Max(Math.Abs(dx), Math.Abs(dy)), Math.Abs(dz)) > MinError);

            double E2 = 0;
            double E3 = 0;
            double E4 = 0;
            double E5 = 0;
            double result = 0;
            E2 = dx * dy + dy * dz + 3 * Math.Pow(dz, 2) + 2 * dz * dx + dx * dz + 2 * dy * dz;
            E3 = Math.Pow(dz, 3) + dx * Math.Pow(dz, 2) + 3 * dx * dy * dz + 2 * dy * Math.Pow(dz, 2) + dy * Math.Pow(dz, 2) + 2 * dx * Math.Pow(dz, 2);
            E4 = dy * Math.Pow(dz, 3) + dx * Math.Pow(dz, 3) + dx * dy * Math.Pow(dz, 2) + 2 * dx * dy * Math.Pow(dz, 2);
            E5 = dx * dy * Math.Pow(dz, 3);

            //http://dlmf.nist.gov/19.36#E2
            result = (1 - (3 / 14) * E2 + (1 / 6) * E3 + (9 / 88) * Math.Pow(E2, 2) - (3 / 22) * E4 - (9 / 52) * E2 * E3 + (3 / 26) * E5 - (1 / 16) * Math.Pow(E2, 3) + (3 / 40) * Math.Pow(E3, 2) + (3 / 20) * E2 * E4 + (45 / 272) * Math.Pow(E2, 2) * E3 - (9 / 68) * (E3 * E4 + E2 * E5));

            result = 3.0 * sum + fac * result / (A * Math.Sqrt(A));
            return result;

        }
    }    
}
