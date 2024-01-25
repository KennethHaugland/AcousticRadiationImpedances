using OxyPlot.Legends;
using OxyPlot.Series;
using OxyPlot;
using Prism.Mvvm;
using System;
using System.Collections.Generic;
using System.Numerics;
using System.Linq;
using Prism.Commands;
using MathNet.Numerics;

namespace RadiationOfSound.ViewModels
{
    public class MainWindowViewModel : BindableBase
    {
        private string _title = "Prism Application";
        public string Title
        {
            get { return _title; }
            set { SetProperty(ref _title, value); }
        }

        private double pA =1d;
        public double A
        {
            get { return pA; }
            set { SetProperty(ref pA, value); }
        }

        private double pB = 3d;
        public double B
        {
            get { return pB; }
            set { SetProperty(ref pB, value); }
        }

        private double pTheta = 45d;
        public double Theta
        {
            get { return pTheta; }
            set { SetProperty(ref pTheta, value); }
        }

        private double pPhi = 45d;
        public double Phi
        {
            get { return pPhi; }
            set { SetProperty(ref pPhi, value); }
        }

        private DelegateCommand _ExecuteCalculateCmd;
        public DelegateCommand CalculateCmd  =>
            _ExecuteCalculateCmd ?? (_ExecuteCalculateCmd = new DelegateCommand(ExecuteCalculateCmd));

        public MainWindowViewModel()
        {
            ExecuteCalculateCmd();
        }


        void ExecuteCalculateCmd()
        {
            double theta = Theta * Math.PI / 180;
            double phi = Phi * Math.PI / 180;
            double a = A;
            double b = B;

            List<Complex> RectangularResults = new List<Complex>();

            List<double> MaxList = new List<double>();
            List<double> MinList = new List<double>();

            List<double> MaxListI = new List<double>();
            List<double> MinListI = new List<double>();

            double UpperMode = 10;
            double deltaWaveNumber = 0.2;

            List<Complex> SphericalRadiator = new List<Complex>();
            List<Complex> NStrip = new List<Complex>();
            List<Complex> WStrip = new List<Complex>();
            List<Complex> RStrip = new List<Complex>();
            for (double k = 0; k <= UpperMode; k += deltaWaveNumber)
            {
                var Rectangular = Acoustics.Radiation.FieldExcited.Rectangular(k, theta, phi, a, b);
                //var WideStrip = Acoustics.Radiation.FieldExcited.WideStrip(k, theta, phi, a);
                var NarrowStrip = Acoustics.Radiation.FieldExcited.NarrowStrip(k, theta, phi, a);
                var WideRStrip = Acoustics.Radiation.FieldExcited.Rectangular(k, 0, 0, a, a*10);
                //double[] StripRealDifference = new double[] { WideStrip.Real, NarrowStrip.Real };
                //double[] StripImagDifference = new double[] { WideStrip.Imaginary, NarrowStrip.Imaginary };

                //MaxList.Add(StripRealDifference.Max());
                //MinList.Add(StripRealDifference.Min());

                //MaxListI.Add(StripImagDifference.Max());
                //MinListI.Add(StripImagDifference.Min());

                RectangularResults.Add(Rectangular);
                //SphericalRadiator.Add(new Complex(0,-1)* Acoustics.SpecialFunctions.SphericalHankel02(1,k)/Acoustics.SpecialFunctions.SphericalHankel02Derivative(1,k));
                NStrip.Add(NarrowStrip);
                RStrip.Add(WideRStrip);
            }

            FieldExcitedPlotModel = OxyPlotHelper.PlotArray(RectangularResults, "");           
            AreaModel = OxyPlotHelper.PlotArray(new List<List<Complex>>() {NStrip,WStrip}, ""); //OxyPlotHelper.PlotAreaGraph(MaxList, MinList, MaxListI, MinListI, "");




            List<Complex> CircularResult = new List<Complex>();
            for (double k = 0; k <= UpperMode; k += deltaWaveNumber)
                CircularResult.Add(Acoustics.Radiation.Pistons.CircularBaffle(k, a));
         
            CircularBaffelPlot = OxyPlotHelper.PlotArray(CircularResult, "Circular piston" );

            List<Complex> EllipticResult = new List<Complex>();
            for (double k = 0; k <= UpperMode; k += deltaWaveNumber)
                EllipticResult.Add(Acoustics.Radiation.Pistons.EllipticBaffle(k, a, b));
   
            EllipticPlotModel = OxyPlotHelper.PlotArray(EllipticResult, "Elliptic piston");

            List<Complex> RetangularPistonResult = new List<Complex>();
            for (double k = 0; k <= UpperMode; k += deltaWaveNumber)
            {
                RetangularPistonResult.Add(Acoustics.Radiation.FieldExcited.Rectangular(k, 0, 0, a, b));
                //   RetangularPistonResult.Add(Acoustics.Radiation.Pistons.RetangularBaffel(k, a, b));
            }
            RectangularPistonPlotmodel = OxyPlotHelper.PlotArray(RetangularPistonResult, "Rectangular piston");

            RaisePropertyChanged();
        }

        private PlotModel pPlotModel = new PlotModel();
        public PlotModel FieldExcitedPlotModel
        {
            get { return pPlotModel; }
            set { SetProperty(ref pPlotModel, value); }
        }

        private PlotModel pCircularBaffelPlot = new PlotModel();
        public PlotModel CircularBaffelPlot
        {
            get { return pCircularBaffelPlot; }
            set { SetProperty(ref pCircularBaffelPlot, value); }
        }

        private PlotModel pEllipticPlotmodel = new PlotModel();
        public PlotModel EllipticPlotModel
        {
            get { return pEllipticPlotmodel; }
            set { SetProperty(ref pEllipticPlotmodel, value); }
        }

        private PlotModel pRectangularPistonPlotmodel = new PlotModel();
        public PlotModel RectangularPistonPlotmodel
        {
            get { return pRectangularPistonPlotmodel; }
            set { SetProperty(ref pRectangularPistonPlotmodel, value); }
        }
        private PlotModel pAreaModel = new PlotModel();
        public PlotModel AreaModel
        {
            get { return pAreaModel; }
            set { SetProperty(ref pAreaModel, value); }
        }

    }
}
