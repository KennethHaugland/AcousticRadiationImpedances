using OxyPlot;
using Prism.Commands;
using Prism.Mvvm;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;

namespace RadiationOfSound.ViewModels
{
    public class CircularPistonsViewModel: BindableBase
    {
        private string _Title = "Circular baffles";
        public string Title
        {
            get { return _Title; }
            set { SetProperty(ref _Title, value); }
        }

        private double pA = 1d;
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


        private double pK_low = 0;
        public double K_Low
        {
            get { return pK_low; }
            set { SetProperty(ref pK_low, value); }
        }
        private double K_hi = 10;
        public double K_High
        {
            get { return K_hi; }
            set { SetProperty(ref K_hi, value); }
        }

        private double _deltaK = 0.2;
        public double deltaK
        {
            get { return _deltaK; }
            set { SetProperty(ref _deltaK, value); }
        }

        public CircularPistonsViewModel()
        {
            ExecuteCalculateCmd();
        }

        private DelegateCommand _ExecuteCalculateCmd;
        public DelegateCommand CalculateCmd =>
            _ExecuteCalculateCmd ?? (_ExecuteCalculateCmd = new DelegateCommand(ExecuteCalculateCmd));

        void ExecuteCalculateCmd()
        {
            double theta = Theta * Math.PI / 180;
            double phi = Phi * Math.PI / 180;
            double a = A;
            double b = B;


            List<Complex> CirclularPiston = new List<Complex>();
            List<Complex> CircularClamped = new List<Complex>();
            List<Complex> CircularSupported = new List<Complex>();

            string[] description = { "Piston", "Clamped", "Supported" };

            for (double k = K_Low; k <= K_High; k += deltaK)
            {
                CirclularPiston.Add(Acoustics.RadiationImpedance.Pistons.CircularBaffle(k, a));
                CircularClamped.Add( Acoustics.RadiationImpedance.Pistons.CircularClampedBaffle(k, a));
                CircularSupported.Add(Acoustics.RadiationImpedance.Pistons.CircularSupportedBaffle(k, a));
            }

            Plotter = OxyPlotHelper.PlotArray(new List<List<Complex>>() { CirclularPiston, CircularClamped, CircularSupported }, "", 2.5, 0, description); 
       

            RaisePropertyChanged();
        }

        private PlotModel pPlotModel = new PlotModel();
        public PlotModel Plotter
        {
            get { return pPlotModel; }
            set { SetProperty(ref pPlotModel, value); }
        }

    }
}
