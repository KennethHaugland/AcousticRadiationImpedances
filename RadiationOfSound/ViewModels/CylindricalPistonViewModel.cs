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
    public class CylindricalPistonViewModel:BindableBase
    {
        private string _Title = "Cylindrical radiator";
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

        public CylindricalPistonViewModel()
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

            List<Complex> Piston = new List<Complex>();

            for (double k = K_Low; k <= K_High; k += deltaK)
            {
                Piston.Add(Acoustics.RadiationImpedance.Objects.Cylindrical(k, a, 6, 0d));
            }

            Plotter = OxyPlotHelper.PlotArray(Piston, "");


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
