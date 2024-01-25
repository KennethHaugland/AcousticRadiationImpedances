using OxyPlot.Series;
using OxyPlot;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Text;
using System.Threading.Tasks;
using OxyPlot.Axes;
using OxyPlot.Legends;

namespace RadiationOfSound.ViewModels
{
    public static class OxyPlotHelper
    {

        public static PlotModel PlotAreaGraph(List<double> MaxArray, List<double> MinArray, List<double> MaxArrayI, List<double> MinArrayI, string degree)
        {
            PlotModel PlotModel = new PlotModel();

            var areaSeriesReal = new AreaSeries()
            {
                StrokeThickness = 2,
                LineStyle = LineStyle.Solid,
                Color = OxyColors.Blue,
                Color2 = OxyColors.Blue,
                Fill = OxyColor.FromRgb(214, 231, 242),
                DataFieldX2 = "X",
                ConstantY2 = 0
            };


            var areaSeriesImg = new AreaSeries()
            {
                StrokeThickness = 2,
                LineStyle = LineStyle.Solid,
                Color = OxyColors.Red,
                Color2 = OxyColors.Red,
                Fill = OxyColor.FromRgb(214, 231, 242),
                DataFieldX2 = "X",
                ConstantY2 = 0
            };

            PlotModel.Axes.Add(new LinearAxis()
            {Position = AxisPosition.Left,
             Title= "Z/ρ₀c₀S",
            });

            PlotModel.Axes.Add(new LinearAxis()
            {
                Position = AxisPosition.Bottom,
                Title = "ka",
            });

            double ind = 0.2;

            for (int i = 0; i < MaxArray.Count; i++)
            {
                areaSeriesReal.Points.Add(new DataPoint(ind * i, MaxArray[i]));
                areaSeriesReal.Points2.Add(new DataPoint(ind * i, MinArray[i]));

                areaSeriesImg.Points.Add(new DataPoint(ind * i, MaxArrayI[i]));
                areaSeriesImg.Points2.Add(new DataPoint(ind * i, MinArrayI[i]));

                ind += 0.2;
            }

            PlotModel.Series.Add(areaSeriesReal);
            PlotModel.Series.Add(areaSeriesImg);

            return PlotModel;

        }

        public static PlotModel PlotArray(List<List<Complex>> array, string degree)
        {

            PlotModel plotModel2 = new PlotModel();

            plotModel2.Legends.Add(new Legend()
            {
                LegendTitle = "Legend",
                LegendPosition = LegendPosition.TopLeft,
            });


            string[] desc = new string[3] { "Narrow strip","Wide strip","Rectangular normal"};

            for (int i = 0; i < array.Count; i++)
            {
                var item2 = array[i];

                

                LineSeries RealPart = new LineSeries();
                RealPart.Title = "Real part - " + degree + " Type: " + desc[i];
                RealPart.MarkerType = MarkerType.None;

                LineSeries ImaginaryPart = new LineSeries();
                ImaginaryPart.Title = "Imaginary part - " + degree + " Type: " + desc[i];
                ImaginaryPart.MarkerType = MarkerType.None;
                ImaginaryPart.LineStyle = LineStyle.Dash;
                double ind = 0;
                foreach (var item in item2)
                {
                    RealPart.Points.Add(new DataPoint(ind, (double)item.Real));
                    ImaginaryPart.Points.Add(new DataPoint(ind, (double)item.Imaginary));
                    ind += 0.2;
                }
                plotModel2.Series.Add(RealPart);
                plotModel2.Series.Add(ImaginaryPart);

            }

            // the y-axis
            plotModel2.Axes.Add(new LinearAxis
            {
                Position = AxisPosition.Left,
                Minimum = 0,
                Maximum = 2,
                Title = "Z/ρ₀c₀S"
            });

            plotModel2.Axes.Add(new LinearAxis()
            {
                Position = AxisPosition.Bottom,
                Title = "ka",
            });

            return plotModel2;
        }

        public static PlotModel PlotArray(List<Complex> array, string degree)
        {

            PlotModel plotModel2 = new PlotModel();

            LineSeries RealPart = new LineSeries();
            RealPart.Title = "Real part - " + degree;
            RealPart.MarkerType = MarkerType.None;

            LineSeries ImaginaryPart = new LineSeries();
            ImaginaryPart.Title = "Imaginary part - " + degree;
            ImaginaryPart.MarkerType = MarkerType.None;
            ImaginaryPart.LineStyle = LineStyle.Dash;

            double ind = 0;
            foreach (var item in array)
            {
                RealPart.Points.Add(new DataPoint(ind, (double)item.Real));
                ImaginaryPart.Points.Add(new DataPoint(ind, (double)item.Imaginary));
                ind += 0.2;
            }

            plotModel2.Series.Add(RealPart);
            plotModel2.Series.Add(ImaginaryPart);
            // the y-axis
            plotModel2.Axes.Add(new OxyPlot.Axes.LinearAxis
            {
                Position = OxyPlot.Axes.AxisPosition.Left,
                Minimum = 0,
                Maximum = 2,
                Title = "Z/ρ₀c₀S"});

            plotModel2.Axes.Add(new LinearAxis()
            {
                Position = AxisPosition.Bottom,
                Title = "ka",
            });

            return plotModel2;
        }
    }
}
