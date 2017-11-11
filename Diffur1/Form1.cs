using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace Diffur1
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }
        public static double f1(double x)
        {
            return 2.0 * x / 11;
        }
        public static double f2(double x)
        {
            return -3.0 * x / 13;
        }
        public static List<double> dudx(double t0, double tmax, int n, int eps)
        {
            double A = 2.0 / 11, B = 3.0 / 13, x0 = A * Math.PI, y0 = B * Math.PI, c = 1.0 / 9;
            double h = (tmax - t0) / (n - 1);

            List<double> y1 = new List<double>();
            List<double> y2 = new List<double>();
            y1.Add(y0);
            y2.Add(x0);
            for (int i = 0; i < n; i++)
            {
                List<double> localpogr = new List<double>();
                double a = c, b2 = 1.0 / (2 * c), b1 = 1 - b2;
                double k11 = h * f1(y2[i]);
                double k21 = h * f2(y1[i]);
                double k12 = h * f1(y2[i] + a * k21);
                double k22 = h * f2(y1[i] + a * k11);

                double y1plus1 = y1[i] + b1 * k11 + b2 * k12;
                double y2plus1 = y2[i] + b1 * k21 + b2 * k22;
                double y1sr=y1[i],y2sr=y2[i];
                if(eps!=0)
                {
                    double tempy1=y1sr;
                    double tempy2=y2sr;    
                    for(int j=0;j<2;j++)
                    {
                        
                        double tempk11 = h * f1(tempy2);
                        double tempk21 = h * f2(tempy1);
                        double tempk12 = h * f1((tempy2) + a * tempk21);
                        double tempk22 = h * f2((tempy1) + a * tempk11);
                        
                        y1sr= tempy1 + b1 * tempk11 + b2 * tempk12;
                        y2sr = tempy2 + b1 * tempk21 + b2 * tempk22;
                    }
                }
                
                localpogr.Add((y1sr - y1plus1)/(1-Math.Pow(2,-2)));
                localpogr.Add((y2sr - y2plus1) / (1 - Math.Pow(2, -2)));
                y1.Add(y1plus1);
                y2.Add(y2plus1);

            }
            List<double> ret = new List<double>();
            ret.Add(y1[y1.Count - 1]);
            ret.Add(y2[y2.Count - 1]);
            return ret;
        }
        public static bool udovlet(double y1, double y2,double y1plus1, double y2plus1,double h,double eps)
        {
            double c = 1.0 / 9, a = c, b2 = 1.0 / (2 * c), b1 = 1 - b2;
            double tempy1 = y1;
            double tempy2 = y2;
            double hloc = h / 2;
            for (int j = 0; j < 2; j++)
            {
                double tempk11 = hloc * f1(tempy2);
                double tempk21 = hloc * f2(tempy1);
                double tempk12 = hloc * f1((tempy2) + a * tempk21);
                double tempk22 = hloc * f2((tempy1) + a * tempk11);

                tempy1 = tempy1 + b1 * tempk11 + b2 * tempk12;
                tempy2 = tempy2 + b1 * tempk21 + b2 * tempk22;
            }
            double abc1 = (Math.Abs(tempy1 - y1plus1) / (1 - Math.Pow(2, -2))); double abc2 = (Math.Abs(tempy2 - y2plus1) / (1 - Math.Pow(2, -2)));
            if ((Math.Abs(tempy1 - y1plus1) / (1 - Math.Pow(2, -2)) < eps) && (Math.Abs(tempy2 - y2plus1) / (1 - Math.Pow(2, -2)) < eps))
                return true;
            else
               return false;
        }
        public static double hepsloc(double y1, double y2,double y1plus1, double y2plus1,double h, double eps)
        {
            double c = 1.0 / 9, a = c, b2 = 1.0 / (2 * c), b1 = 1 - b2;
            double htemp=h;
            if (udovlet(y1, y2,y1plus1,y2plus1, htemp,eps))
                return htemp;
            else
            {
                htemp /= 2;
                double k11 = htemp * f1(y2);
                double k21 = htemp * f2(y1);
                double k12 = htemp * f1(y2 + a * k21);
                double k22 = htemp * f2(y1 + a * k11);
                double y1temp = y1 + b1 * k11 + b2 * k12;
                double y2temp = y2 + b1 * k21 + b2 * k22;
                return hepsloc(y1, y2, y1temp,y2temp, htemp,eps);
            }
        }
        public static List<double> move(double y1,double y2, double h)
        {
            double c = 1.0 / 9, a = c, b2 = 1.0 / (2 * c), b1 = 1 - b2;
            double hloc = h/2;
            double tempy1 = y1;
            double tempy2 = y2;
            for (int j = 0; j < 2; j++)
            {
                double tempk11 = hloc * f1(tempy2);
                double tempk21 = hloc * f2(tempy1);
                double tempk12 = hloc * f1((tempy2) + a * tempk21);
                double tempk22 = hloc * f2((tempy1) + a * tempk11);

                tempy1 = tempy1 + b1 * tempk11 + b2 * tempk12;
                tempy2 = tempy2 + b1 * tempk21 + b2 * tempk22;
            }
            List<double> ret = new List<double>();
            ret.Add(tempy1);
            ret.Add(tempy2);
            return ret;
        }
        public List<double> autoshag(double t0, double tmax, double eps)
        {
            Font f = new Font(FontFamily.GenericSerif, 8);
            Bitmap bitmap = new Bitmap(pictureBox1.Width, pictureBox1.Height);
            Graphics im = Graphics.FromImage(bitmap);
            pictureBox1.Image = bitmap;
            float hw = Convert.ToInt32(pictureBox1.Width) / 2;
            float hh = Convert.ToInt32(pictureBox1.Height) / 2;
            Pen pen = new Pen(Color.Black);
            pen.Width = 2;

            List<koor> graph = new List<koor>();
            im.DrawLine(pen, 0, hw, hh * 2, hw);
            im.DrawLine(pen, hh, 0, hh, hw * 2);
            double rangex = Math.PI;
            double A = 2.0 / 11, B = 3.0 / 13, x0 = A * Math.PI, y0 = B * Math.PI, c = 1.0 / 9;

            List<double> y1 = new List<double>();
            List<double> y2 = new List<double>();
            y1.Add(y0);
            y2.Add(x0);
            double qwe = t0;
            double hloc = (tmax - t0);
            int i = 0;
            while(qwe<tmax)
            {
                double a = c, b2 = 1.0 / (2 * c), b1 = 1 - b2;
                bool rdy = false;
                double y1plus1=0;
                double y2plus1=0;
                while (!rdy)
                {

                    double k11 = hloc * f1(y2[i]);
                    double k21 = hloc * f2(y1[i]);
                    double k12 = hloc * f1(y2[i] + a * k21);
                    double k22 = hloc * f2(y1[i] + a * k11);
                    y1plus1 = y1[i] + b1 * k11 + b2 * k12;
                    y2plus1 = y2[i] + b1 * k21 + b2 * k22;
                    var temp = move(y1[i], y2[i], hloc);
                    var y1loc = (Math.Abs(temp[0] - y1plus1) / (1 - Math.Pow(2, -2)));
                    var y2loc = (Math.Abs(temp[1] - y2plus1) / (1 - Math.Pow(2, -2)));
                    if(y1loc>eps*Math.Pow(2,2) || y2loc>eps*Math.Pow(2,2))
                    {
                        hloc /= 2;
                    }
                    else
                    {
                        if ((y1loc > eps && y1loc < eps * Math.Pow(2, 2)) || (y2loc > eps && y2loc < eps * Math.Pow(2, 2)))
                        {
                            hloc /= 2;
                        }
                        else
                        {
                            if (y1loc < eps && y2loc < eps)
                                rdy = !rdy;
                        }
                    }
                    
                }
                koor point = new koor();
                point.rex = qwe;
                point.rey = hloc;
                graph.Add(point);
                y1.Add(y1plus1);
                y2.Add(y2plus1);
               


                //if(udovlet(y1[i],y2[i],y1plus1,y2plus1,hloc,eps))
                //{
                //    y1.Add(y1plus1);
                //    y2.Add(y2plus1);
                //}
                //else
                //{
                //   hloc=hepsloc(y1[i],y2[i],y1plus1,y2plus1, hloc, eps);
                //   k11 = hloc * f1(y2[i]);
                //   k21 = hloc * f2(y1[i]);
                //   k12 = hloc * f1(y2[i] + a * k21);
                //   k22 = hloc * f2(y1[i] + a * k11);

                //   y1plus1 = y1[i] + b1 * k11 + b2 * k12;
                //   y2plus1 = y2[i] + b1 * k21 + b2 * k22;
                //   y1.Add(y1plus1);
                //   y2.Add(y2plus1);
                //}

                //y1.Add(y1plus1);
                //y2.Add(y2plus1);
                i++;
                qwe += hloc;
            }
            double rangey = 0.25;
            //for (int j = 0; j < graph.Count;j++)
            //{
            //    if (graph[j].rey > rangey)
            //        rangey = graph[j].rey;
            //}
                for (int j = 0; j < 4; j++)
                {
                    im.DrawLine(pen, j * hh / 2, hw - 5, j * hh / 2, hw + 5);
                    im.DrawLine(pen, hh - 5, j * hw / 2, hh + 5, j * hw / 2);
                }
            double g = -rangex;
            for (int j = 0; j <= 4; j++)
            {
                im.DrawString(g.ToString(), f, new SolidBrush(Color.Black), j * hh / 2, hw);
                g += 2 * rangex / 4;
            }
            g = rangey;
            for (int j = 0; j <= 4; j++)
            {
                im.DrawString(g.ToString(), f, new SolidBrush(Color.Black), hh, j * hw / 2);
                g -= 2 * rangey / 4;
            }
            for (int j = 0; j < graph.Count; j++)
            {
                graph[j].imx = hh + (pictureBox1.Height / (2 * rangex)) * graph[j].rex;
                graph[j].imy = hw - (pictureBox1.Width / (2 * rangey)) * graph[j].rey;
            }
            pen.Color = Color.Red;
            for (int j = 1; j < graph.Count; j++)
            {
                Point x = new Point(), y = new Point();
                x.X = Convert.ToInt32(graph[j - 1].imx);
                x.Y = Convert.ToInt32(graph[j - 1].imy);
                y.X = Convert.ToInt32(graph[j].imx);
                y.Y = Convert.ToInt32(graph[j].imy);
                im.DrawLine(pen, x, y);
            }
            List<double> ret = new List<double>();
            ret.Add(y1[y1.Count - 1]);
            ret.Add(y2[y2.Count - 1]);
            return ret;
        }
        public static List<double> globalRunge(double t0, double tmax, int n)
        {
            var y1 = dudx(t0, tmax, n,0);
            var y2 = dudx(t0, tmax, 2 * n,0);
            List<double> ret = new List<double>();
            ret.Add((y2[0] - y1[0]) / (1 - Math.Pow(2, -2)));
            ret.Add((y2[1] - y1[1]) / (1 - Math.Pow(2, -2)));
            return ret;
        }
        public static double hopt(double t0, double tmax, int n, double eps)
        {
            double ret1, ret2;
            double h = (tmax - t0) / (n - 1);
            var otv1 = dudx(t0, tmax, n,0); var otv2 = dudx(t0, tmax, 2 * n,0);
            ret1 = h * Math.Pow(((Math.Pow(2, 2) - 1) * eps) / Math.Abs((otv2[0] - otv1[0])), 1.0 / 2) / 2;
            ret2 = h * Math.Pow(((Math.Pow(2, 2) - 1) * eps) / Math.Abs((otv2[1] - otv1[1])), 1.0 / 2) / 2;
            if (ret1 > ret2)
                return ret2;
            else
                return ret1;
        }

        class koor
        {
            public double rex, rey, imx, imy;
        }
        void draw()
        {
            double t0 = 0, tmax = Math.PI, eps = 10E-6;
            int n = 1000;

            Font f = new Font(FontFamily.GenericSerif, 8);
            Bitmap bitmap = new Bitmap(pictureBox1.Width, pictureBox1.Height);
            Graphics im = Graphics.FromImage(bitmap);
            pictureBox1.Image = bitmap;
            float hw = Convert.ToInt32(pictureBox1.Width) / 2;
            float hh = Convert.ToInt32(pictureBox1.Height) / 2;
            Pen pen = new Pen(Color.Black);
            pen.Width = 2;

            Bitmap bitmap2 = new Bitmap(pictureBox2.Width, pictureBox2.Height);
            Graphics im2 = Graphics.FromImage(bitmap2);
            pictureBox2.Image = bitmap2;
            float hw2 = Convert.ToInt32(pictureBox2.Width) / 2;
            float hh2 = Convert.ToInt32(pictureBox2.Height) / 2;

            List<koor> graph = new List<koor>();
            List<koor> graph2 = new List<koor>();
            im.DrawLine(pen, 0, hw, hh * 2, hw);
            im.DrawLine(pen, hh, 0, hh, hw * 2);
            im2.DrawLine(pen, 0, hw2, hh2 * 2, hw2);
            im2.DrawLine(pen, hh2, 0, hh2, hw2 * 2);

            double rangex = Math.PI, rangey = 0.000075;
            for (int i = 0; i < 4; i++)
            {
                im.DrawLine(pen, i * hh / 2, hw - 5, i * hh / 2, hw + 5);
                im.DrawLine(pen, hh - 5, i * hw / 2, hh + 5, i * hw / 2);
                im2.DrawLine(pen, i * hh2 / 2, hw2 - 5, i * hh2 / 2, hw2 + 5);
                im2.DrawLine(pen, hh2 - 5, i * hw2 / 2, hh2 + 5, i * hw2 / 2);
            }
            double g = -rangex;
            for (int i = 0; i <= 4; i++)
            {
                im.DrawString(g.ToString(), f, new SolidBrush(Color.Black), i * hh / 2, hw);
                im2.DrawString(g.ToString(), f, new SolidBrush(Color.Black), i * hh2 / 2, hw2);
                g += 2 * rangex / 4;
            }
            g = rangey;
            for (int i = 0; i <= 4; i++)
            {
                im.DrawString(g.ToString(), f, new SolidBrush(Color.Black), hh, i * hw / 2);
                im2.DrawString(g.ToString(), f, new SolidBrush(Color.Black), hh2, i * hw2 / 2);
                g -= 2 * rangey / 4;
            }
            double h = hopt(t0, tmax, n, eps);

            int m = Convert.ToInt32(Math.Round(tmax / h))+1;
            double huse = tmax / m;
            for(int i=0;i<m;i++)
            {
                double tloc = t0 + i * huse;
                var pogr = globalRunge(t0, tloc, m);
                koor temp2 = new koor();
                temp2.rex = tloc;
                temp2.rey = pogr[1];
                graph2.Add(temp2);

                koor temp = new koor();
                temp.rex = tloc;
                temp.rey = pogr[0];
                graph.Add(temp);
            }
            pen.Color = Color.Coral;
            for (int i = 0; i < graph.Count; i++)
            {
                graph[i].imx = hh + (pictureBox1.Height / (2 * rangex)) * graph[i].rex;
                graph[i].imy = hw - (pictureBox1.Width / (2 * rangey)) * graph[i].rey;
                graph2[i].imx = hh2 + (pictureBox2.Height / (2 * rangex)) * graph2[i].rex;
                graph2[i].imy = hw2 - (pictureBox2.Width / (2 * rangey)) * graph2[i].rey;
            }
            for (int i = 1; i < graph.Count; i++)
            {
                Point x = new Point(), y = new Point();
                x.X = Convert.ToInt32(graph[i - 1].imx);
                x.Y = Convert.ToInt32(graph[i - 1].imy);
                y.X = Convert.ToInt32(graph[i].imx);
                y.Y = Convert.ToInt32(graph[i].imy);
                im.DrawLine(pen, x, y);

            }
              for (int i = 1; i < graph.Count; i++)
              {
                  pen.Color = Color.Red;
                  Point x2 = new Point(), y2 = new Point();
                  x2.X = Convert.ToInt32(graph2[i - 1].imx);
                  x2.Y = Convert.ToInt32(graph2[i - 1].imy);
                  y2.X = Convert.ToInt32(graph2[i].imx);
                  y2.Y = Convert.ToInt32(graph2[i].imy);
                  im2.DrawLine(pen, x2, y2);
              }
        }
        private void button1_Click(object sender, EventArgs e)
        {
            draw();
           //var x= autoshag(0, Math.PI, 1E-6);
           double qwe = 123;
        }
    }
}
