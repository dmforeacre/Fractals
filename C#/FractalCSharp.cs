/*
    Filename    : fractal.cs
    Author      : Daniel Foreacre
    Course      : CSCI 476-01
    Date        : 11/02/21
    Assignment  : Final project/presentation
    Description : Program to create an image of Newton's fractal using parallel processing to compute
                  each iteration of Newton's method.

*/

/*****************************************************************************/
// System includes

using System;
using System.Numerics;
using System.Collections.Generic;
using System.Windows.Forms;
using System.Drawing;
using System.Drawing.Imaging;
using System.Diagnostics;
using System.Threading.Tasks;

/*****************************************************************************/

namespace Fractal
{
    /*************************************************************************/
    // Class to represent a polynomial as a list of terms
    public class Polynomial
    {
        const double TOLERANCE = .001;
        const int MAX_ITER = 200;
        const int A = 1;

        // Structure to represent a single term
        private struct Term {
            public Complex coef;
            public Complex exp;

            public Term(Complex c, Complex e)
            {
                coef = c;
                exp = e;
            }
        }

        // Structure to represent a single expression for FOILing
        private struct Expr {
            public List<Term> expr;

            public Expr(List<Term> e)
            {
                expr = e;
            }
        }

        private List<Term> polynomial = new List<Term>();
        private int size;

        // Constructor initializes the list and sets size to 0
        public Polynomial(Complex[,] terms)
        {
            this.size = 0;
            for (int i = 0; i < terms.Length/2; ++i)
            {
                this.polynomial.Add(new Term() {coef = terms[i,0], exp = terms[i,1]});
                ++this.size;
            }
        }

        /* Adds a term to the polynomial
           @param       c       Coefficient of the term as a double
           @param       e       Exponent of the term as a double */
        public void AddTerm(Complex c, Complex e)
        {
            this.polynomial.Add(new Term() { coef = c, exp = e});
            ++this.size;
        }

        /* Builds a polynomial from a given list of roots and their multiplicities
           @param       roots       The roots to use in building a polynomial function as a List<Complex>
           @param       mult        The corresponding multiplicities of each root in a parallel List<int> */
        public void BuildPolynomial(List<Complex> roots, List<int> mult)
        {
            List<Expr> expressions = new List<Expr>();
            for (int i = 0; i < roots.Count; ++i)
            {
                //Expr e = new Expr(new Term(1,1), new Term(-(roots[i].Real),0), new Term(-(roots[i].Imaginary),0), mult[i]);
                for (int j = 0; j < mult[i]; ++j)
                {
                    Expr e = new Expr(new List<Term>() {new Term(1,1), new Term(-(roots[i].Real),0), new Term(-(roots[i].Imaginary),0)});
                    expressions.Add(e);
                }
            }
            List<Expr> expansion = new List<Expr>();
            for (int i = 0; i < expressions.Count; ++i)
            {                
                ++this.size;
            }
        }

        /* Gets the number of terms in the polynomial
           @return      Size of the polynomial list as an int */
        public int GetSize()
        {
            return this.size;
        }

        /* Gets the coefficient of the term at the given index
           @param       t       Index to return the term of
           @return      The coefficient of the term at this index */
        public Complex GetCoef(int t)
        {
            return this.polynomial[t].coef;
        }

        /* Gets the coefficient of the term at the given index
           @param       t       Index to return the term of
           @return      The coefficient of the term at this index */
        public Complex GetExp(int t)
        {
            return this.polynomial[t].exp;
        }

        /* Prints the given polynomial to the console */
        public override string ToString()
        {
            string output = "";
            if (this.size == 0)
            {
                output = "Polynomial is empty!";
                return output;
            }
            for (int i = 0; i < this.size; ++i)
            {
                if (this.polynomial[i].exp != 0)
                    output += "(" + this.polynomial[i].coef + "x^" + this.polynomial[i].exp + ")";
                else
                    output += this.polynomial[i].coef;
                if (i != this.size - 1)
                    output += " + ";
            }
            return output;
        }

        /* Evaluates the polynomial at the given point
           @param       point       Point to evaluate as a complex number
           @return      The evaluation of the polynomial at the given point as a complex number */
        public Complex Evaluate(Complex point)
        {
            Complex n = new Complex();
            for (int i = 0; i < this.size; ++i)
                n += this.polynomial[i].coef * Complex.Pow(point, this.polynomial[i].exp);
            return n;
        }

        /* Calculates the first order derivative of a given polynomial 
        @param       poly        Polynomial to take the derivative of as a vector of terms
        @return      A vector of terms representing the first order derivative of the polynomial */
        public Polynomial GetDerivative()
        {
            Complex[,] terms = new Complex[this.GetSize()-1,2];
            for (int i = 0; i < this.GetSize() - 1; ++i)
            {
                terms.SetValue(this.GetCoef(i) * this.GetExp(i),i,0);
                terms.SetValue(this.GetExp(i) - 1,i,1);
            }
            Polynomial deriv = new Polynomial(terms);
            return deriv;
        }        

        /* Recursive function to approximate the closest root of a polynomial, to within the constant
           tolerance value, starting at the given point
           @param       poly        Polynomial that is having its root approximated
           @param       deriv       The derivative of the given polynomial
           @param       point       A complex number representing the current point to be checking
           @return      The closest root to the original point, to within TOLERANCE, as a complex number*/
        public Complex Newton(Polynomial deriv, Complex point, ref int iter)
        {
            if (deriv.Evaluate(point) == 0)
                throw new DivideByZeroException("Division by zero");
            if (iter > MAX_ITER)
                throw new TimeoutException("Max iterations reached");
            Complex newPoint = point - ((this.Evaluate(point) / deriv.Evaluate(point)));
            Complex root;
            if (Complex.Abs(newPoint - point) < TOLERANCE)
                root = newPoint;
            else
            {
                ++iter;
                root = this.Newton(deriv, newPoint, ref iter);
            }
            return root;
        }

        public static Complex Norm(Complex c)
        {
            return new Complex(Math.Round(c.Real*100)/100, Math.Round(c.Imaginary*100)/100);
        }
    }

    /***********************************************************************/
    // Class to package the various parameters needed for thread creation

    public class ThreadParams
    {
        public static Complex[] roots;
        public static Polynomial poly, deriv;
        public static int[,] pixels;
        public static int[,] iterations;
        public static Complex offset;
        public static double scale;

        // Id for thread
        public int id;
    }

    /************************************************************************/
    // Main class, also creates the Windows Form for display

    public class Program : Form
    {
        /********************************************************************/
        // Constants

        const int WIN_X_SIZE = 800,
                  WIN_Y_SIZE = 600;
        const int MAX_THREADS = 8;
        Complex NO_ROOT = new Complex(0,0);
        readonly int[,] COLORS = {{140,95,102},{172,188,165},{232,185,171},{224,152,145},{203,118,158},{44, 17, 30},{49, 33, 35},{78, 94, 69}};
                //{Color.FromArgb(172,188,165),Color.FromArgb(232,185,171),Color.FromArgb(224,152,145),Color.FromArgb(203,118,158),Color.FromArgb(140,95,102)};
                //{};

        /********************************************************************/
        // Class Functions

        // Function to set up the Form
        public void FormLayout()
        {
            this.Name = "Form1";
            this.Text = "Newton's Fractal";
            this.Size = new System.Drawing.Size(Math.Min(WIN_X_SIZE, 1600),Math.Min(WIN_Y_SIZE,1000));
            this.StartPosition = FormStartPosition.CenterScreen;
            //this.MouseWheel += new MouseEventHandler(this.MainForm_MouseWheel);
            this.Paint += new PaintEventHandler(this.MainForm_Paint);
        }

        // Creates reference to current instance of program
        public static Program p = new Program();
        [STAThread]

        /********************************************************************************/

        static void Main(string[] args)
        {   
            //Polynomial poly2 = new Polynomial();
            //poly2.BuildPolynomial(new List<Complex> {new Complex(0,0),new Complex(-1,0),new Complex(0,1),new Complex(0,-1)},
            //                      new List<int> {3, 2, 1, 1});
            p.FormLayout();;
            Application.Run(p);
        }

        /* Function to paint to Windows Form
           @param       sender      Reference to source of call
           @param       e           Additional arguments for paint */
        private void MainForm_Paint(object sender, PaintEventArgs e)
        {
            // Set up polynomial
            // Working polys: 1,5  -6,4  8,3  -4,2  -20,1  16,0
            //                1,5   0,4  0,3   1,2   -1,1   1,0
            //                           1,3  -4,2    1,1  -4,0
            //                           1,3   0,2    0,1  -1,0
            Complex[,] terms = {{new Complex(1,0),new Complex(3,0)},{new Complex(-1,0),new Complex(0,0)}};
            Polynomial poly = new Polynomial(terms);

            // Initialize thread parameters
            Graphics graphics = CreateGraphics();
            Bitmap bitmap = new Bitmap(WIN_X_SIZE, WIN_Y_SIZE);
            ThreadParams.roots = new Complex[]{Polynomial.Norm(new Complex(1,0)),Polynomial.Norm(new Complex(-.5,Math.Sqrt(3)/2)),Polynomial.Norm(new Complex(-.5,-Math.Sqrt(3)/2))};
            ThreadParams.pixels = new int[WIN_X_SIZE, WIN_Y_SIZE];
            ThreadParams.iterations = new int[WIN_X_SIZE, WIN_Y_SIZE];
            ThreadParams.offset = new Complex(0,0);
            ThreadParams.scale = 1;
            ThreadParams.poly = poly;
            ThreadParams.deriv = poly.GetDerivative();

            string type = "Parallel";
            Stopwatch watch = new Stopwatch();
            watch.Start();

            if (type == "Parallel")
            {
                Task[] tasks = new Task[MAX_THREADS];
                // Thread creation loop
                for (int t = 0; t < MAX_THREADS; ++t)
                {                    
                    ThreadParams args = new ThreadParams();
                    args.id = t;
                    Task task = new Task(GetPointsParallel, args);
                    tasks[t] = task;
                    task.Start();                
                }
                Task.WaitAll(tasks);
            }
            else if (type == "Serial")
                GetPointsSerial(); 
            
            watch.Stop();

            // Create bitmap based on calculated roots
            for (int i = 0; i < WIN_X_SIZE; ++i)
                for (int j = 0; j < WIN_Y_SIZE; ++j)
                {
                    int R = Math.Clamp(COLORS[ThreadParams.pixels[i,j],0] - 3*ThreadParams.iterations[i,j],0,255);
                    int G = Math.Clamp(COLORS[ThreadParams.pixels[i,j],1] - 3*ThreadParams.iterations[i,j],0,255);
                    int B = Math.Clamp(COLORS[ThreadParams.pixels[i,j],2] - 3*ThreadParams.iterations[i,j],0,255);
                    bitmap.SetPixel(i, j, Color.FromArgb(R,G,B)); 
                }

            if (WIN_X_SIZE > 1600)
            {
                string fileName = "Fractal" + WIN_X_SIZE + "x" + WIN_Y_SIZE + ".png";
                bitmap.Save(@fileName, ImageFormat.Png);
                MessageBox.Show("Finished in " + watch.ElapsedMilliseconds + " ms\n" + "Bitmap saved to " + fileName);
            }
            else
            {
                graphics.DrawImage(bitmap, 0, 0);
            }

            string output = type + " " + WIN_X_SIZE + "x" + WIN_Y_SIZE + " in C#" +
                            "\nUsing " + MAX_THREADS + "/" + Environment.ProcessorCount + " Threads" + 
                            "\nScale: " + ThreadParams.scale + " Offset: " + ThreadParams.offset.Real + "," + ThreadParams.offset.Imaginary +
                            "\nFunction: " + ThreadParams.poly +
                            "\nIn " + watch.ElapsedMilliseconds + " ms";

            FontFamily myFontFamily = new FontFamily("Arial");
            Font myFont = new Font(myFontFamily, 12);
            graphics.DrawString(output, myFont, Brushes.AntiqueWhite, new PointF(5,5));

            myFontFamily.Dispose();
            myFont.Dispose();
            bitmap.Dispose();
            graphics.Dispose();
        }

        private void MainForm_MouseWheel(object sender, MouseEventArgs e)
        {
            ThreadParams.scale += e.Delta;
        }

        /* Calculates the root (via Newton's method) of each pixel on the screen in parallel
           @param       obj         An object containing necessary ThreadParams for calculation */
        public void GetPointsParallel(Object obj)
        {
            ThreadParams args = (ThreadParams)obj;
            int yStart = (args.id*WIN_Y_SIZE)/MAX_THREADS, 
                yEnd = ((args.id+1)*WIN_Y_SIZE)/MAX_THREADS;

            for (int x = 0; x < WIN_X_SIZE; ++x)
            {
                for (int y = yStart; y < yEnd; ++y)
                {
                    Complex point = new Complex((-(.5 * WIN_X_SIZE)+x+ThreadParams.offset.Real)*ThreadParams.scale, (-(.5 * WIN_Y_SIZE)+y+ThreadParams.offset.Imaginary)*ThreadParams.scale);
                    Complex root;
                    int iter = 0;
                    try                       
                    {
                        root = ThreadParams.poly.Newton(ThreadParams.deriv, point, ref iter);
                    }
                    catch (System.Exception)
                    {
                        root = ThreadParams.roots[0];
                    }
                    ThreadParams.iterations[x,y] = iter;
                    ThreadParams.pixels[x,y] = Array.IndexOf(ThreadParams.roots,Polynomial.Norm(root));
                }                   
            }
        }

        /* Calculates the root (via Newton's method) of each pixel on the screen serially
           @param       obj         An object containing necessary ThreadParams for calculation */
        public void GetPointsSerial()
        {
            for (int x = 0; x < WIN_X_SIZE; ++x)
            {
                for (int y = 0; y < WIN_Y_SIZE; ++y)
                {
                    Complex point = new Complex((-(.5 * WIN_X_SIZE)+x+ThreadParams.offset.Real)*ThreadParams.scale, (-(.5 * WIN_Y_SIZE)+y+ThreadParams.offset.Imaginary)*ThreadParams.scale);
                    Complex root;
                    int iter = 0;
                    try                       
                    {
                        root = ThreadParams.poly.Newton(ThreadParams.deriv, point, ref iter);
                    }
                    catch (System.Exception)
                    {                  
                        root = ThreadParams.roots[0];
                    }
                    ThreadParams.iterations[x,y] = iter;
                    ThreadParams.pixels[x,y] = Array.IndexOf(ThreadParams.roots,Polynomial.Norm(root));
                }                   
            }
        }
    }
}
