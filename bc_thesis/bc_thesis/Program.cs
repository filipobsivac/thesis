using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;
using thesis.Classes;

namespace bc_thesis
{
    class Program
    {
        private static List<ElementParameters> elemParams;

        static void Main(string[] args)
        {
            string path = @"C:\Users\Filip\MUNI\THESIS\bc_thesis\bc_thesis\TestFiles\ElemBond.txt";
            LoadParameters(path);
            Console.WriteLine("Elements with specified bonds:");
            Console.WriteLine("--------------------------------------");
            Console.WriteLine("Symbol\tBond\tA\tB\tKappa");
            foreach (ElementParameters e in elemParams)
            {
                Console.WriteLine($"{e.ElementName}\t{e.BondType}\t{e.A}\t{e.B}\t{e.Kappa}");                
            }
            path = @"C:\Users\Filip\MUNI\THESIS\bc_thesis\bc_thesis\TestFiles\Element.txt";
            LoadParameters(path);
            Console.WriteLine();
            Console.WriteLine("Elements without specified bonds:");
            Console.WriteLine("------------------------------");
            Console.WriteLine("Symbol\tA\tB\tKappa");
            foreach (ElementParameters e in elemParams)
            {
                Console.WriteLine($"{e.ElementName}\t{e.A}\t{e.B}\t{e.Kappa}");
            }
            Console.ReadKey();
        }

        private static void LoadParameters(string fileName)
        {
            elemParams = new List<ElementParameters>();
            StreamReader reader = File.OpenText(fileName);
            string line;
            bool firstLine = true;
            bool bonds = true;
            double kappa = 0;

            while ((line = reader.ReadLine()) != null)
            {
                if (firstLine)
                {
                    firstLine = false;
                    string[] items = line.Split(' ');
                    kappa = double.Parse(items[1], CultureInfo.InvariantCulture);
                    if (items[0].Equals("Element"))
                    {
                        bonds = false;
                    }
                    continue;
                }

                string[] param = line.Split(' ');
                ElementParameters parameters;
                if (bonds)
                {
                    parameters = new ElementParameters(
                        param[0], 
                        double.Parse(param[2], CultureInfo.InvariantCulture), 
                        double.Parse(param[3], CultureInfo.InvariantCulture), 
                        kappa, 
                        int.Parse(param[1]));
                }
                else
                {
                    parameters = new ElementParameters(
                        param[0], 
                        double.Parse(param[1], CultureInfo.InvariantCulture), 
                        double.Parse(param[2], CultureInfo.InvariantCulture), 
                        kappa);
                }

                elemParams.Add(parameters);
            }
        }
    }
}
