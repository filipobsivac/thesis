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
        private static List<Molecule> molecules;

        static void Main(string[] args)
        {
            string path;
                
            path = @"C:\Users\Filip\MUNI\THESIS\bc_thesis\bc_thesis\TestFiles\ElemBond.txt";
            LoadParameters(path);
            Console.WriteLine("Element parameters with specified bonds:");
            Console.WriteLine();
            Console.WriteLine("Symbol\tBond\tA\tB\tKappa");
            foreach (ElementParameters e in elemParams)
            {
                Console.WriteLine($"{e.ElementName}\t{e.BondType}\t{e.A}\t{e.B}\t{e.Kappa}");                
            }
            Console.WriteLine();

            path = @"C:\Users\Filip\MUNI\THESIS\bc_thesis\bc_thesis\TestFiles\Element.txt";
            LoadParameters(path);
            Console.WriteLine("Element parameters without specified bonds:");
            Console.WriteLine();
            Console.WriteLine("Symbol\tA\tB\tKappa");
            foreach (ElementParameters e in elemParams)
            {
                Console.WriteLine($"{e.ElementName}\t{e.A}\t{e.B}\t{e.Kappa}");
            }
            Console.WriteLine();

            path = @"C:\Users\Filip\MUNI\THESIS\bc_thesis\bc_thesis\TestFiles\set01.sdf";
            LoadMolecules(path);
            Console.WriteLine($"Number of molecules: {molecules.Count}");
            Console.WriteLine();

            Console.WriteLine("First Molecule:");
            Console.WriteLine();
            Console.WriteLine("NSC\tAtoms\tBonds");
            var m = molecules.ElementAt(0);
            Console.WriteLine($"{m.NSC}\t{m.NumOfAtoms}\t{m.NumOfBonds}");
            Console.WriteLine();

            Console.WriteLine("Atoms in molecule:");
            Console.WriteLine();
            Console.WriteLine("ID\tSymbol\tX\tY\tZ\tBonds");
            foreach(Atom a in m.Atoms)
            {
                Console.WriteLine($"{a.ID}\t{a.Symbol}\t{a.X}\t{a.Y}\t{a.Z}\t{a.Bonds.Count}");
            }
            Console.WriteLine();

            Console.WriteLine("First atom's bonds:");
            Console.WriteLine();
            Console.WriteLine("1st\t2nd\tType");
            var atom = m.Atoms.ElementAt(0);
            var bonds = atom.Bonds;
            foreach(KeyValuePair<int, int> b in bonds)
            {
                Console.WriteLine($"{atom.ID}\t{b.Key}\t{b.Value}"); 
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

        private static void LoadMolecules(string fileName)
        {
            molecules = new List<Molecule>();
            StreamReader reader = File.OpenText(fileName);
            string line;
            int lineNum = 1;
            Molecule molecule = new Molecule();

            while ((line = reader.ReadLine()) != null)
            {
                if (lineNum == 1)
                {
                    lineNum++;
                    string[] items = line.Split('_');
                    molecule.NSC = int.Parse(items[1]);
                }
                else if(lineNum <= 3)
                {
                    lineNum++;
                }
                else if(lineNum == 4)
                {
                    molecule.NumOfAtoms = int.Parse(line.Substring(0, 3));
                    molecule.NumOfBonds = int.Parse(line.Substring(3, 3));
                    lineNum++;
                }
                else if(lineNum <= molecule.NumOfAtoms + 4)
                {
                    Atom atom = new Atom(
                        lineNum - 4,
                        line.Substring(31, 3),
                        double.Parse(line.Substring(0, 10), CultureInfo.InvariantCulture),
                        double.Parse(line.Substring(10, 10), CultureInfo.InvariantCulture),
                        double.Parse(line.Substring(20, 10), CultureInfo.InvariantCulture));
                    molecule.Atoms.Add(atom);
                    lineNum++;
                }   
                else if (lineNum <= molecule.NumOfAtoms + molecule.NumOfBonds + 4){
                    var firstAtom = molecule.Atoms.Find(a => a.ID == int.Parse(line.Substring(0, 3)));
                    firstAtom.Bonds.Add(
                        int.Parse(line.Substring(3, 3)), 
                        int.Parse(line.Substring(6, 3)));
                    lineNum++;
                }else if (line.Equals("$$$$"))
                {
                    lineNum = 1;
                    molecules.Add(molecule);
                    molecule = new Molecule();
                }else
                {
                    lineNum++;
                }
            }
        }
    }
}
