using MathNet.Numerics.LinearAlgebra;
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
            Console.WriteLine("sdfFileName paramsFileName methodName outputFileName");
            string arguments = Console.ReadLine();
            string[] items = arguments.Split(' ');
            string sdfFilePath = Path.Combine(Environment.CurrentDirectory, items[0]);
            LoadMolecules(sdfFilePath);            
            string paramsFilePath = Path.Combine(Environment.CurrentDirectory, items[1]);
            LoadParameters(paramsFilePath);

            string outputFilePath = Path.Combine(Environment.CurrentDirectory, items[2]);
            SolveEEM(outputFilePath);
            
            #region EEM Matrix solution test
            Console.WriteLine("EEM Matrix:\n");
            var m = molecules.ElementAt(0);
            double[,] matrix = BuildEEMMatrix(m);
            for (int i = 0; i <= m.NumOfAtoms; i++)
            {
                for (int j = 0; j <= m.NumOfAtoms; j++)
                {
                    Console.Write($"{Math.Round(matrix[i, j], 2)}\t");
                }
                Console.WriteLine();
            }
            Console.WriteLine();
            var A = Matrix<double>.Build.DenseOfArray(matrix);
            var B = Vector<double>.Build.Dense(BuildEEMVector(m));
            var X = A.Solve(B);
            int count = 0;
            Console.WriteLine("Results:\n");
            foreach(double a in X)
            {
                if(count != m.NumOfAtoms)
                {
                    if(a >= 0)
                    {
                        Console.WriteLine($"{count + 1}\t{m.Atoms[count].Symbol}\t {Math.Round(a, 6)}");
                    }else
                    {
                        Console.WriteLine($"{count + 1}\t{m.Atoms[count].Symbol}\t{Math.Round(a, 6)}");
                    }                    
                }                
                count++;
            }
            #endregion

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
                else if (lineNum <= 3)
                {
                    lineNum++;
                }
                else if (lineNum == 4)
                {
                    molecule.NumOfAtoms = int.Parse(line.Substring(0, 3));
                    molecule.NumOfBonds = int.Parse(line.Substring(3, 3));
                    lineNum++;
                }
                else if (lineNum <= molecule.NumOfAtoms + 4)
                {
                    Atom atom = new Atom(
                        lineNum - 4,
                        line.Substring(31, 3).Trim(' '),
                        double.Parse(line.Substring(0, 10), CultureInfo.InvariantCulture),
                        double.Parse(line.Substring(10, 10), CultureInfo.InvariantCulture),
                        double.Parse(line.Substring(20, 10), CultureInfo.InvariantCulture));
                    molecule.Atoms.Add(atom);
                    lineNum++;
                }
                else if (lineNum <= molecule.NumOfAtoms + molecule.NumOfBonds + 4) {
                    var firstAtom = molecule.Atoms.Find(a => a.ID == int.Parse(line.Substring(0, 3)));
                    firstAtom.Bonds.Add(
                        int.Parse(line.Substring(3, 3)),
                        int.Parse(line.Substring(6, 3)));
                    lineNum++;
                }
                else if (line.Equals("$$$$"))
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

        private static double CalculateDistance(Atom atom1, Atom atom2)
        {
            return Math.Sqrt(
                  (atom2.X - atom1.X) * (atom2.X - atom1.X) 
                + (atom2.Y - atom1.Y) * (atom2.Y - atom1.Y) 
                + (atom2.Z - atom2.Z) * (atom2.Z - atom1.Z)
                );
        }

        private static double GetElectronegativity(Atom atom)
        {
            switch (atom.Symbol)
            {
                case "H": return 2.20;
                case "Li": return 0.98;
                case "Be": return 1.57;
                case "B": return 2.04;
                case "C": return 2.55;
                case "N": return 3.04;
                case "O": return 3.44;
                case "F": return 3.98;
                case "Na": return 0.93;
                case "Mg": return 1.31;
                case "Al": return 1.61;
                case "Si": return 1.90;
                case "P": return 2.19;
                case "S": return 2.58;
                case "Cl": return 3.16;
                case "K": return 0.82;
                case "Ca": return 1.00;
                case "Sc": return 1.36;
                case "Ti": return 1.54;
                case "V": return 1.64;
                case "Cr": return 1.66;
                case "Mn": return 1.55;
                case "Fe": return 1.83;
                case "Co": return 1.88;
                case "Ni": return 1.91;
                case "Cu": return 1.90;
                case "Zn": return 1.65;
                case "Ga": return 1.81;
                case "Ge": return 2.01;
                case "As": return 2.18;
                case "Se": return 2.55;
                case "Br": return 2.96;
                case "Kr": return 3.00;
                case "Rb": return 0.82;
                case "Sr": return 0.95;
                case "Y": return 1.22;
                case "Zr": return 1.33;
                case "Nb": return 1.6;
                case "Mo": return 2.16;
                case "Tc": return 1.9;
                case "Ru": return 2.2;
                case "Rh": return 2.28;
                case "Pd": return 2.20;
                case "Ag": return 1.93;
                case "Cd": return 1.69;
                case "In": return 1.78;
                case "Sn": return 1.96;
                case "Sb": return 2.05;
                case "Te": return 2.1;
                case "I": return 2.66;
                case "Xe": return 2.60;
                case "Cs": return 0.79;
                case "Ba": return 0.89;
                case "La": return 1.1;
                case "Hf": return 1.3;
                case "Ta": return 1.5;
                case "W": return 2.36;
                case "Re": return 1.9;
                case "Os": return 2.2;
                case "Ir": return 2.20;
                case "Pt": return 2.28;
                case "Au": return 2.54;
                case "Hg": return 2.00;
                case "Tl": return 1.62;
                case "Pb": return 1.87;
                case "Bi": return 2.02;
                case "Po": return 2.0;
                case "At": return 2.2;
                case "Rn": return 2.2;
                case "Fr": return 0.7;
                case "Ra": return 0.9;
                case "Ac": return 1.1;
                case "Ce": return 1.12;
                case "Pr": return 1.13;
                case "Nd": return 1.14;
                case "Pm": return 1.13;
                case "Sm": return 1.17;
                case "Eu": return 1.2;
                case "Gd": return 1.2;
                case "Tb": return 1.1;
                case "Dy": return 1.22;
                case "Ho": return 1.23;
                case "Er": return 1.24;
                case "Tm": return 1.25;
                case "Yb": return 1.1;
                case "Lu": return 1.27;
                case "Th": return 1.3;
                case "Pa": return 1.5;
                case "U": return 1.38;
                case "Np": return 1.36;
                case "Pu": return 1.28;
                case "Am": return 1.13;
                case "Cm": return 1.28;
                case "Bk": return 1.3;
                case "Cf": return 1.3;
                case "Es": return 1.3;
                case "Fm": return 1.3;
                case "Md": return 1.3;
                case "No": return 1.3;
                case "Lr": return 1.3;
                default: return 0;
            }
        }
        
        private static double[,] BuildEEMMatrix(Molecule m)
        {
            int columns = m.NumOfAtoms + 1;
            int rows = m.NumOfAtoms + 1;
            double[,] arr = new double[rows, columns];
            var atom = m.Atoms.ElementAt(0);
            var param = elemParams.First(x => x.ElementName.Equals(atom.Symbol));
            double k = param.Kappa;

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    if (j == columns - 1)
                    {
                        arr[i, j] = -1;
                    }
                    else if (i == j)
                    {
                        var a = m.Atoms.First(x => x.ID == i + 1);
                        var pars = elemParams.FindAll(x => x.ElementName.Equals(a.Symbol));
                        double B = 0;
                        int max = 0;
                        foreach(var p in pars)
                        {
                            if(p.BondType >= max)
                            {
                                B = p.B;
                            }
                        }
                        arr[i, j] = B;
                    }
                    else if (i == rows - 1)
                    {
                        arr[i, j] = 1;
                    }
                    else
                    {
                        var a1 = m.Atoms.First(x => x.ID == i + 1);
                        var a2 = m.Atoms.First(x => x.ID == j + 1);
                        arr[i, j] = (k / (CalculateDistance(a1, a2)));
                    }
                }
            }
            arr[rows - 1, columns - 1] = 0;
            return arr;
        }

        private static double[] BuildEEMVector(Molecule m)
        {
            double[] vector = new double[m.NumOfAtoms + 1];
            int count = 0;
            foreach (Atom a in m.Atoms)
            {
                var pars = elemParams.FindAll(x => x.ElementName.Equals(a.Symbol));
                double A = 0;
                int max = 0;
                foreach (var p in pars)
                {
                    if (p.BondType >= max)
                    {
                        A = p.A;
                    }
                }
                vector[count] = -A;
                count++;
            }
            vector[count] = 0;//Q
            return vector;
        }

        private static void SolveEEM(string outputFilePath)
        {
            using (StreamWriter file = new StreamWriter(outputFilePath))
            {
                foreach(Molecule molecule in molecules)
                {
                    double[,] m = BuildEEMMatrix(molecule);
                    var matrix = Matrix<double>.Build.DenseOfArray(m);
                    var vector = Vector<double>.Build.Dense(BuildEEMVector(molecule));
                    var results = matrix.Solve(vector);
                    file.WriteLine($"NSC_{molecule.NSC}");
                    file.WriteLine(molecule.NumOfAtoms);
                    int count = 0;
                    foreach(double charge in results)
                    {
                        if(count != molecule.NumOfAtoms)
                        {
                            molecule.Atoms[count].Charge = charge;
                            if(charge >= 0)
                            {
                                file.WriteLine($"{count + 1}\t{molecule.Atoms[count].Symbol}\t {charge}");
                            }else
                            {
                                file.WriteLine($"{count + 1}\t{molecule.Atoms[count].Symbol}\t{charge}");
                            }                            
                            count++;
                        }
                    }
                    file.WriteLine("$$$$");
                }                
            }
        }

        private static double[,] BuildDegreeMatrix(Molecule m)
        {
            double[,] arr = new double[m.NumOfAtoms, m.NumOfAtoms];

            for (int i = 0; i < m.Atoms.Count; i++)
            {
                int rank = 0;
                foreach (var b in m.Atoms[i].Bonds)
                {
                    rank += b.Value;
                }
                arr[i, i] = rank;
            }            

            return arr;
        }

        private static double[,] BuildConnectivityMatrix(Molecule m)
        {
            double[,] arr = new double[m.NumOfAtoms, m.NumOfAtoms];

            for(int i = 0; i < m.NumOfAtoms; i++)
            {
                foreach (var b in m.Atoms[i].Bonds)
                {
                    arr[i, b.Key] = b.Value;
                }
            }

            return arr;
        }

        private static double[,] BuildIdentityMatrix(Molecule m)
        {
            double[,] arr = new double[m.NumOfAtoms, m.NumOfAtoms];

            for(int i = 0; i < m.NumOfAtoms; i++)
            {
                arr[i, i] = 1;
            }

            return arr;
        }
    }
}
