using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Diagnostics;
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
        private static bool paramBonds;
        private static double kappa;        
        
        static void Main(string[] args)
        {
            // ..\..\TestFiles\set01.sdf ..\..\TestFiles\ElemBond.txt eem ..\..\outEEM.txt
            // ..\..\TestFiles\set01.sdf ..\..\TestFiles\Element.txt eem ..\..\outEEM.txt
            // ..\..\TestFiles\set01.sdf ..\..\TestFiles\ElemBond.txt mg ..\..\outMG.txt
            // ..\..\TestFiles\set01.sdf ..\..\TestFiles\Element.txt mg ..\..\outMG.txt
            // ..\..\TestFiles\acetonitrile.sdf ..\..\TestFiles\ElemBond.txt mg ..\..\outMG.txt
            // ..\..\TestFiles\acetonitrile.sdf ..\..\TestFiles\Element.txt mg ..\..\outMG.txt 
            if (!CanParseArguments(args))
                return;            
            string sdfFilePath = Path.GetFullPath(Path.Combine(Environment.CurrentDirectory, args[0]));
            LoadMolecules(sdfFilePath);
            string paramsFilePath = Path.GetFullPath(Path.Combine(Environment.CurrentDirectory, args[1]));
            LoadParameters(paramsFilePath);
            string outputFilePath = Path.GetFullPath(Path.Combine(Environment.CurrentDirectory, args[3]));
            switch (args[2])
            {
                case "eem": SolveEEM(outputFilePath); break;
                case "mg": SolveMG(outputFilePath); break;
                case "og": break;
                default: break;
            }

            #region Tests
            /*
            var m = molecules.ElementAt(0);

            Console.WriteLine($"Atom 6: {m.Atoms[5].X}, {m.Atoms[5].Y}, {m.Atoms[5].Z}");
            Console.WriteLine($"Atom 7: {m.Atoms[6].X}, {m.Atoms[6].Y}, {m.Atoms[6].Z}");
            Console.WriteLine($"Distance: {CalculateDistance(m.Atoms[5], m.Atoms[6])}");
            Console.WriteLine();

            #region EEM Matrix solution test
            Console.WriteLine("EEM Matrix:\n");
            double[,] matrix = BuildEEMMatrix(m);
            for (int i = 0; i <= m.NumOfAtoms; i++)
            {
                for (int j = 0; j <= m.NumOfAtoms; j++)
                {
                    Console.Write("{0,6:F2}", matrix[i, j]);
                }
                Console.WriteLine();
            }
            Console.WriteLine();
            var m1 = Matrix<double>.Build.DenseOfArray(matrix);
            var v1 = Vector<double>.Build.Dense(BuildEEMVector(m));
            Console.WriteLine("EEM Vector:\n");
            int count = 0;
            foreach (double a in v1)
            {
                Console.Write("{0,-4}", count + 1);
                Console.WriteLine("{0,7:F4}", a);
                count++;
            }
            Console.WriteLine();
            var res = m1.Solve(v1);
            count = 0;
            Console.WriteLine("Results:\n");
            foreach (double a in res)
            {
                if (count != m.NumOfAtoms)
                {
                    Console.Write("{0,-4}", count + 1);
                    Console.WriteLine("{0,13:F10}", a);
                }
                count++;
            }
            #endregion
            
            #region MG Matrix solution test   
            
            var D = Matrix<double>.Build.DenseOfArray(BuildDegreeMatrix(m));
            var A = Matrix<double>.Build.DenseOfArray(BuildConnectivityMatrix(m));
            var I = Matrix<double>.Build.DenseOfArray(BuildIdentityMatrix(m));
            var V = Vector<double>.Build.Dense(BuildMGVector(m));
            double avgEN = CountGeometricAverageEN(V);
            Matrix<double> S = D - A + I;
            var X = S.Solve(V);
            var deltaX = X - V;
            var results = deltaX * (1.0 / CountGeometricAverageEN(V));

            Console.WriteLine("D matrix:\n");
            for (int i = 0; i < m.NumOfAtoms; i++)
            {
                for (int j = 0; j < m.NumOfAtoms; j++)
                {
                    Console.Write("{0,3}", D[i, j]);
                }
                Console.WriteLine();
            }

            Console.WriteLine();
            Console.WriteLine("A matrix:\n");
            for (int i = 0; i < m.NumOfAtoms; i++)
            {
                for (int j = 0; j < m.NumOfAtoms; j++)
                {
                    Console.Write("{0,3}", A[i, j]);
                }
                Console.WriteLine();
            }

            Console.WriteLine();
            Console.WriteLine("I matrix:\n");
            for (int i = 0; i < m.NumOfAtoms; i++)
            {
                for (int j = 0; j < m.NumOfAtoms; j++)
                {
                    Console.Write("{0,3}", I[i, j]);
                }
                Console.WriteLine();
            }

            Console.WriteLine();
            Console.WriteLine("S matrix:\n");
            for (int i = 0; i < m.NumOfAtoms; i++)
            {
                for (int j = 0; j < m.NumOfAtoms; j++)
                {
                    Console.Write("{0,3}", S[i, j]);
                }
                Console.WriteLine();
            }

            Console.WriteLine("Vector:\n");
            foreach(double val in V)
            {
                Console.Write("{0,6}", val);
            }
            Console.WriteLine();

            Console.WriteLine("X vector:\n");
            foreach (double val in X)
            {
                Console.Write("{0,6:F2}", val);
            }

            Console.WriteLine();
            Console.WriteLine("deltaX matrix:\n");
            foreach (double val in deltaX)
            {
                Console.Write("{0,6:F2}", val);
            }
            Console.WriteLine();            

            Console.WriteLine();
            Console.WriteLine($"Geometric average of electronegativities: {avgEN}\n");
       
            Console.WriteLine("Results:\n");
            int count = 0;
            foreach (double charge in results)
            {
                if (count != m.NumOfAtoms)
                {
                    Console.Write("{0,-4}", count + 1);
                    Console.WriteLine("{0,9:F6}", charge);
                }
                count++;
            }            
            #endregion
            */
            #endregion
        }

        private static bool CanParseArguments(string[] items)
        {
            if(items.Length != 4)
            {
                Console.WriteLine("All four arguments must be entered.");
                return false;
            }
            else if (!(items[2].Equals("eem") || items[2].Equals("mg") || items[2].Equals("og")))
            {
                Console.WriteLine("Incorrect method name, only \"eem\", \"mg\" or \"og\" supported.");
                return false;
            }
            return true;
        }

        private static void LoadParameters(string fileName)
        {
            try
            {
                elemParams = new List<ElementParameters>();
                using (StreamReader reader = File.OpenText(fileName))
                {                    
                    string line;
                    bool firstLine = true;

                    while ((line = reader.ReadLine()) != null)
                    {
                        if (firstLine)
                        {
                            firstLine = false;
                            string[] items = line.Split(' ');
                            kappa = double.Parse(items[1], CultureInfo.InvariantCulture);
                            if (items[0].Equals("Element"))
                            {
                                paramBonds = false;
                            }
                            else
                            {
                                paramBonds = true;
                            }
                            continue;
                        }

                        string[] param = line.Split(' ');
                        ElementParameters parameters;
                        if (paramBonds)
                        {
                            parameters = new ElementParameters(
                                param[0],
                                double.Parse(param[2], CultureInfo.InvariantCulture),
                                double.Parse(param[3], CultureInfo.InvariantCulture),
                                int.Parse(param[1]));
                        }
                        else
                        {
                            parameters = new ElementParameters(
                                param[0],
                                double.Parse(param[1], CultureInfo.InvariantCulture),
                                double.Parse(param[2], CultureInfo.InvariantCulture)
                                );
                        }

                        elemParams.Add(parameters);
                    }
                }
            } catch (Exception ex)
            {
                Console.WriteLine("Could not load parameters. Exception: " + ex.Message);
            }
        }        

        private static void LoadMolecules(string fileName)
        {
            molecules = new List<Molecule>();
            try
            {
                using (StreamReader reader = File.OpenText(fileName))
                {
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
                            lineNum++;
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
                        else if (lineNum <= molecule.NumOfAtoms + molecule.NumOfBonds + 4)
                        {
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
                        }
                        else
                            lineNum++;
                    }
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine("Could not load molecules. Exception: " + ex.Message);
            }
        }

        private static double CalculateDistance(Atom atom1, Atom atom2)
        {
            return Math.Sqrt(
                  (atom2.X - atom1.X) * (atom2.X - atom1.X)
                + (atom2.Y - atom1.Y) * (atom2.Y - atom1.Y) 
                + (atom2.Z - atom1.Z) * (atom2.Z - atom1.Z)
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
        
        //if loaded parameters file does not have specified bond types for parameters, returns 0
        private static int GetHighestBondType(Molecule m, Atom a)
        {
            int bond = 0;

            if (!paramBonds) return bond;

            foreach(var b in a.Bonds)
                if (b.Value > bond)
                    bond = b.Value;

            foreach(var atom in m.Atoms)
                foreach(KeyValuePair<int, int> b in atom.Bonds)
                    if (b.Key == a.ID && b.Value > bond)
                        bond = b.Value;

            return bond;
        }

        private static double[,] BuildEEMMatrix(Molecule m)
        {
            int columns = m.NumOfAtoms + 1;
            int rows = m.NumOfAtoms + 1;
            double[,] arr = new double[rows, columns];

            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    if (j == columns - 1)
                        arr[i, j] = -1;
                    else if (i == rows - 1)
                        arr[i, j] = 1;
                    else if (i == j)
                    {
                        var a = m.Atoms.First(x => x.ID == i + 1);
                        var pars = elemParams.First(x =>
                            x.ElementName.Equals(a.Symbol) &&
                            x.BondType == GetHighestBondType(m, a));
                        arr[i, j] = pars.B;
                    }                    
                    else
                    {
                        var a1 = m.Atoms.First(x => x.ID == i + 1);
                        var a2 = m.Atoms.First(x => x.ID == j + 1);
                        arr[i, j] = (kappa / (CalculateDistance(a1, a2)));
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
                var pars = elemParams.First(x =>
                    x.ElementName.Equals(a.Symbol) &&
                    x.BondType == GetHighestBondType(m, a));
                vector[count] = -pars.A;
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
                    var matrix = Matrix<double>.Build.DenseOfArray(BuildEEMMatrix(molecule));
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
                            file.Write("{0,-4}", count + 1);
                            file.Write("{0,-3}", molecule.Atoms[count].Symbol);
                            file.WriteLine("{0,9:F6}", charge);
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

            foreach(var a in m.Atoms)
            {
                foreach(var b in a.Bonds)
                {
                    arr[b.Key - 1, b.Key - 1] += b.Value;
                    arr[a.ID - 1, a.ID - 1] += b.Value;
                }
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
                    arr[i, b.Key - 1] = b.Value;
                    arr[b.Key - 1, i] = b.Value;                    
                }
            }

            return arr;
        }

        private static double[,] BuildIdentityMatrix(Molecule m)
        {
            double[,] arr = new double[m.NumOfAtoms, m.NumOfAtoms];

            for(int i = 0; i < m.NumOfAtoms; i++)
                arr[i, i] = 1;

            return arr;
        }

        private static double[] BuildMGVector(Molecule m)
        {
            double[] arr = new double[m.NumOfAtoms];

            for(int i = 0; i < m.NumOfAtoms; i++)
                arr[i] = GetElectronegativity(m.Atoms[i]);

            return arr;
        }

        private static double CountGeometricAverageEN(Vector<double> vector)
        {
            double avgEN = 1;

            for(int i = 0; i < vector.Count; i++)
                avgEN *= vector[i];
            avgEN = Math.Pow(avgEN, 1.0 / vector.Count);

            return avgEN;
        }

        private static void SolveMG(string outputFilePath)
        {
            using (StreamWriter file = new StreamWriter(outputFilePath))
            {
                foreach (Molecule molecule in molecules)
                {
                    var a = Matrix<double>.Build.DenseOfArray(BuildConnectivityMatrix(molecule));
                    var d = Matrix<double>.Build.DenseOfArray(BuildDegreeMatrix(molecule));
                    var i = Matrix<double>.Build.DenseOfArray(BuildIdentityMatrix(molecule));
                    var vector = Vector<double>.Build.Dense(BuildMGVector(molecule));

                    Matrix<double> s = d - a + i;
                    var x = s.Solve(vector);                    
                    var charges = (x - vector) * (1.0 / CountGeometricAverageEN(vector));

                    file.WriteLine($"NSC_{molecule.NSC}");
                    file.WriteLine(molecule.NumOfAtoms);
                    int count = 0;
                    foreach (double charge in charges)
                    {
                        if (count != molecule.NumOfAtoms)
                        {
                            molecule.Atoms[count].Charge = charge;
                            file.Write("{0,-4}", count + 1);
                            file.Write("{0,-3}", molecule.Atoms[count].Symbol);
                            file.WriteLine("{0,9:F6}", charge);
                            count++;
                        }
                    }
                    file.WriteLine("$$$$");
                }
            }
        }
    }
}
