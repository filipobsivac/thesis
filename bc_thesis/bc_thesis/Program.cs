using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Reflection;
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
            // ..\..\TestFiles\set01.sdf ..\..\TestFiles\ElemBond.txt mgc ..\..\outMGC.txt
            // ..\..\TestFiles\set01.sdf ..\..\TestFiles\Element.txt mgc ..\..\outMGC.txt
            // ..\..\TestFiles\acetonitrile.sdf ..\..\TestFiles\ElemBond.txt mgc ..\..\outMGC.txt
            // ..\..\TestFiles\acetonitrile.sdf ..\..\TestFiles\Element.txt mgc ..\..\outMGC.txt 
            // ..\..\outEEM.txt ..\..\outMGC.txt stats ..\..\outSTATS.txt
            if (!CanParseArguments(args))
                return;                        
            string firstFilePath = Path.GetFullPath(Path.Combine(Environment.CurrentDirectory, args[0]));
            string secondFilePath = Path.GetFullPath(Path.Combine(Environment.CurrentDirectory, args[1]));
            string outputFilePath = Path.GetFullPath(Path.Combine(Environment.CurrentDirectory, args[3]));
            if (!args[2].Equals("stats"))
            {
                LoadMolecules(firstFilePath);
                LoadParameters(secondFilePath);
            }
            switch (args[2])
            {
                case "eem": SolveEEM(outputFilePath); break;
                case "mgc": SolveMGC(outputFilePath); break;
                case "stats": GenerateStatistics(firstFilePath, secondFilePath, outputFilePath); break;
                case "ogc": break;
                default: break;
            }

            #region OGC Tests
            foreach (var atom in molecules[0].Atoms)
            {
                Console.Write($"{atom.Symbol} {atom.ID} - ");
                foreach(var b in atom.Bonds)
                    Console.Write($"{b.Key}|{b.Value} ");
                Console.WriteLine();
            }

            Console.WriteLine();

            var m = BuildOGCMolecule(molecules[0]);
            foreach (var atom in m.Atoms)
            {
                Console.Write($"{atom.Symbol} {atom.OrbitalID} - ");
                foreach(var b in atom.OrbitalBonds)
                    Console.Write($"{b.Key}|{b.Value} ");
                Console.WriteLine();
            }

            Console.WriteLine();

            var d = BuildDegreeMatrixOGC(m);
            for (int i = 0; i < m.NumOfAtoms; i++)
            {
                for (int j = 0; j < m.NumOfAtoms; j++)
                    Console.Write($"{d[i,j]} ");
                Console.WriteLine();
            }

            Console.WriteLine();

            var a = BuildConnectivityMatrixOGC(m);
            for (int i = 0; i < m.NumOfAtoms; i++)
            {
                for (int j = 0; j < m.NumOfAtoms; j++)
                    Console.Write($"{a[i, j]} ");
                Console.WriteLine();
            }

            Console.WriteLine();

            var id = BuildIdentityMatrix(m);
            for (int i = 0; i < m.NumOfAtoms; i++)
            {
                for (int j = 0; j < m.NumOfAtoms; j++)
                    Console.Write($"{id[i, j]} ");
                Console.WriteLine();
            }

            Console.WriteLine();

            var vector = BuildOGCVector(m);
            Matrix<double> s = d - a + id;
            var equalizedEN = s.Solve(vector);
            for (int i = 0; i < m.NumOfAtoms; i++)
            {
                var orbType = m.Atoms[i].OrbitalBonds.First(x => !x.Value.Equals("x")).Value;
                var orbString = "n";
                if (!orbType.Equals("n"))
                {
                    var bondedOrbID = m.Atoms[i].OrbitalBonds.FirstOrDefault(x => !x.Value.Equals("x")).Key;
                    var bondedOrbSymbol = m.Atoms.Find(x => x.OrbitalID == bondedOrbID).Symbol;
                    orbString = $"{orbType}({m.Atoms[i].Symbol}-{bondedOrbSymbol})";
                }                
                Console.WriteLine($"{orbString}\t{vector[i]}\t{equalizedEN[i]}");
            }

            Console.ReadKey();
            #endregion
        }

        private static bool CanParseArguments(string[] items)
        {
            if(items.Length != 4)
            {
                Console.WriteLine("Arguments must be entered in the following format:");
                Console.WriteLine("\"<sdf file> <parameters file> <method> <output file>\" or \"<first molecules file> <second molecules file> \"stats\" <output file>\"");
                return false;
            }
            else if (!(items[2].Equals("eem") || items[2].Equals("mgc") || items[2].Equals("ogc") || items[2].Equals("stats")))
            {
                Console.WriteLine("Incorrect method name. Only \"eem\", \"mgc\", \"ogc\" or \"stats\" supported.");
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
                            Atom atom = new Atom();                            
                            atom.ID = lineNum - 4;
                            atom.Symbol = line.Substring(31, 3).Trim(' ');
                            atom.X = double.Parse(line.Substring(0, 10), CultureInfo.InvariantCulture);
                            atom.Y = double.Parse(line.Substring(10, 10), CultureInfo.InvariantCulture);
                            atom.Z = double.Parse(line.Substring(20, 10), CultureInfo.InvariantCulture);
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
            try
            {
                using (StreamReader reader = File.OpenText(@"..\..\Tables\ElementEN.csv"))
                {
                    string line;                    
                    while ((line = reader.ReadLine()) != null)
                    {
                        string[] items = line.Split(';');
                        if (items[0].Equals(atom.Symbol))
                            return double.Parse(items[1], CultureInfo.InvariantCulture);
                    }
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine("Could not read Element EN .csv file. Exception: " + ex.Message);
            }
            return 0;
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
                foreach(var b in atom.Bonds)
                    if (b.Key == a.ID && b.Value > bond)
                        bond = b.Value;

            return bond;
        }

        //needs non-oriented graph
        private static int GetHighestBondTypeOGC(Atom atom)
        {
            int bondType = 1;

            switch (atom.Symbol)
            {
                case "H":
                    return 1;
                case "C":
                    bondType = 1;
                    foreach (var bond in atom.Bonds)
                    {
                        if (bond.Value == 2)
                            bondType = 2;
                        if (bond.Value == 3)
                            bondType = 3;
                    }
                    return bondType;
                case "N":
                    bondType = 1;
                    foreach (var bond in atom.Bonds)
                    {
                        if (bond.Value == 2)
                            bondType = 2;
                        if (bond.Value == 3)
                            bondType = 3;
                    }
                    return bondType;
                case "O":
                    bondType = 1;
                    foreach (var bond in atom.Bonds)
                    {
                        if (bond.Value == 2)
                            bondType = 2;
                    }
                    return bondType;
                case "S":
                    bondType = 1;
                    foreach (var bond in atom.Bonds)
                    {
                        if (bond.Value == 2)
                            bondType = 2;
                    }
                    return bondType;
                default:
                    return 0;
            }
        }

        private static Matrix<double> BuildEEMMatrix(Molecule m)
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
            return Matrix<double>.Build.DenseOfArray(arr);
        }

        private static Vector<double> BuildEEMVector(Molecule m)
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
            return Vector<double>.Build.Dense(vector);
        }

        private static void SolveEEM(string outputFilePath)
        {
            try
            {
                using (StreamWriter file = new StreamWriter(outputFilePath))
                {
                    foreach (Molecule molecule in molecules)
                    {
                        var matrix = BuildEEMMatrix(molecule);
                        var vector = BuildEEMVector(molecule);
                        var results = matrix.Solve(vector);
                        SaveResults(file, molecule, results);
                    }
                }
            }catch(Exception ex)
            {
                Console.WriteLine("Could not save results. Exception: " + ex.Message);
            }
        }        

        private static Matrix<double> BuildDegreeMatrix(Molecule m)
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

            return Matrix<double>.Build.DenseOfArray(arr);
        }

        private static Matrix<double> BuildDegreeMatrixOGC(Molecule m)
        {
            double[,] arr = new double[m.NumOfAtoms, m.NumOfAtoms];

            int x = 0;
            foreach(var a in m.Atoms)
            {
                int n = 0;
                foreach(var b in a.OrbitalBonds)
                    if (!b.Value.Equals("n"))
                        n++;
                arr[x, x] = n;
                x++;
            }

            return Matrix<double>.Build.DenseOfArray(arr);
        }

        private static Matrix<double> BuildConnectivityMatrix(Molecule m)
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

            return Matrix<double>.Build.DenseOfArray(arr);
        }

        private static Matrix<double> BuildConnectivityMatrixOGC(Molecule m)
        {
            double[,] arr = new double[m.NumOfAtoms, m.NumOfAtoms];

            foreach(var a in m.Atoms)
            {
                foreach(var b in a.OrbitalBonds)
                {
                    if (!b.Value.Equals("n"))
                    {
                        var o = m.Atoms.Find(x => x.OrbitalID.Equals(b.Key));
                        arr[a.OGCID - 1, o.OGCID - 1] = 1;
                    }
                }
            }

            return Matrix<double>.Build.DenseOfArray(arr);
        }

        private static Matrix<double> BuildIdentityMatrix(Molecule m)
        {
            double[,] arr = new double[m.NumOfAtoms, m.NumOfAtoms];

            for(int i = 0; i < m.NumOfAtoms; i++)
                arr[i, i] = 1;

            return Matrix<double>.Build.DenseOfArray(arr);
        }

        private static Vector<double> BuildMGCVector(Molecule m)
        {
            double[] arr = new double[m.NumOfAtoms];

            for(int i = 0; i < m.NumOfAtoms; i++)
                arr[i] = GetElectronegativity(m.Atoms[i]);

            return Vector<double>.Build.Dense(arr);
        }

        private static Vector<double> BuildOGCVector(Molecule m)
        {
            double[] arr = new double[m.NumOfAtoms];

            foreach (var a in m.Atoms)
            {
                var bond = a.OrbitalBonds.First(x => !x.Value.Equals("x")).Value;
                var en = a.OrbitalCharges[bond];
                arr[a.OGCID - 1] = en;
            }

            return Vector<double>.Build.Dense(arr);
        }

        private static void SolveMGC(string outputFilePath)
        {
            try
            {
                using (StreamWriter file = new StreamWriter(outputFilePath))
                {
                    foreach (Molecule molecule in molecules)
                    {
                        var a = BuildConnectivityMatrix(molecule);
                        var d = BuildDegreeMatrix(molecule);
                        var i = BuildIdentityMatrix(molecule);
                        var vector = BuildMGCVector(molecule);
                        double avgEN = 1;
                        for (int j = 0; j < vector.Count; j++)
                            avgEN *= vector[j];
                        avgEN = Math.Pow(avgEN, 1.0 / vector.Count);
                        Matrix<double> s = d - a + i;
                        var x = s.Solve(vector);
                        var results = (x - vector) * (1.0 / avgEN);

                        SaveResults(file, molecule, results);
                    }
                }
            }catch(Exception ex)
            {
                Console.WriteLine("Could not save results. Exception: " + ex.Message);
            }
        }

        private static void SolveOGC(string outputFilePath)
        {
            try
            {
                using (StreamWriter file = new StreamWriter(outputFilePath))
                {
                    foreach (Molecule molecule in molecules)
                    {
                        var m = BuildOGCMolecule(molecule);
                        var a = BuildConnectivityMatrixOGC(m);
                        var d = BuildDegreeMatrixOGC(m);
                        var i = BuildIdentityMatrix(m);
                        var vector = BuildOGCVector(m);
                        Matrix<double> s = d - a + i;
                        var x = s.Solve(vector);
                        //TODO
                    }
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine("Could not save results. Exception: " + ex.Message);
            }
        }

        private static void SaveResults(StreamWriter file, Molecule molecule, Vector<double> results)
        {
            file.WriteLine($"NSC_{molecule.NSC}");
            file.WriteLine(molecule.NumOfAtoms);
            int count = 0;
            foreach (double charge in results)
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

        private static List<Molecule> LoadMoleculesFromOutputFile(string filePath)
        {
            List<Molecule> set = new List<Molecule>();
            try
            {
                using (StreamReader reader = File.OpenText(filePath))
                {
                    string line;
                    int lineNum = 1;
                    Molecule molecule = new Molecule();
                    while ((line = reader.ReadLine()) != null)
                    {
                        if (lineNum == 1)
                        {
                            molecule = new Molecule();
                            int nsc = int.Parse(line.Substring(4));
                            molecule.NSC = nsc;
                            lineNum++;
                        }
                        else if (lineNum == 2)
                        {
                            molecule.NumOfAtoms = int.Parse(line);
                            lineNum++;
                        }
                        else if (lineNum <= molecule.NumOfAtoms + 2)
                        {
                            Atom atom = new Atom();
                            atom.ID = int.Parse(line.Substring(0, 4));
                            atom.Symbol = line.Substring(4, 3);
                            atom.Charge = double.Parse(line.Substring(7));
                            molecule.Atoms.Add(atom);
                            lineNum++;
                        }
                        else
                        {
                            set.Add(molecule);
                            lineNum = 1;
                        }
                    }
                }
            }catch(Exception ex)
            {
                Console.WriteLine("Could not load molecules. Exception: " + ex.Message);
            }
            return set;
        }

        private static void GenerateStatistics(string firstFilePath, string secondFilePath, string outputFilePath)
        {
            try
            {
                List<Molecule> firstSet = LoadMoleculesFromOutputFile(firstFilePath);
                List<Molecule> secondSet = LoadMoleculesFromOutputFile(secondFilePath);

                using (StreamWriter file = new StreamWriter(outputFilePath))
                {
                    double avg_d_avg = 0;
                    double avg_d_max = 0;
                    double avg_rmsd = 0;
                    double avg_pearson = 0;

                    for(int i = 0; i < firstSet.Count; i++)
                    {
                        double d_max = 0;

                        double d_avg = 0;
                        double sum_avg = 0;

                        double rmsd = 0;
                        double sum_rmsd = 0;

                        double pearson = 0;
                        double sum_xy = 0;
                        double sum_xsq = 0;
                        double sum_ysq = 0;
                        double sum_x = 0;
                        double sum_y = 0;

                        for (int j = 0; j < firstSet.ElementAt(i).NumOfAtoms; j++)
                        {
                            double x = firstSet.ElementAt(i).Atoms.ElementAt(j).Charge;
                            double y = secondSet.ElementAt(i).Atoms.ElementAt(j).Charge;

                            if(Math.Abs(x - y) > d_max)
                                d_max = Math.Abs(x - y);

                            sum_avg += Math.Abs(x - y);
                            sum_rmsd += (x - y) * (x - y);
                            sum_x += x;
                            sum_y += y;
                            sum_xy += x * y;
                            sum_xsq += x * x;
                            sum_ysq += y * y;
                        }

                        d_avg = sum_avg / (double)firstSet.ElementAt(i).NumOfAtoms;
                        rmsd = Math.Sqrt(sum_rmsd / (double)firstSet.ElementAt(i).NumOfAtoms);
                        int n = firstSet.ElementAt(i).NumOfAtoms;
                        pearson = (n * sum_xy - sum_x * sum_y) / ( Math.Sqrt(n * sum_xsq - sum_x * sum_x) * Math.Sqrt(n * sum_ysq - sum_y * sum_y) );

                        avg_d_avg += d_avg;
                        avg_d_max += d_max;//perhaps just save highest value instead of average?
                        avg_rmsd += rmsd;
                        avg_pearson += pearson;

                        file.WriteLine($"NSC_{firstSet.ElementAt(i).NSC}");
                        file.WriteLine($"Largest absolute difference: {d_max}");
                        file.WriteLine($"Average absolute difference: {d_avg}");
                        file.WriteLine($"Root-mean-square deviation: {rmsd}");
                        file.WriteLine($"Pearson correlation coefficient: {pearson}");
                        file.WriteLine("$$$$");
                    }

                    string avgStatsFilePath = Path.ChangeExtension(outputFilePath, "AvgStats.txt");                    
                    using (StreamWriter avgStatsFile = new StreamWriter(avgStatsFilePath))
                    {
                        avgStatsFile.WriteLine("Average values for the entire set of molecules:");
                        avgStatsFile.WriteLine($"Largest absolute difference: {avg_d_max / (double)firstSet.Count}");
                        avgStatsFile.WriteLine($"Average absolute difference: {avg_d_avg / (double)firstSet.Count}");
                        avgStatsFile.WriteLine($"Root-mean-square deviation: {avg_rmsd / (double)firstSet.Count}");
                        avgStatsFile.WriteLine($"Pearson correlation coefficient: {avg_pearson / (double)firstSet.Count}");
                    }

                    GNUPlot(firstSet, secondSet, outputFilePath);
                }
            }
            catch(Exception ex)
            {
                Console.WriteLine("Could not save results. Exception: " + ex.Message);
            }
        }

        private static void GNUPlot(List<Molecule> firstSet, List<Molecule> secondSet, string outputFilePath)
        {
            Process plotProcess = new Process();
            plotProcess.StartInfo.FileName = Path.GetFullPath(Path.Combine(Environment.CurrentDirectory, @"..\..\gnuplot\bin\gnuplot.exe"));
            plotProcess.StartInfo.RedirectStandardInput = true;
            plotProcess.StartInfo.UseShellExecute = false;
            plotProcess.Start();

            NumberFormatInfo nfi = new NumberFormatInfo();
            nfi.NumberDecimalSeparator = ".";

            string inputFilePath = Path.ChangeExtension(outputFilePath, "InputFile.txt");
            string graphFilePath = Path.ChangeExtension(outputFilePath, @"Graph.png");

            using (StreamWriter inputFile = new StreamWriter(inputFilePath))
                for (int i = 0; i < firstSet.Count; i++)
                    for (int j = 0; j < firstSet.ElementAt(i).NumOfAtoms; j++)
                    {
                        string a = firstSet.ElementAt(i).Atoms.ElementAt(j).Charge.ToString(nfi);
                        string b = secondSet.ElementAt(i).Atoms.ElementAt(j).Charge.ToString(nfi);
                        int c = GetAtomColor(firstSet.ElementAt(i).Atoms.ElementAt(j).Symbol);
                        inputFile.WriteLine($"{a} {b} {c}");
                    }                        

            using (StreamWriter gnuplotFile = plotProcess.StandardInput)
            {
                gnuplotFile.WriteLine("set terminal png");
                gnuplotFile.WriteLine("set title \"Partial atomic charge correlation graph\"");
                gnuplotFile.WriteLine("set xlabel \"First set of atom charges\"");
                gnuplotFile.WriteLine("set ylabel \"Second set of atom charges\"");
                gnuplotFile.WriteLine($"set output '{graphFilePath}'");
                gnuplotFile.WriteLine("set label \"H\" at graph 0.05,0.95 tc \"black\"");
                gnuplotFile.WriteLine("set label \"C\" at graph 0.05,0.90 tc \"purple\"");
                gnuplotFile.WriteLine("set label \"N\" at graph 0.05,0.85 tc \"red\"");
                gnuplotFile.WriteLine("set label \"O\" at graph 0.05,0.80 tc \"orange\"");
                gnuplotFile.WriteLine("set label \"S\" at graph 0.05,0.75 tc \"yellow\"");
                gnuplotFile.WriteLine("unset colorbox");
                gnuplotFile.WriteLine($"plot '{inputFilePath}' u 1:2:3 with points palette notitle");
            }
        }

        private static int GetAtomColor(string symbol)
        {
            switch (symbol.Trim())
            {
                case "H": return 1;
                case "C": return 2;
                case "N": return 3;
                case "O": return 4;
                case "S": return 5;
                default: return 0;
            }
        }

        private static void SetOrbitalCharges(Molecule molecule)
        {
            int bondType = 1;
            foreach (var atom in molecule.Atoms)
            {
                switch (atom.Symbol)
                {
                    case "H":
                        atom.OrbitalCharges.Add("s", GetOrbitalCharge("H", "s_"));
                        break;
                    case "C":
                        bondType = GetHighestBondTypeOGC(atom);
                        switch (bondType)
                        {
                            case 1:                                
                                atom.OrbitalCharges.Add("s", GetOrbitalCharge("C", "te_ te te te"));
                                break;
                            case 2:
                                atom.OrbitalCharges.Add("s", GetOrbitalCharge("C", "tr_ tr tr pp"));
                                atom.OrbitalCharges.Add("p", GetOrbitalCharge("C", "tr tr tr pp_"));                                
                                break;
                            case 3:
                                atom.OrbitalCharges.Add("s", GetOrbitalCharge("C", "di_ di pp pp"));
                                atom.OrbitalCharges.Add("p", GetOrbitalCharge("C", "di di pp_ pp"));
                                break;
                            default: break;
                        }
                        break;
                    case "N":
                        bondType = GetHighestBondTypeOGC(atom);
                        switch (bondType)
                        {
                            case 1:                                
                                atom.OrbitalCharges.Add("s", GetOrbitalCharge("N", "te2 te_ te te"));
                                atom.OrbitalCharges.Add("n", GetOrbitalCharge("N", "te2_ te te te"));
                                break;
                            case 2:
                                atom.OrbitalCharges.Add("s", GetOrbitalCharge("N", "tr2 tr_ tr pp"));
                                atom.OrbitalCharges.Add("p", GetOrbitalCharge("N", "tr2 tr tr pp_"));                                
                                atom.OrbitalCharges.Add("n", GetOrbitalCharge("N", "tr2_ tr tr pp"));
                                break;
                            case 3:
                                atom.OrbitalCharges.Add("s", GetOrbitalCharge("N", "di2 di_ pp pp"));
                                atom.OrbitalCharges.Add("p", GetOrbitalCharge("N", "di2 di pp_ pp"));
                                atom.OrbitalCharges.Add("n", GetOrbitalCharge("N", "di2_ di pp pp"));
                                break;
                            default: break;
                        }
                        break;
                    case "O":
                        bondType = GetHighestBondTypeOGC(atom);
                        switch (bondType)
                        {
                            case 1:
                                atom.OrbitalCharges.Add("s", GetOrbitalCharge("O", "te2 te2 te_ te"));
                                atom.OrbitalCharges.Add("n", GetOrbitalCharge("O", "te2_ te2 te te"));                                
                                break;
                            case 2:
                                atom.OrbitalCharges.Add("s", GetOrbitalCharge("O", "tr2 tr2 tr_ pp"));
                                atom.OrbitalCharges.Add("p", GetOrbitalCharge("O", "tr2 tr2 tr pp_"));
                                atom.OrbitalCharges.Add("n", GetOrbitalCharge("O", "tr2_ tr2 tr pp"));
                                break;
                            default: break;
                        }
                        break;
                    case "S":
                        bondType = GetHighestBondTypeOGC(atom);
                        switch (bondType)
                        {
                            case 1:
                                atom.OrbitalCharges.Add("s", GetOrbitalCharge("S", "te2 te2 te_ te"));
                                atom.OrbitalCharges.Add("n", GetOrbitalCharge("S", "te2_ te2 te te"));                                
                                break;
                            case 2:
                                atom.OrbitalCharges.Add("s", GetOrbitalCharge("S", "tr2 tr2 tr_ pp"));
                                atom.OrbitalCharges.Add("p", GetOrbitalCharge("S", "tr2 tr2 tr pp_"));
                                atom.OrbitalCharges.Add("n", GetOrbitalCharge("S", "tr2_ tr2 tr pp"));
                                break;
                            default: break;
                        }
                        break;
                    default: break;
                }
            }
        }

        private static double GetOrbitalCharge(string symbol, string state)
        {
            try
            {
                using (StreamReader reader = File.OpenText(@"..\..\Tables\OrbitalHardnessAndEN.csv"))
                {
                    string line;
                    bool skip = true;
                    while ((line = reader.ReadLine()) != null)
                    {
                        if (skip)
                        {
                            skip = false;
                            continue;
                        }
                        string[] items = line.Split(';');
                        if (items[0].Equals(symbol) && items[1].Equals(state))
                            return double.Parse(items[2], CultureInfo.InvariantCulture);
                    }
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine("Could not get orbital charges. Exception: " + ex.Message);
            }
            return 0;
        }

        private static Atom GetFreeOrbital(Molecule m, int id)
        {
            return m.Atoms.Find( a => (a.ID == id) && !a.OrbitalBonds.Any( b => !b.Value.Equals("x") ) );
        }

        private static void MakeBondsNonOriented(Molecule m)
        {
            foreach (var a in m.Atoms)
            {
                foreach(var b in a.Bonds)
                {
                    var bondedAtom = m.Atoms.Find(x => x.ID == b.Key);
                    if(!bondedAtom.Bonds.Any(x => x.Key == a.ID))
                        bondedAtom.Bonds.Add(a.ID, b.Value);
                }
            }
        }
                
        private static Molecule BuildOGCMolecule(Molecule m)
        {
            MakeBondsNonOriented(m);
            SetOrbitalCharges(m);
            Molecule ogcMolecule = new Molecule();            
            ogcMolecule.NSC = m.NSC;
            ogcMolecule.NumOfAtoms = 0;
            ogcMolecule.NumOfBonds = m.NumOfBonds;

            int x = 1;
            foreach (var a in m.Atoms)
            {
                int n;
                if (a.Symbol.Trim().Equals("H"))
                    n = 1;
                else
                    n = 4;
                for (int i = 1; i <= n; i++)
                {
                    Atom atom = new Atom();
                    atom.OGCID = x;
                    atom.ID = a.ID;
                    atom.Symbol = a.Symbol;
                    atom.Bonds = a.Bonds;
                    atom.OrbitalCharges = a.OrbitalCharges;
                    atom.OrbitalID = $"{atom.ID}x{i}";
                    ogcMolecule.Atoms.Add(atom);
                    ogcMolecule.NumOfAtoms++;
                    x++;
                }                
            }

            foreach(var a in m.Atoms)
            {
                var orbitals = ogcMolecule.Atoms.FindAll(orb => orb.ID == a.ID);
                foreach (var o in orbitals)
                {
                    var orbs = orbitals.FindAll(orb => orb.OrbitalID != o.OrbitalID);
                    foreach(var orb in orbs)
                        o.OrbitalBonds.Add(orb.OrbitalID, "x");
                }

                if (a.Symbol.Trim().Equals("N"))
                {
                    var firstOrbital = GetFreeOrbital(ogcMolecule, a.ID);
                    firstOrbital.OrbitalBonds.Add(firstOrbital.ID.ToString(), "n");
                }

                if (a.Symbol.Trim().Equals("O"))
                {
                    for(int i = 0; i < 2; i++)
                    {
                        var firstOrbital = GetFreeOrbital(ogcMolecule, a.ID);
                        firstOrbital.OrbitalBonds.Add(firstOrbital.ID.ToString(), "n");
                    }                    
                }

                foreach (var b in a.Bonds)
                {
                    var firstOrbital = GetFreeOrbital(ogcMolecule, b.Key);
                    if(firstOrbital == null)
                        continue;
                    if(b.Value == 1)
                    {
                        var sOrbital = GetFreeOrbital(ogcMolecule, a.ID);
                        firstOrbital.OrbitalBonds.Add(sOrbital.OrbitalID, "s");
                        sOrbital.OrbitalBonds.Add(firstOrbital.OrbitalID, "s");
                    }
                    if(b.Value == 2)
                    {
                        var pOrbital = GetFreeOrbital(ogcMolecule, a.ID);
                        firstOrbital.OrbitalBonds.Add(pOrbital.OrbitalID, "p");
                        pOrbital.OrbitalBonds.Add(firstOrbital.OrbitalID, "p");
                        firstOrbital = GetFreeOrbital(ogcMolecule, b.Key);
                        var sOrbital = GetFreeOrbital(ogcMolecule, a.ID);
                        firstOrbital.OrbitalBonds.Add(sOrbital.OrbitalID, "s");
                        sOrbital.OrbitalBonds.Add(firstOrbital.OrbitalID, "s");
                    }
                    if(b.Value == 3)
                    {
                        var pOrbital = GetFreeOrbital(ogcMolecule, a.ID);
                        firstOrbital.OrbitalBonds.Add(pOrbital.OrbitalID, "p");
                        pOrbital.OrbitalBonds.Add(firstOrbital.OrbitalID, "p");
                        firstOrbital = GetFreeOrbital(ogcMolecule, b.Key);
                        var sOrbital = GetFreeOrbital(ogcMolecule, a.ID);
                        firstOrbital.OrbitalBonds.Add(sOrbital.OrbitalID, "s");
                        sOrbital.OrbitalBonds.Add(firstOrbital.OrbitalID, "s");
                        firstOrbital = GetFreeOrbital(ogcMolecule, b.Key);
                        var p2Orbital = GetFreeOrbital(ogcMolecule, a.ID);
                        firstOrbital.OrbitalBonds.Add(p2Orbital.OrbitalID, "p");
                        p2Orbital.OrbitalBonds.Add(firstOrbital.OrbitalID, "p");
                    }
                }
            }

            return ogcMolecule;
        }
    }
}
