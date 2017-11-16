using bc_thesis.Classes;
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
            // eem ..\..\TestFiles\set01.sdf ..\..\outEEMb.txt ..\..\TestFiles\ElemBond.txt 
            // eem ..\..\TestFiles\DTP_large.sdf ..\..\outEEMb.txt ..\..\TestFiles\ElemBond.txt 
            // eem ..\..\TestFiles\set01.sdf ..\..\outEEM.txt ..\..\TestFiles\Element.txt
            // eem ..\..\TestFiles\DTP_large.sdf ..\..\outEEM.txt ..\..\TestFiles\Element.txt 
            // mgc ..\..\TestFiles\set01.sdf ..\..\outMGC.txt
            // ogc ..\..\TestFiles\set01.sdf ..\..\outOGC.txt
            // mgc ..\..\TestFiles\acetonitrile.sdf ..\..\outMGC.txt
            // ogc ..\..\TestFiles\acetonitrile.sdf ..\..\outOGC.txt
            // stats ..\..\outEEMb.txt ..\..\outEEM.txt ..\..\outSTATS.txt y
            // stats ..\..\outEEMb.txt ..\..\outMGC.txt ..\..\outSTATS.txt y
            // stats ..\..\outEEMb.txt ..\..\outOGC.txt ..\..\outSTATS.txt y
            // stats ..\..\outMGC.txt ..\..\outOGC.txt ..\..\outSTATS.txt y
            if (!CanParseArguments(args))
                return;
            if (args[0].Equals("stats"))
            {
                string firstFilePath = Path.GetFullPath(Path.Combine(Environment.CurrentDirectory, args[1]));
                string secondFilePath = Path.GetFullPath(Path.Combine(Environment.CurrentDirectory, args[2]));
                string outputFilePath = Path.GetFullPath(Path.Combine(Environment.CurrentDirectory, args[3]));
                if (args[4].Equals("y"))
                    paramBonds = true;
                if (args[4].Equals("n"))
                    paramBonds = false;
                GenerateStatistics(firstFilePath, secondFilePath, outputFilePath);
            }
            else
            {
                string moleculesFilePath = Path.GetFullPath(Path.Combine(Environment.CurrentDirectory, args[1]));
                LoadMoleculesV3000(moleculesFilePath);
                string outputFilePath = Path.GetFullPath(Path.Combine(Environment.CurrentDirectory, args[2]));
                if (args[0].Equals("eem"))
                {
                    string parametersFilePath = Path.GetFullPath(Path.Combine(Environment.CurrentDirectory, args[3]));
                    LoadParameters(parametersFilePath);
                    SolveEEM(outputFilePath);
                }
                if (args[0].Equals("mgc"))
                    SolveMGC(outputFilePath);
                if (args[0].Equals("ogc"))
                    SolveOGC(outputFilePath);                    
            }

            #region OGC Tests
            /*
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
            for (int i = 0; i < m.Atoms.Count; i++)
            {
                for (int j = 0; j < m.Atoms.Count; j++)
                    Console.Write($"{d[i,j]} ");
                Console.WriteLine();
            }

            Console.WriteLine();

            var a = BuildConnectivityMatrixOGC(m);
            for (int i = 0; i < m.Atoms.Count; i++)
            {
                for (int j = 0; j < m.Atoms.Count; j++)
                    Console.Write($"{a[i, j]} ");
                Console.WriteLine();
            }

            Console.WriteLine();

            var id = BuildIdentityMatrix(m);
            for (int i = 0; i < m.Atoms.Count; i++)
            {
                for (int j = 0; j < m.Atoms.Count; j++)
                    Console.Write($"{id[i, j]} ");
                Console.WriteLine();
            }

            Console.WriteLine();

            var vector = BuildOGCVector(m);
            Matrix<double> s = d - a + id;
            var equalizedEN = s.Solve(vector);
            var deltaEN = equalizedEN - vector;   
            double dividend = 0;
            double divisor = 0;
            int n = 0;
            foreach (var orbital in m.Atoms)
            {
                var bond = orbital.OrbitalBonds.First(x => !x.Value.Equals("x")).Value;
                var orbHardness = orbital.OrbitalHardnesses.First(x => x.Key.Equals(bond)).Value;
                dividend += deltaEN[n] / orbHardness;
                divisor += (deltaEN[n] * deltaEN[n]) / (Math.Pow(orbHardness, 3) * GetCovalentRadius(orbital.Symbol.Trim() + orbital.HighestBondType.ToString()));
                n++;
            }
            double dm = dividend / divisor;
            n = 0;
            foreach (var orbital in m.Atoms)
            {
                var bond = orbital.OrbitalBonds.First(x => !x.Value.Equals("x")).Value;
                var orbHardness = orbital.OrbitalHardnesses.First(x => x.Key.Equals(bond)).Value;
                orbital.OrbitalCharge = (deltaEN[n] / orbHardness) - (((deltaEN[n] * deltaEN[n]) * dm) / (Math.Pow(orbHardness, 3) * GetCovalentRadius(orbital.Symbol.Trim() + orbital.HighestBondType.ToString())));                                
                n++;
            }

            for (int i = 0; i < m.Atoms.Count; i++)
            {
                var orbType = m.Atoms[i].OrbitalBonds.First(x => !x.Value.Equals("x")).Value;
                var orbHardness = m.Atoms[i].OrbitalHardnesses.First(x => x.Key.Equals(orbType)).Value;
                var covRad = GetCovalentRadius(m.Atoms[i].Symbol + m.Atoms[i].HighestBondType.ToString());
                var orbString = "n";
                if (!orbType.Equals("n"))
                {
                    var bondedOrbID = m.Atoms[i].OrbitalBonds.FirstOrDefault(x => !x.Value.Equals("x")).Key;
                    var bondedOrbSymbol = m.Atoms.Find(x => x.OrbitalID == bondedOrbID).Symbol;
                    orbString = $"{orbType}({m.Atoms[i].Symbol}-{bondedOrbSymbol})";
                }
                Console.WriteLine($"{orbString}\t{vector[i]}\t{Math.Round(equalizedEN[i], 3)}\t{orbHardness}\t{covRad}\t{Math.Round(m.Atoms[i].OrbitalCharge, 3)}");
            }
            Console.WriteLine();
            Console.WriteLine($"Dm = {dm}");
            Console.WriteLine();

            double[] arr = new double[m.NumOfAtoms];
            foreach (var atom in m.Atoms)
            {
                var orbitals = m.Atoms.FindAll(o => o.ID == atom.ID);
                arr[atom.ID - 1] = orbitals.Sum(x => x.OrbitalCharge);                
            }

            for(int i = 0; i < molecules[0].Atoms.Count; i++)
                Console.WriteLine($"{molecules[0].Atoms[i].Symbol} = {Math.Round(arr[i], 3)}");
            */
            #endregion

            Console.ReadKey();
        }

        private static void PrintHelp()
        {
            Console.WriteLine("Arguments must be entered in the following format:");
            Console.WriteLine("\"<method> <molecules file> <output file> <parameters file - optional>\" or \"\"stats\" <first molecules file> <second molecules file> <output file> <EEM bonds>\"");
            Console.WriteLine("<method> - Method name. Only \"eem\", \"mgc\" or \"ogc\" supported.");
            Console.WriteLine("<molecules file> - Path and file name of the set of molecules.");
            Console.WriteLine("<output file> - Desired path and file name for output file.");
            Console.WriteLine("<parameters file> - Required only for EEM.");
            Console.WriteLine("<EEM bonds> - \"y\" or \"n\" to specify whether stats should be generated for atom types according to their bond types or not.");
        }

        private static bool CanParseArguments(string[] items)
        {
            if (items.Length == 0)
            {
                PrintHelp();
                return false;
            }
            else if (!(items[0].Equals("eem") || items[0].Equals("mgc") || items[0].Equals("ogc") || items[0].Equals("stats")))
            {
                Console.WriteLine("Incorrect method name. Only \"eem\", \"mgc\", \"ogc\" or \"stats\" supported.");
                PrintHelp();
                return false;
            }      
            else if(items[0].Equals("stats") && items.Length != 5)
            {
                Console.WriteLine("Incorrect number of arguments for statistics generation.");
                PrintHelp();
                return false;
            }  
            else if(((items[0].Equals("mgc") || items[0].Equals("ogc")) && items.Length != 3) || (items[0].Equals("eem") && items.Length != 4))
            {
                Console.WriteLine("Incorrect number of arguments for this method.");
                PrintHelp();
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
                                paramBonds = false;
                            else
                                paramBonds = true;
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
        
        private static void LoadMoleculesV2000(string fileName)
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
                            molecule.Code = line;
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
                            atom.Symbol = line.Substring(31, 3).Trim();
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
                            foreach (var atom in molecule.Atoms)
                                atom.HighestBondType = GetHighestBondType(molecule, atom);
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

        private static void LoadMoleculesV3000(string fileName)
        {
            molecules = new List<Molecule>();
            try
            {
                using (StreamReader reader = File.OpenText(fileName))
                {
                    string line;
                    string code = "START";
                    Molecule molecule = new Molecule();

                    while ((line = reader.ReadLine()) != null)
                    {                    
                        if (code.Equals("PARSE"))
                        {
                            if (line.Contains("COUNTS"))
                                code = "COUNTS";
                            if (line.Contains("BEGIN ATOM"))
                            {
                                code = "ATOMS";
                                continue;
                            }        
                            if(line.Contains("BEGIN BOND"))
                            {
                                code = "BONDS";
                                continue;
                            }
                            if(line.Contains("$$$$"))
                            {
                                code = "START";
                                continue;
                            }
                        }
                        if (code.Equals("START"))
                        {
                            molecule.Code = line;
                            code = "PARSE";
                        } 
                        if(code.Equals("COUNTS"))                        
                        {
                            string[] items = line.Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                            molecule.NumOfAtoms = int.Parse(items[3]);
                            molecule.NumOfBonds = int.Parse(items[4]);
                            code = "PARSE";
                        }
                        if (code.Equals("ATOMS"))
                        {
                            if(line.Contains("END ATOM"))
                            {
                                code = "PARSE";
                                continue;
                            }
                            string[] items = line.Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                            Atom atom = new Atom();
                            atom.ID = int.Parse(items[2]);
                            atom.Symbol = items[3].Trim();
                            atom.X = double.Parse(items[4], CultureInfo.InvariantCulture);
                            atom.Y = double.Parse(items[5], CultureInfo.InvariantCulture);
                            atom.Z = double.Parse(items[6], CultureInfo.InvariantCulture);
                            molecule.Atoms.Add(atom);
                        }
                        if (code.Equals("BONDS"))
                        {
                            if(line.Contains("END BOND"))
                            {
                                foreach (var atom in molecule.Atoms)
                                    atom.HighestBondType = GetHighestBondType(molecule, atom);
                                molecules.Add(molecule);
                                molecule = new Molecule();
                                code = "PARSE";
                                continue;
                            }                                
                            string[] items = line.Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                            var firstAtom = molecule.Atoms.Find(a => a.ID == int.Parse(items[4]));
                            firstAtom.Bonds.Add(int.Parse(items[5]), int.Parse(items[3]));
                        }
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
                        {
                            if (items[1].Equals("no data"))
                                throw new UnsupportedElementException($"Element {atom.Symbol} does not have an electronegativity value.");
                            return double.Parse(items[1], CultureInfo.InvariantCulture);
                        }
                            
                    }
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine("Could not read Element EN .csv file. Exception: " + ex.Message);
            }
            return 0;
        }        
                
        private static int GetHighestBondType(Molecule m, Atom a)
        {
            int bond = 0;

            foreach(var b in a.Bonds)
                if (b.Value > bond)
                    bond = b.Value;

            foreach(var atom in m.Atoms)
                foreach(var b in atom.Bonds)
                    if (b.Key == a.ID && b.Value > bond)
                        bond = b.Value;

            return bond;
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
                        int bondType = 0;
                        if (paramBonds)
                            bondType = a.HighestBondType;
                        var pars = elemParams.FirstOrDefault(x =>
                            x.ElementName.Equals(a.Symbol) &&
                            x.BondType == bondType);
                        if (pars == null)
                            throw new UnsupportedElementException($"No parameters found for {a.Symbol}.");
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
                int bondType = 0;
                if (paramBonds)
                    bondType = a.HighestBondType;
                var pars = elemParams.FirstOrDefault(x =>
                    x.ElementName.Equals(a.Symbol) &&
                    x.BondType == bondType);
                if (pars == null)
                    throw new UnsupportedElementException($"No parameters found for {a.Symbol}.");
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
                        try
                        {
                            var matrix = BuildEEMMatrix(molecule);
                            var vector = BuildEEMVector(molecule);
                            var results = matrix.Solve(vector);
                            SaveResults(file, molecule, results);
                        }
                        catch (UnsupportedElementException ex)
                        {
                            SaveIncorrectResults(file, molecule);
                        }
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
            double[,] arr = new double[m.Atoms.Count, m.Atoms.Count];

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
            double[,] arr = new double[m.Atoms.Count, m.Atoms.Count];

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
            double[,] arr = new double[m.Atoms.Count, m.Atoms.Count];

            for(int i = 0; i < m.Atoms.Count; i++)
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
            double[] arr = new double[m.Atoms.Count];

            foreach (var a in m.Atoms)
            {
                var bond = a.OrbitalBonds.First(x => !x.Value.Equals("x")).Value;
                var en = a.OrbitalENs[bond];
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
                        try
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
                        catch (UnsupportedElementException ex)
                        {
                            SaveIncorrectResults(file, molecule);
                        }
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
                        try
                        {
                            var m = BuildOGCMolecule(molecule);
                            var a = BuildConnectivityMatrixOGC(m);
                            var d = BuildDegreeMatrixOGC(m);
                            var i = BuildIdentityMatrix(m);
                            var vector = BuildOGCVector(m);
                            Matrix<double> s = d - a + i;
                            var equalizedEN = s.Solve(vector);
                            var deltaEN = equalizedEN - vector;

                            double dividend = 0;
                            double divisor = 0;
                            int n = 0;
                            foreach (var orbital in m.Atoms)
                            {
                                var bond = orbital.OrbitalBonds.First(x => !x.Value.Equals("x")).Value;
                                var orbHardness = orbital.OrbitalHardnesses.First(x => x.Key.Equals(bond)).Value;
                                dividend += deltaEN[n] / orbHardness;
                                divisor += (deltaEN[n] * deltaEN[n]) / (Math.Pow(orbHardness, 3) * GetCovalentRadius(orbital.Symbol.Trim() + orbital.HighestBondType.ToString()));
                                n++;
                            }
                            double dm = dividend / divisor;

                            n = 0;
                            foreach (var orbital in m.Atoms)
                            {
                                var bond = orbital.OrbitalBonds.First(x => !x.Value.Equals("x")).Value;
                                var orbHardness = orbital.OrbitalHardnesses.First(x => x.Key.Equals(bond)).Value;
                                orbital.OrbitalCharge = (deltaEN[n] / orbHardness) - (((deltaEN[n] * deltaEN[n]) * dm) / (Math.Pow(orbHardness, 3) * GetCovalentRadius(orbital.Symbol.Trim() + orbital.HighestBondType.ToString())));
                                n++;
                            }

                            double[] arr = new double[m.NumOfAtoms];
                            foreach (var atom in m.Atoms)
                            {
                                var orbitals = m.Atoms.FindAll(o => o.ID == atom.ID);
                                double x = 0;
                                foreach (var val in orbitals)
                                    x += val.OrbitalCharge;
                                arr[atom.ID - 1] = x;                            
                            }
                            Vector<double> results = Vector<double>.Build.Dense(arr);

                            SaveResults(file, molecule, results);
                        }catch(UnsupportedElementException ex)
                        {
                            SaveIncorrectResults(file, molecule);
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine(ex);
                //Console.WriteLine("Could not save results. Exception: " + ex.Message);
            }
        }

        private static double GetCovalentRadius(string symbol)
        {
            try
            {
                using (StreamReader reader = File.OpenText(@"..\..\Tables\CovalentRadius.csv"))
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
                        if (items[0].Equals(symbol))
                            return double.Parse(items[1], CultureInfo.InvariantCulture);
                    }
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine("Could not get orbital charges. Exception: " + ex.Message);
            }
            throw new UnsupportedElementException($"Could not find covalent radius values for {symbol}.");
        }

        private static void SaveResults(StreamWriter file, Molecule molecule, Vector<double> results)
        {
            file.WriteLine($"{molecule.Code}");
            file.WriteLine(molecule.NumOfAtoms);
            int count = 0;
            foreach (double charge in results)
            {
                if (count != molecule.NumOfAtoms)
                {
                    file.Write("{0,-4}", count + 1);
                    file.Write("{0,-3}", molecule.Atoms[count].Symbol);
                    file.Write("{0,9:F6}", charge);
                    file.WriteLine("{0,2}", molecule.Atoms[count].HighestBondType);
                    count++;
                }
            }
            file.WriteLine("$$$$");
        }

        private static void SaveIncorrectResults(StreamWriter file, Molecule molecule)
        {
            file.WriteLine($"{molecule.Code}");
            file.WriteLine("X");
            file.WriteLine("Molecule contains unsupported elements, cannot compute charges.");
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
                        if (lineNum < 1)
                            lineNum++;
                        else if (lineNum == 1)
                        {
                            molecule = new Molecule();                            
                            molecule.Code = line;
                            lineNum++;
                        }
                        else if (lineNum == 2)
                        {
                            if (line.Equals("X"))
                                lineNum = -1;
                            else
                            {
                                molecule.NumOfAtoms = int.Parse(line);
                                lineNum++;
                            }                            
                        }
                        else if (lineNum <= molecule.NumOfAtoms + 2)
                        {
                            Atom atom = new Atom();
                            atom.ID = int.Parse(line.Substring(0, 4));
                            atom.Symbol = line.Substring(4, 3);
                            atom.Charge = double.Parse(line.Substring(7, 9));
                            atom.HighestBondType = int.Parse(line.Substring(17));
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
                                
                //TODO - should this be here?
                //make sure sets contain same molecules in case of different numbers of molecules due to i.e. computing EEM charges with bonds and without (with bonds there are usually less results)
                new List<Molecule>(firstSet).ForEach(m => { if (secondSet.Find(x => x.Code == m.Code) == null) { firstSet.Remove(m); } });
                new List<Molecule>(secondSet).ForEach( m => { if (firstSet.Find(x => x.Code == m.Code) == null) { secondSet.Remove(m); } } );
                
                using (StreamWriter file = new StreamWriter(outputFilePath))
                {
                    double avg_d_avg = 0;
                    double avg_d_max = 0;
                    double avg_rmsd = 0;
                    double avg_pearson = 0;

                    Dictionary<string, List<double>> xValues = new Dictionary<string, List<double>>() { { "molecule", new List<double>() } };
                    Dictionary<string, List<double>> yValues = new Dictionary<string, List<double>>() { { "molecule", new List<double>() } };
                    
                    for (int i = 0; i < firstSet.Count; i++)
                    {
                        double d_max = 0;
                        double d_avg = 0;
                        double sum_avg = 0;
                        double rmsd = 0;
                        double sum_rmsd = 0;
                        double pearson = 0;

                        xValues["molecule"].RemoveRange(0, xValues["molecule"].Count);
                        yValues["molecule"].RemoveRange(0, yValues["molecule"].Count);

                        for (int j = 0; j < firstSet.ElementAt(i).NumOfAtoms; j++)
                        {
                            double x = firstSet.ElementAt(i).Atoms.ElementAt(j).Charge;
                            double y = secondSet.ElementAt(i).Atoms.ElementAt(j).Charge;

                            string symbol = firstSet.ElementAt(i).Atoms.ElementAt(j).Symbol.Trim();
                            if (paramBonds)
                                symbol += firstSet.ElementAt(i).Atoms.ElementAt(j).HighestBondType.ToString();

                            if (!xValues.ContainsKey(symbol))
                            {
                                xValues.Add(symbol, new List<double>());
                                yValues.Add(symbol, new List<double>());                                
                            }
                            xValues[symbol].Add(x);
                            yValues[symbol].Add(y);
                            xValues["molecule"].Add(x);
                            yValues["molecule"].Add(y);

                            if (Math.Abs(x - y) > d_max)
                                d_max = Math.Abs(x - y);

                            sum_avg += Math.Abs(x - y);
                            sum_rmsd += (x - y) * (x - y);
                        }

                        d_avg = sum_avg / (double)firstSet.ElementAt(i).NumOfAtoms;
                        rmsd = Math.Sqrt(sum_rmsd / (double)firstSet.ElementAt(i).NumOfAtoms);
                        pearson = MathNet.Numerics.Statistics.Correlation.Pearson(xValues["molecule"], yValues["molecule"]);

                        avg_d_avg += d_avg;
                        avg_d_max += d_max;
                        avg_rmsd += rmsd;
                        avg_pearson += pearson;

                        file.WriteLine(firstSet.ElementAt(i).Code);
                        file.WriteLine($"Largest absolute difference: {d_max}");
                        file.WriteLine($"Average absolute difference: {d_avg}");
                        file.WriteLine($"Root-mean-square deviation: {rmsd}");
                        file.WriteLine($"Pearson correlation coefficient: {pearson}");
                        file.WriteLine("$$$$");
                    }

                    string avgStatsFilePath = Path.Combine(Path.GetDirectoryName(outputFilePath), $"{Path.GetFileNameWithoutExtension(outputFilePath)}AvgStats.txt");
                    using (StreamWriter avgStatsFile = new StreamWriter(avgStatsFilePath))
                    {
                        avgStatsFile.WriteLine("Average values for the entire set of molecules:");
                        avgStatsFile.WriteLine($"Largest absolute difference: {avg_d_max / (double)firstSet.Count}");
                        avgStatsFile.WriteLine($"Average absolute difference: {avg_d_avg / (double)firstSet.Count}");
                        avgStatsFile.WriteLine($"Root-mean-square deviation: {avg_rmsd / (double)firstSet.Count}");
                        avgStatsFile.WriteLine($"Pearson correlation coefficient: {avg_pearson / (double)firstSet.Count}");
                        avgStatsFile.WriteLine();

                        foreach (var atomType in xValues.Keys.ToList())
                        {
                            if (atomType.Equals("molecule"))
                                continue;

                            double d_max = 0;
                            double d_avg = 0;
                            double sum_avg = 0;
                            double rmsd = 0;
                            double sum_rmsd = 0;
                            double pearson = 0;

                            for (int i = 0; i < xValues[atomType].Count; i++)
                            {
                                double x = xValues[atomType][i];
                                double y = yValues[atomType][i];

                                if (Math.Abs(x - y) > d_max)
                                    d_max = Math.Abs(x - y);

                                sum_avg += Math.Abs(x - y);
                                sum_rmsd += (x - y) * (x - y);
                            }

                            d_avg = sum_avg / (double)xValues[atomType].Count;
                            rmsd = Math.Sqrt(sum_rmsd / (double)xValues[atomType].Count);
                            pearson = MathNet.Numerics.Statistics.Correlation.Pearson(xValues[atomType], yValues[atomType]);

                            avgStatsFile.WriteLine($"Values for {atomType}:");
                            avgStatsFile.WriteLine($"Largest absolute difference: {d_max}");
                            avgStatsFile.WriteLine($"Average absolute difference: {d_avg}");
                            avgStatsFile.WriteLine($"Root-mean-square deviation: {rmsd}");
                            avgStatsFile.WriteLine($"Pearson correlation coefficient: {pearson}");
                            avgStatsFile.WriteLine(xValues[atomType].Count);
                        }
                    }

                    List<string> symbols = xValues.Keys.ToList();
                    symbols.Remove("molecule");
                    GNUPlot(firstSet, secondSet, outputFilePath, symbols);
                }
            }
            catch (Exception ex)
            {                
                Console.WriteLine("Could not save results. Exception: " + ex.Message);
            }
        }

        private static void GNUPlot(List<Molecule> firstSet, List<Molecule> secondSet, string outputFilePath, List<string> symbols)
        {
            Process plotProcess = new Process();
            plotProcess.StartInfo.FileName = Path.GetFullPath(Path.Combine(Environment.CurrentDirectory, @"..\..\gnuplot\bin\gnuplot.exe"));
            plotProcess.StartInfo.RedirectStandardInput = true;
            plotProcess.StartInfo.UseShellExecute = false;
            plotProcess.Start();

            NumberFormatInfo nfi = new NumberFormatInfo();
            nfi.NumberDecimalSeparator = ".";

            string inputFilePath = Path.Combine(Path.GetDirectoryName(outputFilePath), $"{Path.GetFileNameWithoutExtension(outputFilePath)}GNUPlotInputFile.txt");
            string graphFilePath = Path.Combine(Path.GetDirectoryName(outputFilePath), $"{Path.GetFileNameWithoutExtension(outputFilePath)}Graph.png");            

            using (StreamWriter inputFile = new StreamWriter(inputFilePath))
                for (int i = 0; i < firstSet.Count; i++)
                    for (int j = 0; j < firstSet.ElementAt(i).NumOfAtoms; j++)
                    {
                        string a = firstSet.ElementAt(i).Atoms.ElementAt(j).Charge.ToString(nfi);
                        string b = secondSet.ElementAt(i).Atoms.ElementAt(j).Charge.ToString(nfi);
                        var atom = firstSet.ElementAt(i).Atoms.ElementAt(j);
                        int c;
                        if (paramBonds)
                            c = symbols.IndexOf(atom.Symbol.Trim() + atom.HighestBondType.ToString());
                        else
                            c = symbols.IndexOf(atom.Symbol.Trim());
                        inputFile.WriteLine($"{a} {b} {c}");
                    }

            using (StreamWriter gnuplotFile = plotProcess.StandardInput)
            {
                gnuplotFile.WriteLine("set terminal png");
                gnuplotFile.WriteLine("set title \"Partial atomic charge correlation graph\"");
                gnuplotFile.WriteLine("set xlabel \"First set of atom charges\"");
                gnuplotFile.WriteLine("set ylabel \"Second set of atom charges\"");
                gnuplotFile.WriteLine($"set output '{graphFilePath}'");
                gnuplotFile.WriteLine($"set palette model RGB");
                gnuplotFile.Write("set palette defined(");
                int i = 0;
                foreach(var s in symbols)
                {
                    if(i == symbols.Count-1)
                        gnuplotFile.WriteLine($"{i} \"{GetAtomColorName(i)}\")");
                    else
                        gnuplotFile.Write($"{i} \"{GetAtomColorName(i)}\", ");
                    i++;
                }
                gnuplotFile.WriteLine($"set cbrange [0:{symbols.Count-1}]");
                double pos = 1.0;
                i = 0;
                foreach (var s in symbols)
                {
                    pos -= 0.05;
                    gnuplotFile.WriteLine($"set label \"{s}\" at graph 0.05,{pos.ToString(nfi)} tc \"{GetAtomColorName(i)}\"");
                    i++;
                }                
                gnuplotFile.WriteLine("unset colorbox");
                gnuplotFile.WriteLine($"plot '{inputFilePath}' u 1:2:3 with points palette notitle");
            }
        }

        private static string GetAtomColorName(int n)
        {
            switch (n)
            {
                case 0: return "black";
                case 1: return "red";
                case 2: return "green";
                case 3: return "blue";
                case 4: return "purple";
                case 5: return "orange";
                case 6: return "yellow";
                case 7: return "greenyellow";
                case 8: return "aquamarine";
                case 9: return "cyan";
                case 10: return "magenta";
                case 11: return "aqua";
                default: return "error";
            }
        }

        private static void SetOrbitalENsAndHardnesses(Molecule molecule)
        {
            int bondType = 1;
            foreach (var atom in molecule.Atoms)
            {
                switch (atom.Symbol)
                {
                    case "H":
                        atom.OrbitalENs.Add("s", GetOrbitalEN("H", "s_"));
                        atom.OrbitalHardnesses.Add("s", GetOrbitalHardness("H", "s_"));
                        break;
                    case "C":
                        bondType = atom.HighestBondType;
                        switch (bondType)
                        {
                            case 1:                                
                                atom.OrbitalENs.Add("s", GetOrbitalEN("C", "te_ te te te"));
                                atom.OrbitalHardnesses.Add("s", GetOrbitalHardness("C", "te_ te te te"));
                                break;
                            case 2:
                                atom.OrbitalENs.Add("s", GetOrbitalEN("C", "tr_ tr tr pp"));
                                atom.OrbitalENs.Add("p", GetOrbitalEN("C", "tr tr tr pp_"));
                                atom.OrbitalHardnesses.Add("s", GetOrbitalHardness("C", "tr_ tr tr pp"));
                                atom.OrbitalHardnesses.Add("p", GetOrbitalHardness("C", "tr tr tr pp_"));
                                break;
                            case 3:
                                atom.OrbitalENs.Add("s", GetOrbitalEN("C", "di_ di pp pp"));
                                atom.OrbitalENs.Add("p", GetOrbitalEN("C", "di di pp_ pp"));
                                atom.OrbitalHardnesses.Add("s", GetOrbitalHardness("C", "di_ di pp pp"));
                                atom.OrbitalHardnesses.Add("p", GetOrbitalHardness("C", "di di pp_ pp"));
                                break;
                            default: break;
                        }
                        break;
                    case "N":
                        bondType = atom.HighestBondType;
                        switch (bondType)
                        {
                            case 1:                                
                                atom.OrbitalENs.Add("s", GetOrbitalEN("N", "te2 te_ te te"));
                                atom.OrbitalENs.Add("n", GetOrbitalEN("N", "te2_ te te te"));
                                atom.OrbitalHardnesses.Add("s", GetOrbitalHardness("N", "te2 te_ te te"));
                                atom.OrbitalHardnesses.Add("n", GetOrbitalHardness("N", "te2_ te te te"));
                                break;
                            case 2:
                                atom.OrbitalENs.Add("s", GetOrbitalEN("N", "tr2 tr_ tr pp"));
                                atom.OrbitalENs.Add("p", GetOrbitalEN("N", "tr2 tr tr pp_"));                                
                                atom.OrbitalENs.Add("n", GetOrbitalEN("N", "tr2_ tr tr pp"));
                                atom.OrbitalHardnesses.Add("s", GetOrbitalHardness("N", "tr2 tr_ tr pp"));
                                atom.OrbitalHardnesses.Add("p", GetOrbitalHardness("N", "tr2 tr tr pp_"));
                                atom.OrbitalHardnesses.Add("n", GetOrbitalHardness("N", "tr2_ tr tr pp"));
                                break;
                            case 3:
                                atom.OrbitalENs.Add("s", GetOrbitalEN("N", "di2 di_ pp pp"));
                                atom.OrbitalENs.Add("p", GetOrbitalEN("N", "di2 di pp_ pp"));
                                atom.OrbitalENs.Add("n", GetOrbitalEN("N", "di2_ di pp pp"));
                                atom.OrbitalHardnesses.Add("s", GetOrbitalHardness("N", "di2 di_ pp pp"));
                                atom.OrbitalHardnesses.Add("p", GetOrbitalHardness("N", "di2 di pp_ pp"));
                                atom.OrbitalHardnesses.Add("n", GetOrbitalHardness("N", "di2_ di pp pp"));
                                break;
                            default: break;
                        }
                        break;
                    case "O":
                        bondType = atom.HighestBondType;
                        switch (bondType)
                        {
                            case 1:
                                atom.OrbitalENs.Add("s", GetOrbitalEN("O", "te2 te2 te_ te"));
                                atom.OrbitalENs.Add("n", GetOrbitalEN("O", "te2_ te2 te te"));
                                atom.OrbitalHardnesses.Add("s", GetOrbitalHardness("O", "te2 te2 te_ te"));
                                atom.OrbitalHardnesses.Add("n", GetOrbitalHardness("O", "te2_ te2 te te"));
                                break;
                            case 2:
                                atom.OrbitalENs.Add("s", GetOrbitalEN("O", "tr2 tr2 tr_ pp"));
                                atom.OrbitalENs.Add("p", GetOrbitalEN("O", "tr2 tr2 tr pp_"));
                                atom.OrbitalENs.Add("n", GetOrbitalEN("O", "tr2_ tr2 tr pp"));
                                atom.OrbitalHardnesses.Add("s", GetOrbitalHardness("O", "tr2 tr2 tr_ pp"));
                                atom.OrbitalHardnesses.Add("p", GetOrbitalHardness("O", "tr2 tr2 tr pp_"));
                                atom.OrbitalHardnesses.Add("n", GetOrbitalHardness("O", "tr2_ tr2 tr pp"));
                                break;
                            default: break;
                        }
                        break;
                    case "S":
                        bondType = atom.HighestBondType;
                        switch (bondType)
                        {
                            case 1:
                                atom.OrbitalENs.Add("s", GetOrbitalEN("S", "te2 te2 te_ te"));
                                atom.OrbitalENs.Add("n", GetOrbitalEN("S", "te2_ te2 te te"));
                                atom.OrbitalHardnesses.Add("s", GetOrbitalHardness("S", "te2 te2 te_ te"));
                                atom.OrbitalHardnesses.Add("n", GetOrbitalHardness("S", "te2_ te2 te te"));
                                break;
                            case 2:
                                atom.OrbitalENs.Add("s", GetOrbitalEN("S", "tr2 tr2 tr_ pp"));
                                atom.OrbitalENs.Add("p", GetOrbitalEN("S", "tr2 tr2 tr pp_"));
                                atom.OrbitalENs.Add("n", GetOrbitalEN("S", "tr2_ tr2 tr pp"));
                                atom.OrbitalHardnesses.Add("s", GetOrbitalHardness("S", "tr2 tr2 tr_ pp"));
                                atom.OrbitalHardnesses.Add("p", GetOrbitalHardness("S", "tr2 tr2 tr pp_"));
                                atom.OrbitalHardnesses.Add("n", GetOrbitalHardness("S", "tr2_ tr2 tr pp"));
                                break;
                            default: break;
                        }
                        break;
                    default: throw new UnsupportedElementException($"Cannot set Orbital EN or hardness for {atom.Symbol + atom.HighestBondType.ToString()}.");
                }
            }
        }

        private static double GetOrbitalEN(string symbol, string state)
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
            throw new UnsupportedElementException($"Could not find orbital electronegativity for {symbol}.");
        }

        private static double GetOrbitalHardness(string symbol, string state)
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
                            return double.Parse(items[3], CultureInfo.InvariantCulture);
                    }
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine("Could not get orbital hardnesses. Exception: " + ex.Message);
            }
            throw new UnsupportedElementException($"Could not find orbital hardness for {symbol}.");
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
            SetOrbitalENsAndHardnesses(m);
            Molecule ogcMolecule = new Molecule();            
            ogcMolecule.Code = m.Code;
            ogcMolecule.NumOfAtoms = m.NumOfAtoms;
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
                    atom.OrbitalENs = a.OrbitalENs;
                    atom.OrbitalHardnesses = a.OrbitalHardnesses;
                    atom.OrbitalID = $"{atom.ID}x{i}";
                    atom.HighestBondType = a.HighestBondType;
                    ogcMolecule.Atoms.Add(atom);
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

                if (a.Symbol.Trim().Equals("O") || a.Symbol.Trim().Equals("S"))
                {
                    for(int i = 0; i < 2; i++)
                    {
                        var firstOrbital = GetFreeOrbital(ogcMolecule, a.ID);
                        firstOrbital.OrbitalBonds.Add(firstOrbital.ID.ToString(), "n");
                    }                    
                }

                foreach (var b in a.Bonds)
                {
                    var freeOrbital = GetFreeOrbital(ogcMolecule, b.Key);
                    if(freeOrbital == null)
                        continue;
                    if(b.Value == 1)
                    {
                        var sOrbital = GetFreeOrbital(ogcMolecule, a.ID);
                        freeOrbital.OrbitalBonds.Add(sOrbital.OrbitalID, "s");
                        sOrbital.OrbitalBonds.Add(freeOrbital.OrbitalID, "s");
                    }
                    if(b.Value == 2)
                    {
                        var pOrbital = GetFreeOrbital(ogcMolecule, a.ID);
                        freeOrbital.OrbitalBonds.Add(pOrbital.OrbitalID, "p");
                        pOrbital.OrbitalBonds.Add(freeOrbital.OrbitalID, "p");
                        freeOrbital = GetFreeOrbital(ogcMolecule, b.Key);
                        var sOrbital = GetFreeOrbital(ogcMolecule, a.ID);
                        freeOrbital.OrbitalBonds.Add(sOrbital.OrbitalID, "s");
                        sOrbital.OrbitalBonds.Add(freeOrbital.OrbitalID, "s");
                    }
                    if(b.Value == 3)
                    {
                        var pOrbital = GetFreeOrbital(ogcMolecule, a.ID);
                        freeOrbital.OrbitalBonds.Add(pOrbital.OrbitalID, "p");
                        pOrbital.OrbitalBonds.Add(freeOrbital.OrbitalID, "p");
                        freeOrbital = GetFreeOrbital(ogcMolecule, b.Key);
                        var sOrbital = GetFreeOrbital(ogcMolecule, a.ID);
                        freeOrbital.OrbitalBonds.Add(sOrbital.OrbitalID, "s");
                        sOrbital.OrbitalBonds.Add(freeOrbital.OrbitalID, "s");
                        freeOrbital = GetFreeOrbital(ogcMolecule, b.Key);
                        var p2Orbital = GetFreeOrbital(ogcMolecule, a.ID);
                        freeOrbital.OrbitalBonds.Add(p2Orbital.OrbitalID, "p");
                        p2Orbital.OrbitalBonds.Add(freeOrbital.OrbitalID, "p");
                    }
                }
            }

            return ogcMolecule;
        }
    }
}
