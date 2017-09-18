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
            #region Data loading test
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
            Console.WriteLine();
            #endregion

            #region Distance calculating test
            Console.WriteLine("Distance between the first and second atoms:\n");
            Console.WriteLine("First atom coordinates:");
            Console.WriteLine("X\tY\tZ");
            Console.WriteLine($"{atom.X}\t{atom.Y}\t{atom.Z}\n");
            Console.WriteLine("Second atom coordinates:");
            Console.WriteLine("X\tY\tZ");
            var atom2 = m.Atoms.ElementAt(1);
            Console.WriteLine($"{atom2.X}\t{atom2.Y}\t{atom2.Z}\n");
            double distance = CalculateDistance(atom, atom2);
            Console.WriteLine($"Distance between the two atoms: {distance}\n");
            #endregion

            #region First matrix solution method test
            Console.WriteLine("First matrix solution method:\n");
            Console.WriteLine("EEM Matrix:\n");
            double[,] matrix = BuildEEMMatrix(m);
            for(int i = 0; i <= m.NumOfAtoms; i++)
            {
                for(int j = 0; j <= m.NumOfAtoms + 2; j++)
                {
                    Console.Write($"{Math.Round(matrix[i,j], 2)}\t");
                }
                Console.WriteLine();
            }            
            Console.WriteLine(ResolveMatrix(matrix, m.NumOfAtoms + 1, m.NumOfAtoms + 1));
            Console.WriteLine();
            #endregion

            #region Second matrix solution method test
            Console.WriteLine("Second matrix solution method:\n");
            Console.WriteLine("EEM Matrix:\n");
            double[,] matrix2 = BuildEEMMatrix2(m);
            for (int i = 0; i <= m.NumOfAtoms; i++)
            {
                for (int j = 0; j <= m.NumOfAtoms + 1; j++)
                {
                    Console.Write($"{Math.Round(matrix2[i, j], 2)}\t");
                }
                Console.WriteLine();
            }
            Console.WriteLine();
            Console.WriteLine(ResolveMatrix2(matrix2));
            Console.WriteLine();
            #endregion

            #region Third matrix solution method test
            Console.WriteLine("Third matrix solution method:\n");
            Console.WriteLine("EEM Matrix:\n");
            double[,] matrix3 = BuildEEMMatrix3(m);
            for (int i = 0; i <= m.NumOfAtoms; i++)
            {
                for (int j = 0; j <= m.NumOfAtoms; j++)
                {
                    Console.Write($"{Math.Round(matrix3[i, j], 2)}\t");
                }
                Console.WriteLine();
            }
            Console.WriteLine();
            var A = Matrix<double>.Build.DenseOfArray(matrix3);
            var B = Vector<double>.Build.Dense(BuildEEMVector(m));
            var X = A.Solve(B);
            int count = 0;
            foreach(double a in X)
            {
                if(count == m.NumOfAtoms)
                {
                    Console.Write($"X = {a}\n");
                }else
                {
                    Console.Write($"x{count} = {a}\n");
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
                        line.Substring(31, 3).Trim(' '),
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

        //TODO Q value (perhaps X value as well?)
        private static double[,] BuildEEMMatrix(Molecule m)
        {
            int columns = m.NumOfAtoms + 3;
            int rows = m.NumOfAtoms + 1;
            double[,] arr = new double[rows, columns];
            var atom = m.Atoms.ElementAt(0);
            var param = elemParams.First(x => x.ElementName.Equals(atom.Symbol));
            double k = param.Kappa;

            for (int i = 0; i < rows - 1; i++)
            {
                for (int j = 0; j < columns - 1; j++)
                {
                    if (j == columns - 2)
                    {
                        var a = m.Atoms.First(x => x.ID == i + 1);
                        var p = elemParams.First(x => x.ElementName.Equals(a.Symbol));
                        arr[i, j] = -(p.A);
                    }
                    else if (j == columns - 3)
                    {
                        arr[i, j] = -1;
                    }
                    else if(i == j)
                    {
                        var a = m.Atoms.First(x => x.ID == i + 1);
                        var p = elemParams.First(x => x.ElementName.Equals(a.Symbol));
                        arr[i, j] = p.B;
                    }else
                    {
                        var a1 = m.Atoms.First(x => x.ID == i + 1);
                        var a2 = m.Atoms.First(x => x.ID == j + 1);
                        arr[i, j] = k / (CalculateDistance(a1, a2));
                    }
                }
                arr[i, columns - 1] = 0;
            }
            for(int i = 0; i < columns - 3; i++)
            {
                arr[rows - 1, i] = 1;
            }
            arr[rows - 1, columns - 3] = 0;
            arr[rows - 1, columns - 2] = 1;//should be Q, for now 1
            arr[rows - 1, columns - 1] = 0;
            return arr;
        }

        //TODO Q value (perhaps X value as well?)
        private static double[,] BuildEEMMatrix2(Molecule m)
        {
            int columns = m.NumOfAtoms + 2;
            int rows = m.NumOfAtoms + 1;
            double[,] arr = new double[rows, columns];
            var atom = m.Atoms.ElementAt(0);
            var param = elemParams.First(x => x.ElementName.Equals(atom.Symbol));
            double k = param.Kappa;

            for (int i = 0; i < rows - 1; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    if (j == columns - 1)
                    {
                        var a = m.Atoms.First(x => x.ID == i + 1);
                        var p = elemParams.First(x => x.ElementName.Equals(a.Symbol));
                        arr[i, j] = -(p.A);
                    }
                    else if (j == columns - 2)
                    {
                        arr[i, j] = -1;
                    }
                    else if (i == j)
                    {
                        var a = m.Atoms.First(x => x.ID == i + 1);
                        var p = elemParams.First(x => x.ElementName.Equals(a.Symbol));
                        arr[i, j] = p.B;
                    }
                    else
                    {
                        var a1 = m.Atoms.First(x => x.ID == i + 1);
                        var a2 = m.Atoms.First(x => x.ID == j + 1);
                        arr[i, j] = k / (CalculateDistance(a1, a2));
                    }
                }
            }
            for (int i = 0; i < columns - 2; i++)
            {
                arr[rows - 1, i] = 1;
            }
            arr[rows - 1, columns - 2] = 0;
            arr[rows - 1, columns - 1] = 1;//should be Q, for now 1
            return arr;
        }
        
        private static double[,] BuildEEMMatrix3(Molecule m)
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
                        var p = elemParams.First(x => x.ElementName.Equals(a.Symbol));
                        arr[i, j] = p.B;
                    }
                    else if (i == rows - 1)
                    {
                        arr[i, j] = 1;
                    }
                    else
                    {
                        var a1 = m.Atoms.First(x => x.ID == i + 1);
                        var a2 = m.Atoms.First(x => x.ID == j + 1);
                        arr[i, j] = k / (CalculateDistance(a1, a2));
                    }
                }
            }
            arr[rows - 1, columns - 1] = 0;
            return arr;
        }

        //TODO Q value (perhaps X value as well?)
        private static double[] BuildEEMVector(Molecule m)
        {
            double[] vector = new double[m.NumOfAtoms + 1];
            int count = 0;
            foreach (Atom a in m.Atoms)
            {
                var p = elemParams.First(x => x.ElementName.Equals(a.Symbol));
                vector[count] = -(p.A);
                count++;
            }
            vector[count] = 1;//should be Q, for now 1
            return vector;
        }

        private static string ResolveMatrix(double[,] matrix, int num_rows, int num_cols)
        {
            const double tiny = 0.00001;
            string txt = "";

            double[,] arr = matrix;
            double[,] orig_arr = matrix;

            // Start solving.
            for (int r = 0; r < num_rows - 1; r++)
            {
                // Zero out all entries in column r after this row.
                // See if this row has a non-zero entry in column r.
                if (Math.Abs(arr[r, r]) < tiny)
                {
                    // Too close to zero. Try to swap with a later row.
                    for (int r2 = r + 1; r2 < num_rows; r2++)
                    {
                        if (Math.Abs(arr[r2, r]) > tiny)
                        {
                            // This row will work. Swap them.
                            for (int c = 0; c <= num_cols; c++)
                            {
                                double tmp = arr[r, c];
                                arr[r, c] = arr[r2, c];
                                arr[r2, c] = tmp;
                            }
                            break;
                        }
                    }
                }

                // If this row has a non-zero entry in column r, use it.
                if (Math.Abs(arr[r, r]) > tiny)
                {
                    // Zero out this column in later rows.
                    for (int r2 = r + 1; r2 < num_rows; r2++)
                    {
                        double factor = -arr[r2, r] / arr[r, r];
                        for (int c = r; c <= num_cols; c++)
                        {
                            arr[r2, c] = arr[r2, c] + factor * arr[r, c];
                        }
                    }
                }
            }

            // See if we have a solution.
            if (arr[num_rows - 1, num_cols - 1] == 0)
            {
                // We have no solution.
                // See if all of the entries in this row are 0.
                bool all_zeros = true;
                for (int c = 0; c <= num_cols + 1; c++)
                {
                    if (arr[num_rows - 1, c] != 0)
                    {
                        all_zeros = false;
                        break;
                    }
                }
                if (all_zeros)
                {
                    txt = "The solution is not unique";
                }
                else
                {
                    txt = "There is no solution";
                }
            }
            else
            {
                // Backsolve.
                for (int r = num_rows - 1; r >= 0; r--)
                {
                    double tmp = arr[r, num_cols];
                    for (int r2 = r + 1; r2 < num_rows; r2++)
                    {
                        tmp -= arr[r, r2] * arr[r2, num_cols + 1];
                    }
                    arr[r, num_cols + 1] = tmp / arr[r, r];
                }

                // Display the results.
                txt = "       Values:";
                for (int r = 0; r < num_rows; r++)
                {
                    txt += "\r\nx" + r.ToString() + " = " +
                        arr[r, num_cols + 1].ToString();
                }

                // Verify.
                txt += "\r\n    Check:";
                for (int r = 0; r < num_rows; r++)
                {
                    double tmp = 0;
                    for (int c = 0; c < num_cols; c++)
                    {
                        tmp += orig_arr[r, c] * arr[c, num_cols + 1];
                    }
                    txt += "\r\n" + tmp.ToString();
                }

                txt = txt.Substring("\r\n".Length + 1);
            }

            return txt;
        }

        private static bool ResolveMatrix2(double[,] M)
        {
            // input checks
            int rowCount = M.GetUpperBound(0) + 1;
            if (M == null || M.Length != rowCount * (rowCount + 1))
                throw new ArgumentException("The algorithm must be provided with a (n x n+1) matrix.");
            if (rowCount < 1)
                throw new ArgumentException("The matrix must at least have one row.");

            // pivoting
            for (int col = 0; col + 1 < rowCount; col++) if (M[col, col] == 0)
                // check for zero coefficients
                {
                    // find non-zero coefficient
                    int swapRow = col + 1;
                    for (; swapRow < rowCount; swapRow++) if (M[swapRow, col] != 0) break;

                    if (M[swapRow, col] != 0) // found a non-zero coefficient?
                    {
                        // yes, then swap it with the above
                        double[] tmp = new double[rowCount + 1];
                        for (int i = 0; i < rowCount + 1; i++)
                        { tmp[i] = M[swapRow, i]; M[swapRow, i] = M[col, i]; M[col, i] = tmp[i]; }
                    }
                    else return false; // no, then the matrix has no unique solution
                }

            // elimination
            for (int sourceRow = 0; sourceRow + 1 < rowCount; sourceRow++)
            {
                for (int destRow = sourceRow + 1; destRow < rowCount; destRow++)
                {
                    double df = M[sourceRow, sourceRow];
                    double sf = M[destRow, sourceRow];
                    for (int i = 0; i < rowCount + 1; i++)
                        M[destRow, i] = M[destRow, i] * df - M[sourceRow, i] * sf;
                }
            }

            // back-insertion
            for (int row = rowCount - 1; row >= 0; row--)
            {
                double f = M[row, row];
                if (f == 0) return false;

                for (int i = 0; i < rowCount + 1; i++) M[row, i] /= f;
                for (int destRow = 0; destRow < row; destRow++)
                { M[destRow, rowCount] -= M[destRow, row] * M[row, rowCount]; M[destRow, row] = 0; }
            }
            return true;
        }
    }
}
