using EMOG.Structures;
using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static EMOG.Methods.MiscFunctions;

namespace EMOG.Methods
{
    class EEM
    {
        private List<Molecule> molecules;

        private bool paramBonds;
        private double kappa;
        private List<ElementParameters> elemParams;

        private void LoadParameters(string fileName)
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
                        paramBonds = items[0].Equals("ElemBond") ? true : false;
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
        }

        private double CalculateDistance(Atom atom1, Atom atom2)
        {
            return Math.Sqrt(
                  (atom2.X - atom1.X) * (atom2.X - atom1.X)
                + (atom2.Y - atom1.Y) * (atom2.Y - atom1.Y)
                + (atom2.Z - atom1.Z) * (atom2.Z - atom1.Z)
                );
        }

        private Matrix<double> BuildEEMMatrix(Molecule m)
        {
            double[,] arr = new double[m.NumOfAtoms + 1, m.NumOfAtoms + 1];

            for (int i = 0; i < m.NumOfAtoms + 1; i++)
            {
                for (int j = 0; j < m.NumOfAtoms + 1; j++)
                {
                    if (j == m.NumOfAtoms)
                        arr[i, j] = -1;
                    else if (i == m.NumOfAtoms)
                        arr[i, j] = 1;
                    else if (i == j)
                    {
                        var a = m.Atoms.First(x => x.ID == i + 1);
                        //if parameters file has specified bond types, assign the highest one, otherwise 0
                        int bondType = paramBonds ? a.HighestBondType : 0;
                        var pars = elemParams.FirstOrDefault(x =>
                            x.ElementName.Equals(a.Symbol) &&
                            x.BondType == bondType);
                        if (pars == null)
                        {
                            string s = $"No parameters found for {a.Symbol}";
                            if (paramBonds) s += a.HighestBondType.ToString();
                            throw new UnsupportedAtomException($"{s}.");
                        }
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
            arr[m.NumOfAtoms, m.NumOfAtoms] = 0;//Q
            return Matrix<double>.Build.DenseOfArray(arr);
        }

        private Vector<double> BuildEEMVector(Molecule m)
        {
            double[] vector = new double[m.NumOfAtoms + 1];
            int count = 0;
            foreach (Atom a in m.Atoms)
            {
                //if parameters file has specified bond types, assign the highest one, otherwise 0
                int bondType = paramBonds ? a.HighestBondType : 0;
                var pars = elemParams.FirstOrDefault(x =>
                    x.ElementName.Equals(a.Symbol) &&
                    x.BondType == bondType);
                if (pars == null)
                {
                    string s = $"No parameters found for {a.Symbol}";
                    if (paramBonds) s += a.HighestBondType.ToString();
                    throw new UnsupportedAtomException($"{s}.");
                }
                vector[count] = -pars.A;
                count++;
            }
            vector[count] = 0;//Q
            return Vector<double>.Build.Dense(vector);
        }

        public void SolveEEM(string parametersFilePath, string outputFilePath, string moleculesFilePath)
        {
            try
            {
                parametersFilePath = Path.GetFullPath(Path.Combine(Environment.CurrentDirectory, parametersFilePath));
                LoadParameters(parametersFilePath);
            }catch(Exception ex)
            {
                Console.WriteLine("Could not load parameters. Exception: " + ex.Message);
                return;
            }

            try
            {
                moleculesFilePath = Path.GetFullPath(Path.Combine(Environment.CurrentDirectory, moleculesFilePath));
                molecules = LoadMolecules(moleculesFilePath);
            }catch(Exception ex)
            {
                Console.WriteLine("Could not load molecules. Exception: " + ex.Message);
                return;
            }            

            try
            {
                outputFilePath = Path.GetFullPath(Path.Combine(Environment.CurrentDirectory, outputFilePath));
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
                        catch (UnsupportedAtomException ex)
                        {
                            SaveIncorrectResults(file, molecule, ex);
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine("Could not save results. Exception: " + ex.Message);
            }
        }
    }
}
