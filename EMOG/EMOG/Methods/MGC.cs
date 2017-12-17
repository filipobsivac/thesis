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
    class MGC
    {
        private List<Molecule> molecules;
        private const string ElementENTablePath = @"..\..\Data\ElementEN.csv";

        private double GetElectronegativity(Atom atom)
        {
            try
            {
                using (StreamReader reader = File.OpenText(ElementENTablePath))
                {
                    string line;
                    while ((line = reader.ReadLine()) != null)
                    {
                        string[] items = line.Split(';');
                        if (items[0].Equals(atom.Symbol))
                        {
                            if (items[1].Equals("no data"))
                                throw new UnsupportedAtomException($"Element {atom.Symbol} does not have an EN value in database.");
                            return double.Parse(items[1], CultureInfo.InvariantCulture);
                        }

                    }
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine("Could not get element EN. Exception: " + ex.Message);
            }
            throw new UnsupportedAtomException($"Could not find element {atom.Symbol} in EN database.");
        }

        private Matrix<double> BuildDegreeMatrix(Molecule m)
        {
            double[,] arr = new double[m.NumOfAtoms, m.NumOfAtoms];

            foreach (var a in m.Atoms)
            {
                foreach (var b in a.Bonds)
                {
                    arr[b.Key - 1, b.Key - 1] += b.Value;
                    arr[a.ID - 1, a.ID - 1] += b.Value;
                }
            }

            return Matrix<double>.Build.DenseOfArray(arr);
        }

        private Matrix<double> BuildConnectivityMatrix(Molecule m)
        {
            double[,] arr = new double[m.NumOfAtoms, m.NumOfAtoms];

            for (int i = 0; i < m.NumOfAtoms; i++)
            {
                foreach (var b in m.Atoms[i].Bonds)
                {
                    arr[i, b.Key - 1] = b.Value;
                    arr[b.Key - 1, i] = b.Value;
                }
            }

            return Matrix<double>.Build.DenseOfArray(arr);
        }

        private Matrix<double> BuildIdentityMatrix(Molecule m)
        {
            double[,] arr = new double[m.Atoms.Count, m.Atoms.Count];

            for (int i = 0; i < m.Atoms.Count; i++)
                arr[i, i] = 1;

            return Matrix<double>.Build.DenseOfArray(arr);
        }

        private Vector<double> BuildMGCVector(Molecule m)
        {
            double[] arr = new double[m.NumOfAtoms];

            for (int i = 0; i < m.NumOfAtoms; i++)
                arr[i] = GetElectronegativity(m.Atoms[i]);

            return Vector<double>.Build.Dense(arr);
        }

        public void SolveMGC(string outputFilePath, string moleculesFilePath)
        {
            try
            {
                moleculesFilePath = Path.GetFullPath(Path.Combine(Environment.CurrentDirectory, moleculesFilePath));
                molecules = LoadMolecules(moleculesFilePath);
            }
            catch(Exception ex)
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
                            var a = BuildConnectivityMatrix(molecule);
                            var d = BuildDegreeMatrix(molecule);
                            var i = BuildIdentityMatrix(molecule);
                            var vector = BuildMGCVector(molecule);
                            double avgEN = MathNet.Numerics.Statistics.Statistics.Mean(vector);
                            Matrix<double> s = d - a + i;
                            var x = s.Solve(vector);
                            var results = (x - vector) * (1.0 / avgEN);

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
