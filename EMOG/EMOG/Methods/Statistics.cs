using EMOG.Structures;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EMOG.Methods
{
    class Statistics
    {
        private static string GNUPlotExePath;
        private bool paramBonds { get; set; }

        public void GenerateStatistics(string firstFilePath, string secondFilePath, string outputFilePath, string GnuplotExePath, string parameterBonds)
        {
            if (parameterBonds.Equals("y"))
                paramBonds = true;
            if (parameterBonds.Equals("n"))
                paramBonds = false;

            List<Molecule> firstSet;
            List<Molecule> secondSet;
            bool origParamBonds = paramBonds;
            bool firstSetParamBonds = paramBonds;
            bool secondSetParamBonds = paramBonds;
            try
            {
                firstFilePath = Path.GetFullPath(Path.Combine(Environment.CurrentDirectory, firstFilePath));
                secondFilePath = Path.GetFullPath(Path.Combine(Environment.CurrentDirectory, secondFilePath));

                origParamBonds = paramBonds;
                firstSet = LoadMoleculesFromOutputFile(firstFilePath);
                firstSetParamBonds = paramBonds;
                secondSet = LoadMoleculesFromOutputFile(secondFilePath);
                secondSetParamBonds = paramBonds;
                //if at least one set has bond types specified, use them (in case they should be used)
                if (origParamBonds && (firstSetParamBonds || secondSetParamBonds))
                    paramBonds = true;
                else
                    paramBonds = false;
            }
            catch(Exception ex)
            {
                Console.WriteLine("Could not load molecules. Exception: " + ex.Message);
                return;
            }

            NumberFormatInfo nfi = new NumberFormatInfo();
            nfi.NumberDecimalSeparator = ".";

            try
            {
                outputFilePath = Path.GetFullPath(Path.Combine(Environment.CurrentDirectory, outputFilePath));
                GNUPlotExePath = Path.GetFullPath(GnuplotExePath);

                //make sure both sets contain the same molecules in case different methods were used that do not work for the same atom types
                //remove duplicates as well
                new List<Molecule>(firstSet).ForEach(
                    m => {
                        if (secondSet.FindAll(
                     x => x.Code.Equals(m.Code)).Count != 1)
                        {
                            firstSet.RemoveAll(mol => mol.Code.Equals(m.Code));
                            secondSet.RemoveAll(mol => mol.Code.Equals(m.Code));
                        }
                    }
                );
                new List<Molecule>(secondSet).ForEach(
                    m => {
                        if (firstSet.FindAll(
                     x => x.Code.Equals(m.Code)).Count != 1)
                        {
                            secondSet.RemoveAll(mol => mol.Code.Equals(m.Code));
                            firstSet.RemoveAll(mol => mol.Code.Equals(m.Code));
                        }
                    }
                );

                //sort lists by molecule codes in case sets have molecules in different orders
                firstSet = firstSet.OrderBy(m => m.Code).ToList();
                secondSet = secondSet.OrderBy(m => m.Code).ToList();

                if (firstSet.Count == 0 || secondSet.Count == 0)
                {
                    Console.WriteLine("Cannot generate statistics - sets have no common molecules.");
                    return;
                }

                using (StreamWriter file = new StreamWriter(outputFilePath))
                {
                    file.WriteLine("Molecule;Largest absolute difference;Average absolute difference;Root-mean-square deviation;Pearson correlation coefficient;Pearson squared");

                    //variables used to calculate average values for the entire set of molecules
                    double avg_d_avg = 0;
                    double avg_d_max = 0;
                    double avg_rmsd = 0;
                    double avg_pearson = 0;
                    double avg_pearson_sq = 0;

                    //variables used to store charges for each atom type
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

                            string symbol = firstSet.ElementAt(i).Atoms.ElementAt(j).Symbol;
                            //add bond type in case it was specified to generate statistics with bond types
                            if (paramBonds)
                            {
                                var set = firstSetParamBonds ? firstSet : secondSet;
                                symbol += set.ElementAt(i).Atoms.ElementAt(j).HighestBondType.ToString();
                            }

                            //add new atom type
                            if (!xValues.ContainsKey(symbol))
                            {
                                xValues.Add(symbol, new List<double>());
                                yValues.Add(symbol, new List<double>());
                            }
                            //save values for atom type and current molecule
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
                        avg_pearson_sq += pearson * pearson;

                        file.WriteLine($"{firstSet.ElementAt(i).Code};{Math.Round(d_max, 4).ToString(nfi)};{Math.Round(d_avg, 4).ToString(nfi)};{Math.Round(rmsd, 4).ToString(nfi)};{Math.Round(pearson, 4).ToString(nfi)};{Math.Round(pearson * pearson, 4).ToString(nfi)}");
                    }

                    string avgStatsFilePath = Path.Combine(Path.GetDirectoryName(outputFilePath), $"{Path.GetFileNameWithoutExtension(outputFilePath)}AvgStats.txt");
                    using (StreamWriter avgStatsFile = new StreamWriter(avgStatsFilePath))
                    {
                        avgStatsFile.WriteLine("Average statistics for the entire set of molecules:");
                        avgStatsFile.WriteLine("D_MAX    D_AVG    RMSD     PEARSON  PEARSON^2");
                        avgStatsFile.Write(string.Format("{0,6:F4}", avg_d_max / (double)firstSet.Count).Replace(',', '.'));
                        avgStatsFile.Write(string.Format("{0,9:F4}", avg_d_avg / (double)firstSet.Count).Replace(',', '.'));
                        avgStatsFile.Write(string.Format("{0,9:F4}", avg_rmsd / (double)firstSet.Count).Replace(',', '.'));
                        avgStatsFile.Write(string.Format("{0,9:F4}", avg_pearson / (double)firstSet.Count).Replace(',', '.'));
                        avgStatsFile.WriteLine(string.Format("{0,9:F4}", avg_pearson_sq / (double)firstSet.Count).Replace(',', '.'));
                        avgStatsFile.WriteLine();
                        avgStatsFile.WriteLine("Statistics for individual atom types: ");
                        avgStatsFile.WriteLine("Atom     D_MAX    D_AVG    RMSD     PEARSON  PEARSON^2");

                        //calculate values for each atom type
                        foreach (var atomType in xValues.Keys.ToList())
                        {
                            //skip values for last molecule
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

                            avgStatsFile.Write("{0,-6}", atomType);
                            avgStatsFile.Write(string.Format("{0,9:F4}", d_max).Replace(',', '.'));
                            avgStatsFile.Write(string.Format("{0,9:F4}", d_avg).Replace(',', '.'));
                            avgStatsFile.Write(string.Format("{0,9:F4}", rmsd).Replace(',', '.'));
                            avgStatsFile.Write(string.Format("{0,9:F4}", pearson).Replace(',', '.'));
                            avgStatsFile.Write(string.Format("{0,9:F4}", pearson * pearson).Replace(',', '.'));
                            avgStatsFile.WriteLine();
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

        private List<Molecule> LoadMoleculesFromOutputFile(string filePath)
        {
            List<Molecule> set = new List<Molecule>();

            using (StreamReader reader = File.OpenText(filePath))
            {
                string line;
                int lineNum = 1;
                Molecule molecule = new Molecule();
                while ((line = reader.ReadLine()) != null)
                {
                    //check whether file is the standard output file or OpenBabel .pqr file
                    if (line.Contains("MODEL        1"))
                    {
                        //bond types are not supported for OpenBabel files
                        paramBonds = false;
                        return LoadMoleculesFromOpenBabelFile(reader);
                    }
                    if (lineNum < 1)
                        lineNum++;
                    else if (lineNum == 1)
                    {
                        molecule = new Molecule();
                        molecule.Code = line.Trim();
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
                        string[] items = line.Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                        atom.ID = int.Parse(items[0]);
                        atom.Symbol = items[1].Trim();
                        atom.Charge = double.Parse(items[2], CultureInfo.InvariantCulture);
                        //check whether file has specified bond types in the fourth column and update paramBonds flag accordingly
                        if (items.Length == 4)
                            atom.HighestBondType = int.Parse(items[3]);
                        else
                            paramBonds = false;
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

            return set;
        }

        private List<Molecule> LoadMoleculesFromOpenBabelFile(StreamReader reader)
        {
            List<Molecule> set = new List<Molecule>();

            string line;
            string code = "PARSE";
            Molecule molecule = new Molecule();
            while ((line = reader.ReadLine()) != null)
            {
                if (code.Equals("PARSE"))
                {
                    if (line.Contains("COMPND"))
                    {
                        molecule = new Molecule();
                        string[] items = line.Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                        molecule.Code = items[1];
                    }
                    if (line.Contains("AUTHOR"))
                        code = "ATOMS";
                    if (line.Contains("ENDMDL"))
                    {
                        molecule.NumOfAtoms = molecule.Atoms.Count;
                        set.Add(molecule);
                    }
                    continue;
                }
                if (code.Equals("ATOMS"))
                {
                    if (line.Contains("MASTER"))
                    {
                        code = "PARSE";
                        continue;
                    }
                    if (line.Contains("CONECT"))
                    {
                        string[] items = line.Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                        var atom = molecule.Atoms.Find(a => a.ID == int.Parse(items[1]));
                        //skip in case atom has no bonds specified
                        if (!(items.Length > 2))
                            continue;
                        for (int i = 2; i < items.Length; i++)
                            atom.Bonds.Add(int.Parse(items[i]), -1);//default value -1 because OpenBabel .pqr files do not specify bond types
                    }
                    else
                    {
                        string[] items = line.Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                        Atom atom = new Atom();
                        atom.ID = int.Parse(items[1]);
                        //check whether chainID is included and longer format is used
                        if (items.Length > 11)
                        {
                            atom.Charge = double.Parse(items[9], CultureInfo.InvariantCulture);
                            atom.Symbol = items[11].Trim();
                        }
                        else
                        {
                            atom.Charge = double.Parse(items[8], CultureInfo.InvariantCulture);
                            atom.Symbol = items[10].Trim();
                        }
                        molecule.Atoms.Add(atom);
                    }
                }
            }

            return set;
        }

        private void GNUPlot(List<Molecule> firstSet, List<Molecule> secondSet, string outputFilePath, List<string> symbols)
        {
            Process plotProcess = new Process();
            plotProcess.StartInfo.FileName = Path.GetFullPath(Path.Combine(Environment.CurrentDirectory, GNUPlotExePath));
            plotProcess.StartInfo.RedirectStandardInput = true;
            plotProcess.StartInfo.UseShellExecute = false;
            plotProcess.Start();

            NumberFormatInfo nfi = new NumberFormatInfo();
            nfi.NumberDecimalSeparator = ".";

            string inputFilePath = Path.Combine(Path.GetDirectoryName(outputFilePath), $"{Path.GetFileNameWithoutExtension(outputFilePath)}GNUPlotInputFile.txt");
            string graphFilePath = Path.Combine(Path.GetDirectoryName(outputFilePath), $"{Path.GetFileNameWithoutExtension(outputFilePath)}Graph.pdf");

            using (StreamWriter inputFile = new StreamWriter(inputFilePath))
                for (int i = 0; i < firstSet.Count; i++)
                    for (int j = 0; j < firstSet.ElementAt(i).NumOfAtoms; j++)
                    {
                        string a = firstSet.ElementAt(i).Atoms.ElementAt(j).Charge.ToString(nfi);
                        string b = secondSet.ElementAt(i).Atoms.ElementAt(j).Charge.ToString(nfi);
                        var atom = firstSet.ElementAt(i).Atoms.ElementAt(j);
                        string s = atom.Symbol;
                        if (paramBonds) s += atom.HighestBondType.ToString();
                        //third column is for point color
                        inputFile.WriteLine($"{a} {b} {symbols.IndexOf(s)}");
                    }

            using (StreamWriter gnuplotFile = plotProcess.StandardInput)
            {
                gnuplotFile.WriteLine("set terminal pdf size 3.75,2.25");
                gnuplotFile.WriteLine("set title \"Partial atomic charge correlation graph\"");
                gnuplotFile.WriteLine("set xlabel \"First set of charges\"");
                gnuplotFile.WriteLine("set ylabel \"Second set of charges\"");
                gnuplotFile.WriteLine("set xrange [-1.5:1.5]");
                gnuplotFile.WriteLine("set yrange [-1.5:1.5]");
                gnuplotFile.WriteLine($"set output '{graphFilePath}'");
                gnuplotFile.WriteLine($"set palette model RGB");
                gnuplotFile.Write("set palette defined(");
                //assign colours to numbers in input file
                int i = 0;
                string str = "";
                foreach (var s in symbols)
                {
                    str += $"{i} \"{GetAtomColorName(i)}\", ";
                    i++;
                }
                gnuplotFile.WriteLine(str.Substring(0, str.Length - 2) + ")");
                //set labels for each atom type with corresponding colours
                double pos = 1.0;
                i = 0;
                foreach (var s in symbols)
                {
                    pos -= 0.05;
                    gnuplotFile.WriteLine($"set label \"{s}\" at graph 0.05,{pos.ToString(nfi)} front tc \"{GetAtomColorName(i)}\" font \"Verdana,8\"");
                    i++;
                }
                gnuplotFile.WriteLine("set style fill transparent solid 0.75 noborder");
                gnuplotFile.WriteLine("set style circle radius 0.015");
                gnuplotFile.WriteLine("set arrow from -1.5,-1.5 to 1.5,1.5 nohead lc \"grey\"");
                gnuplotFile.WriteLine("unset colorbox");
                gnuplotFile.WriteLine("set size square");
                gnuplotFile.WriteLine($"plot '{inputFilePath}' u 1:2:3 w circles palette notitle");
            }
        }

        //return color name corresponding to input number
        private string GetAtomColorName(int n)
        {
            switch (n)
            {
                case 0: return "#000000";
                case 1: return "#ff0000";
                case 2: return "#008000";
                case 3: return "#0000ff";
                case 4: return "#800080";
                case 5: return "#ffa500";
                case 6: return "#ffff00";
                case 7: return "#adff2f";
                case 8: return "#7fffd4";
                case 9: return "#00ffff";
                case 10: return "#ff00ff";
                case 11: return "#00ffff";
                case 12: return "#008080";
                case 13: return "#a52a2a";
                case 14: return "#b8860b";
                case 15: return "#87cefa";
                case 16: return "#da70d6";
                case 17: return "#008b8b";
                default: throw new UnsupportedAtomException("Too many atom types to display labels for.");
            }
        }
    }
}
