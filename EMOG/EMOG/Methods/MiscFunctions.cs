using EMOG.Structures;
using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace EMOG.Methods
{
    static class MiscFunctions
    {
        public static void SaveResults(StreamWriter file, Molecule molecule, Vector<double> results)
        {
            file.WriteLine($"{molecule.Code}");
            file.WriteLine(molecule.NumOfAtoms);
            int count = 0;
            foreach (double charge in results)
            {
                if (count != molecule.NumOfAtoms)
                {
                    file.Write("{0,5}", count + 1);
                    file.Write("{0,3}", molecule.Atoms[count].Symbol);
                    file.Write(string.Format("{0,12:F6}", charge).Replace(',', '.'));
                    file.WriteLine("{0,2}", molecule.Atoms[count].HighestBondType);
                    count++;
                }
            }
            file.WriteLine("$$$$");
        }

        public static void SaveIncorrectResults(StreamWriter file, Molecule molecule, Exception ex)
        {
            file.WriteLine($"{molecule.Code}");
            file.WriteLine("X");
            file.WriteLine($"Cannot compute charges for this molecule: {ex.Message}");
            file.WriteLine("$$$$");
        }

        public static void PrintHelp()
        {
            Console.WriteLine("Arguments must be entered in the following format:");
            Console.WriteLine("\"<method> <molecules file> <output file> <parameters file>\" or \"stats <first charges file> <second charges file> <output file> <GNUPlot EXE file> <EEM bonds flag>\"");
            Console.WriteLine();
            Console.WriteLine("<method> - Method name. Only \"eem\", \"mgc\" or \"ogc\" supported.");
            Console.WriteLine("<molecules file> - Path and file name of the set of molecules.");
            Console.WriteLine("<charges file> - Path and file name of the set of partial charges.");
            Console.WriteLine("<output file> - Desired path and file name for output file.");
            Console.WriteLine("<GNUPlot .exe file> - Path to GNUPlot executable file.");
            Console.WriteLine("<parameters file> - Path and file name of the parameters file. Required only for EEM.");
            Console.WriteLine("<EEM bonds flag> - \"y\" or \"n\" to specify whether statisticss should be generated for atoms according to their bond types or not. Not supported for OpenBabel output files.");
        }

        public static bool CanParseArguments(string[] items)
        {
            if (items.Length == 0)
            {
                Console.WriteLine("No arguments.");
                return false;
            }
            else if (!(items[0].Equals("eem") || items[0].Equals("mgc") || items[0].Equals("ogc") || items[0].Equals("stats")))
            {
                Console.WriteLine("Incorrect method name. Only \"eem\", \"mgc\", \"ogc\" or \"stats\" supported.");
                return false;
            }
            else if (items[0].Equals("stats") && items.Length != 6)
            {
                Console.WriteLine("Incorrect number of arguments for statistics generation.");
                return false;
            }
            else if (items[0].Equals("stats") && !(items[5].Equals("n") || items[5].Equals("y")))
            {
                Console.WriteLine("Please specify whether you want to generate statistics for atom types according to their bond types with either \"y\" or \"n\"");
                return false;
            }
            else if (((items[0].Equals("mgc") || items[0].Equals("ogc")) && items.Length != 3) || (items[0].Equals("eem") && items.Length != 4))
            {
                Console.WriteLine($"Incorrect number of arguments for {items[0]} method.");
                return false;
            }
            return true;
        }

        public static string FileVersion(string fileName)
        {
            try
            {
                using (StreamReader reader = File.OpenText(fileName))
                {
                    string line;
                    int lineNum = 1;
                    while ((line = reader.ReadLine()) != null)
                    {
                        if (lineNum < 4) lineNum++;
                        else
                        {
                            if (line.Contains("V2000")) return "V2000";
                            if (line.Contains("V3000")) return "V3000";
                        }
                    }
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine("Could not read .sdf file. Exception: " + ex.Message);                
            }
            return "error";
        }

        //return the highest bond type of an atom's bonds
        public static int GetHighestBondType(Molecule m, Atom a)
        {
            int bond = 0;

            foreach (var b in a.Bonds)
                if (b.Value > bond)
                    bond = b.Value;

            foreach (var atom in m.Atoms)
                foreach (var b in atom.Bonds)
                    if (b.Key == a.ID && b.Value > bond)
                        bond = b.Value;

            return bond;
        }

        public static List<Molecule> LoadMoleculesV2000(StreamReader reader)
        {
            List<Molecule> molecules = new List<Molecule>();

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
                    var atom = molecule.Atoms.Find(a => a.ID == int.Parse(line.Substring(0, 3)));
                    atom.Bonds.Add(
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

            return molecules;
        }

        public static List<Molecule> LoadMoleculesV3000(StreamReader reader)
        {
            List<Molecule> molecules = new List<Molecule>();
            
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
                    if (line.Contains("BEGIN BOND"))
                    {
                        code = "BONDS";
                        continue;
                    }
                    if (line.Contains("$$$$"))
                    {
                        code = "START";
                        continue;
                    }
                }
                if (code.Equals("START"))
                {
                    molecule.Code = line.Trim();
                    code = "PARSE";
                }
                if (code.Equals("COUNTS"))
                {
                    string[] items = line.Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
                    molecule.NumOfAtoms = int.Parse(items[3]);
                    molecule.NumOfBonds = int.Parse(items[4]);
                    code = "PARSE";
                }
                if (code.Equals("ATOMS"))
                {
                    if (line.Contains("END ATOM"))
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
                    if (line.Contains("END BOND"))
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

            return molecules;
        }

        public static List<Molecule> LoadMolecules(string fileName)
        {
            List<Molecule> molecules = new List<Molecule>();
            
            using (StreamReader reader = File.OpenText(fileName))
            {
                if (FileVersion(fileName).Equals("V2000"))
                    molecules = LoadMoleculesV2000(reader);
                else if (FileVersion(fileName).Equals("V3000"))
                    molecules = LoadMoleculesV3000(reader);
                else
                    throw new Exception("Unsupported file format - only SDF file versions V2000 and V3000 supported.");
            }

            return molecules;
        }
    }
}
