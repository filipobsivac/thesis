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
    class OGC
    {
        private List<Molecule> molecules;

        private const string CovalentRadiusTablePath = @"..\..\Data\CovalentRadius.csv";
        private const string OrbitalHardnessAndENTablePath = @"..\..\Data\OrbitalHardnessAndEN.csv";

        private Vector<double> BuildOGCVector(Molecule m)
        {
            double[] arr = new double[m.Atoms.Count];

            foreach (var a in m.Atoms)
            {
                //don't count bonds with orbitals from the same atom (marked by "x")
                var bond = a.OrbitalBonds.First(x => !x.Value.Equals("x")).Value;
                var en = a.OrbitalENs[bond];
                arr[a.OrbitalID - 1] = en;
            }

            return Vector<double>.Build.Dense(arr);
        }

        private Matrix<double> BuildConnectivityMatrixOGC(Molecule m)
        {
            double[,] arr = new double[m.Atoms.Count, m.Atoms.Count];

            foreach (var a in m.Atoms)
                foreach (var b in a.OrbitalBonds)
                    //don't count lone electron pairs (marked by "n")
                    if (!b.Value.Equals("n"))
                    {
                        var o = m.Atoms.Find(x => x.OrbitalID.Equals(b.Key));
                        arr[a.OrbitalID - 1, o.OrbitalID - 1] = 1;
                    }

            return Matrix<double>.Build.DenseOfArray(arr);
        }

        private Matrix<double> BuildDegreeMatrixOGC(Molecule m)
        {
            double[,] arr = new double[m.Atoms.Count, m.Atoms.Count];

            int x = 0;
            foreach (var a in m.Atoms)
            {
                int n = 0;
                foreach (var b in a.OrbitalBonds)
                    //don't count lone electron pairs (marked by "n")
                    if (!b.Value.Equals("n"))
                        n++;
                arr[x, x] = n;
                x++;
            }

            return Matrix<double>.Build.DenseOfArray(arr);
        }

        private Matrix<double> BuildIdentityMatrixOGC(Molecule m)
        {
            double[,] arr = new double[m.Atoms.Count, m.Atoms.Count];

            for (int i = 0; i < m.Atoms.Count; i++)
                arr[i, i] = 1;

            return Matrix<double>.Build.DenseOfArray(arr);
        }

        private double GetCovalentRadius(string symbol)
        {
            try
            {
                using (StreamReader reader = File.OpenText(CovalentRadiusTablePath))
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
                Console.WriteLine("Could not get covalent radius. Exception: " + ex.Message);
            }
            throw new UnsupportedAtomException($"Could not find covalent radius values for {symbol}.");
        }

        //set orbital EN and hardness to each atom according to its bond type
        private void SetOrbitalENsAndHardnesses(Molecule molecule)
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
                    default: throw new UnsupportedAtomException($"Cannot set Orbital EN and hardness for {atom.Symbol + atom.HighestBondType.ToString()}.");
                }
            }
        }

        private double GetOrbitalEN(string symbol, string state)
        {
            try
            {
                using (StreamReader reader = File.OpenText(OrbitalHardnessAndENTablePath))
                {
                    string line;
                    while ((line = reader.ReadLine()) != null)
                    {
                        string[] items = line.Split(';');
                        if (items[0].Equals(symbol) && items[1].Equals(state))
                            return double.Parse(items[2], CultureInfo.InvariantCulture);
                    }
                }
            }
            catch (Exception ex)
            {
                Console.WriteLine("Could not get orbital EN. Exception: " + ex.Message);
            }
            throw new UnsupportedAtomException($"Could not find orbital EN for {symbol} - {state}.");
        }

        private double GetOrbitalHardness(string symbol, string state)
        {
            try
            {
                using (StreamReader reader = File.OpenText(OrbitalHardnessAndENTablePath))
                {
                    string line;
                    while ((line = reader.ReadLine()) != null)
                    {
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
            throw new UnsupportedAtomException($"Could not find orbital hardness for {symbol} - {state}.");
        }

        //return the first free orbital of an atom which is not yet bonded with any other orbital (only has bonds with related orbitals marked by "x")
        private Atom GetFreeOrbital(Molecule m, int id)
        {
            return m.Atoms.Find(a => (a.ID == id) && !a.OrbitalBonds.Any(b => !b.Value.Equals("x")));
        }

        //add all bonds in which an atom appears to its own bonds        
        private void MakeBondsNonOriented(Molecule m)
        {
            foreach (var a in m.Atoms)
            {
                foreach (var b in a.Bonds)
                {
                    var bondedAtom = m.Atoms.Find(x => x.ID == b.Key);
                    if (!bondedAtom.Bonds.Any(x => x.Key == a.ID))
                        bondedAtom.Bonds.Add(a.ID, b.Value);
                }
            }
        }

        //only allow supported atom bond types
        private bool AreBondsSupported(Atom a)
        {
            if (a.Symbol.Equals("H") && a.Bonds.Sum(b => b.Value) != 1)
                return false;
            if (a.Symbol.Equals("C"))
            {
                if (a.Bonds.Sum(b => b.Value) != 4)
                    return false;
                if (a.Bonds.Count == 2)
                    if (a.Bonds.First().Value == 2 && a.Bonds.Last().Value == 2)
                        return false;
            }
            if (a.Symbol.Equals("N") && a.Bonds.Sum(b => b.Value) != 3)
                return false;
            if ((a.Symbol.Equals("O") || a.Symbol.Equals("S")) && a.Bonds.Sum(b => b.Value) != 2)
                return false;
            return true;
        }

        private Molecule BuildOGCMolecule(Molecule m)
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
                if (!AreBondsSupported(a))
                {
                    string s = "Bonds = ";
                    foreach (var b in a.Bonds)
                        s += $"{m.Atoms.Find(atom => atom.ID == b.Key).Symbol}{b.Value} ";
                    throw new UnsupportedAtomException($"Atom {a.Symbol} has an unsupported bond type. {s}");
                }

                //H always has only 1 orbital, otherwise there are 4
                int n = a.Symbol.Equals("H") ? 1 : 4;
                for (int i = 1; i <= n; i++)
                {
                    Atom atom = new Atom();
                    atom.OrbitalID = x;
                    atom.ID = a.ID;
                    atom.Symbol = a.Symbol;
                    atom.Bonds = a.Bonds;
                    atom.OrbitalENs = a.OrbitalENs;
                    atom.OrbitalHardnesses = a.OrbitalHardnesses;
                    atom.HighestBondType = a.HighestBondType;
                    ogcMolecule.Atoms.Add(atom);
                    x++;
                }
            }

            foreach (var atom in m.Atoms)
            {
                //connect all orbitals from the same atom
                var atomOrbitals = ogcMolecule.Atoms.FindAll(orb => orb.ID == atom.ID);
                foreach (var orbital in atomOrbitals)
                {
                    var relatedOrbitals = atomOrbitals.FindAll(orb => orb.OrbitalID != orbital.OrbitalID);
                    relatedOrbitals.ForEach(orb => orbital.OrbitalBonds.Add(orb.OrbitalID, "x"));
                }

                //set an orbital as a lone electron pair
                if (atom.Symbol.Equals("N"))
                {
                    var freeOrbital = GetFreeOrbital(ogcMolecule, atom.ID);
                    freeOrbital.OrbitalBonds.Add(freeOrbital.OrbitalID, "n");
                }

                //set two orbitals as lone electron pairs
                if (atom.Symbol.Equals("O") || atom.Symbol.Equals("S"))
                {
                    for (int i = 0; i < 2; i++)
                    {
                        var freeOrbital = GetFreeOrbital(ogcMolecule, atom.ID);
                        freeOrbital.OrbitalBonds.Add(freeOrbital.OrbitalID, "n");
                    }
                }

                //add orbital bonds
                foreach (var bond in atom.Bonds)
                {
                    //skip if orbital bonds are already set for this atom bond
                    if (GetFreeOrbital(ogcMolecule, bond.Key) == null) continue;

                    Atom sourceOrbital;
                    Atom targetOrbital;

                    if (bond.Value == 1)
                    {
                        targetOrbital = GetFreeOrbital(ogcMolecule, bond.Key);
                        sourceOrbital = GetFreeOrbital(ogcMolecule, atom.ID);
                        targetOrbital.OrbitalBonds.Add(sourceOrbital.OrbitalID, "s");
                        sourceOrbital.OrbitalBonds.Add(targetOrbital.OrbitalID, "s");
                    }
                    if (bond.Value == 2)
                    {
                        targetOrbital = GetFreeOrbital(ogcMolecule, bond.Key);
                        sourceOrbital = GetFreeOrbital(ogcMolecule, atom.ID);
                        targetOrbital.OrbitalBonds.Add(sourceOrbital.OrbitalID, "p");
                        sourceOrbital.OrbitalBonds.Add(targetOrbital.OrbitalID, "p");

                        targetOrbital = GetFreeOrbital(ogcMolecule, bond.Key);
                        sourceOrbital = GetFreeOrbital(ogcMolecule, atom.ID);
                        targetOrbital.OrbitalBonds.Add(sourceOrbital.OrbitalID, "s");
                        sourceOrbital.OrbitalBonds.Add(targetOrbital.OrbitalID, "s");
                    }
                    if (bond.Value == 3)
                    {
                        targetOrbital = GetFreeOrbital(ogcMolecule, bond.Key);
                        sourceOrbital = GetFreeOrbital(ogcMolecule, atom.ID);
                        targetOrbital.OrbitalBonds.Add(sourceOrbital.OrbitalID, "p");
                        sourceOrbital.OrbitalBonds.Add(targetOrbital.OrbitalID, "p");

                        targetOrbital = GetFreeOrbital(ogcMolecule, bond.Key);
                        sourceOrbital = GetFreeOrbital(ogcMolecule, atom.ID);
                        targetOrbital.OrbitalBonds.Add(sourceOrbital.OrbitalID, "s");
                        sourceOrbital.OrbitalBonds.Add(targetOrbital.OrbitalID, "s");

                        targetOrbital = GetFreeOrbital(ogcMolecule, bond.Key);
                        sourceOrbital = GetFreeOrbital(ogcMolecule, atom.ID);
                        targetOrbital.OrbitalBonds.Add(sourceOrbital.OrbitalID, "p");
                        sourceOrbital.OrbitalBonds.Add(targetOrbital.OrbitalID, "p");
                    }
                }
            }

            return ogcMolecule;
        }

        public void SolveOGC(string outputFilePath, string moleculesFilePath)
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
                            var m = BuildOGCMolecule(molecule);
                            var a = BuildConnectivityMatrixOGC(m);
                            var d = BuildDegreeMatrixOGC(m);
                            var i = BuildIdentityMatrixOGC(m);
                            var vector = BuildOGCVector(m);
                            Matrix<double> s = d - a + i;
                            var equalizedEN = s.Solve(vector);
                            var deltaEN = equalizedEN - vector;

                            //compute Dm
                            double dividend = 0;
                            double divisor = 0;
                            int n = 0;
                            foreach (var orbital in m.Atoms)
                            {
                                var bond = orbital.OrbitalBonds.First(x => !x.Value.Equals("x")).Value;
                                var orbHardness = orbital.OrbitalHardnesses.First(x => x.Key.Equals(bond)).Value;
                                dividend += deltaEN[n] / orbHardness;
                                divisor += (deltaEN[n] * deltaEN[n]) / (Math.Pow(orbHardness, 3) * GetCovalentRadius(orbital.Symbol + orbital.HighestBondType.ToString()));
                                n++;
                            }
                            double dm = dividend / divisor;

                            //compute orbital charges
                            n = 0;
                            foreach (var orbital in m.Atoms)
                            {
                                var bond = orbital.OrbitalBonds.First(x => !x.Value.Equals("x")).Value;
                                var orbHardness = orbital.OrbitalHardnesses.First(x => x.Key.Equals(bond)).Value;
                                orbital.OrbitalCharge = (deltaEN[n] / orbHardness) - (((deltaEN[n] * deltaEN[n]) * dm) / (Math.Pow(orbHardness, 3) * GetCovalentRadius(orbital.Symbol + orbital.HighestBondType.ToString())));
                                n++;
                            }

                            //sum orbital charges for each atom
                            double[] arr = new double[m.NumOfAtoms];
                            foreach (var atom in m.Atoms)
                            {
                                var orbitals = m.Atoms.FindAll(o => o.ID == atom.ID);
                                arr[atom.ID - 1] = orbitals.Sum(o => o.OrbitalCharge);
                            }
                            Vector<double> results = Vector<double>.Build.Dense(arr);

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
