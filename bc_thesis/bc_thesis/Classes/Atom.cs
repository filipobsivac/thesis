using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace thesis.Classes
{
    public class Atom
    {
        //atom number in molecule
        public int ID { get; set; }
        public int OGCID { get; set; }
        public string Symbol { get; set; }
        public double X { get; set; }
        public double Y { get; set; }
        public double Z { get; set; }
        public double Charge { get; set; }
        public double OrbitalCharge { get; set; }
        public int HighestBondType { get; set; }
        // <type: s/p/n, charge>
        public Dictionary<string, double> OrbitalENs { get; set; }
        // <type: s/p/n, hardness>
        public Dictionary<string, double> OrbitalHardnesses { get; set; }
        // <bonded atom ID, bond type>       
        public Dictionary<int, int> Bonds { get; set; }
        public string OrbitalID { get; set; }
        // <bonded orbital ID, orbital charge type: s/p/n/x}>
        public Dictionary<string, string> OrbitalBonds { get; set; }

        public Atom()
        {            
            Bonds = new Dictionary<int, int>();
            OrbitalENs = new Dictionary<string, double>();
            OrbitalBonds = new Dictionary<string, string>();
            OrbitalHardnesses = new Dictionary<string, double>();
        }
    }
}
