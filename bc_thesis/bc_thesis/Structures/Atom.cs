using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace bc_thesis.Structures
{
    public class Atom
    {
        //atom ID in molecule
        public int ID { get; set; }
        //atom orbital ID
        public int OrbitalID { get; set; }
        public string Symbol { get; set; }
        public double X { get; set; }
        public double Y { get; set; }
        public double Z { get; set; }
        public double Charge { get; set; }
        public double OrbitalCharge { get; set; }
        public int HighestBondType { get; set; }

        // <bonded atom ID, bond type>       
        public Dictionary<int, int> Bonds { get; set; }

        // <bonded OrbitalID, orbital charge type: s/p/n/x}>
        //s - sigma, p - pi, n - lone electron pair, x - bond to orbital from the same atom
        public Dictionary<int, string> OrbitalBonds { get; set; }

        // <type: s/p/n, charge>
        public Dictionary<string, double> OrbitalENs { get; set; }

        // <type: s/p/n, hardness>
        public Dictionary<string, double> OrbitalHardnesses { get; set; }        

        public Atom()
        {            
            Bonds = new Dictionary<int, int>();
            OrbitalENs = new Dictionary<string, double>();
            OrbitalBonds = new Dictionary<int, string>();
            OrbitalHardnesses = new Dictionary<string, double>();
        }
    }
}
