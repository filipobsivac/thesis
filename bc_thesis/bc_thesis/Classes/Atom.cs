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
        public string Symbol { get; set; }
        public double X { get; set; }
        public double Y { get; set; }
        public double Z { get; set; }
        public double Charge { get; set; }

        // <bonded atom ID, bond type>       
        public Dictionary<int, int> Bonds { get; set; }

        public Atom()
        {            
            Bonds = new Dictionary<int, int>();
        }
    }
}
