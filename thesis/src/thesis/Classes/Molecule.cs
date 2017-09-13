using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace thesis.Classes
{
    public class Molecule
    {
        public int ID { get; set; }
        public List<Atom> Atoms { get; set; }

        Molecule()
        {
            Atoms = new List<Atom>();
        }
    }
}
