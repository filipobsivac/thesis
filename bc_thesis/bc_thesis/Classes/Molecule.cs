﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace thesis.Classes
{
    public class Molecule
    {
        public string Code { get; set; }
        public int NumOfAtoms { get; set; }
        public int NumOfOrbitals { get; set; }
        public int NumOfBonds { get; set; }
        public List<Atom> Atoms { get; set; }

        public Molecule()
        {
            Atoms = new List<Atom>();
        }
    }
}
