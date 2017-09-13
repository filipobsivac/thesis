using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace thesis.Classes
{
    public class ParameterSet
    {
        public string AtomType { get; set; }
        public double Kappa { get; set; }
        public List<ElementParameters> ElementParams;        

        public ParameterSet()
        {
            ElementParams = new List<ElementParameters>();
        }
    }
}
