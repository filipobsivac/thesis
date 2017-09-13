using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace thesis.Classes
{
    public class ElementParameters
    {
        public string ElementName { get; set; }

        //0 if not specified
        public int BondType { get; set; }
        public double A { get; set; }
        public double B { get; set; }
        public double Kappa { get; set; }

        public ElementParameters(int bondType = 0)
        {
            BondType = bondType;
        }
    }
}
