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

        public ElementParameters(string elemName, double a, double b, int bondType = 0)
        {
            ElementName = elemName;
            A = a;
            B = b;
            BondType = bondType;
        }
    }
}
