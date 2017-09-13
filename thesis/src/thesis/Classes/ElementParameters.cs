using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace thesis.Classes
{
    public class ElementParameters
    {
        public string Name { get; set; }

        //key = bond type
        public Dictionary<int, double> A { get; set; }
        public Dictionary<int, double> B { get; set; }

        public ElementParameters()
        {
            A = new Dictionary<int, double>();
            B = new Dictionary<int, double>();            
        }
    }
}
