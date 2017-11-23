using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace bc_thesis.Structures
{
    class UnsupportedAtomException : Exception
    {
        public UnsupportedAtomException()
        {
        }

        public UnsupportedAtomException(string message) : base(message)
        {
        }

        public UnsupportedAtomException(string message, Exception inner) : base(message, inner)
        { 
        }
    }
}
