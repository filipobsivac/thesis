using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace bc_thesis.Classes
{
    class UnsupportedElementException : Exception
    {
        public UnsupportedElementException()
        {
        }

        public UnsupportedElementException(string message) : base(message)
        {
        }

        public UnsupportedElementException(string message, Exception inner) : base(message, inner)
        { 
        }
    }
}
