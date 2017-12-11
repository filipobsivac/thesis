using EMOG.Structures;
using MathNet.Numerics.LinearAlgebra;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using EMOG.Methods;
using static EMOG.Methods.MiscFunctions;

namespace EMOG
{
    class Program
    {
        static void Main(string[] args)
        {
            if (!CanParseArguments(args))
            {
                PrintHelp();
                return;
            }                
            switch (args[0])
            {
                case "stats": new Statistics().GenerateStatistics(args[1], args[2], args[3], args[4], args[5]); break;
                case "eem": new EEM().SolveEEM(args[3], args[2], args[1]); break;
                case "mgc": new MGC().SolveMGC(args[2], args[1]); break;
                case "ogc": new OGC().SolveOGC(args[2], args[1]); break;
                default: PrintHelp(); break;
            }         
        }
    }
}
