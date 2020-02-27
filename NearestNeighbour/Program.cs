using CsvHelper;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace NearestNeighbour
{
    public class Point
    {
        public int Id { get; set; }
        public double X { get; set; }
        public double Y { get; set; }

    }
    public class Program
    {
        static void Main(string[] args)
        {
            using (var reader = new StreamReader("Origin.csv"))
            using (var csv = new CsvReader(reader, CultureInfo.InvariantCulture))
            {
                var records = csv.GetRecords<Point>().ToList();

                Console.WriteLine("Carmen do the magic - safe");
                Console.ReadKey();
            }

        }
    }
}
