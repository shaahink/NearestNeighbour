using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using CsvHelper;
using Google.OrTools.LinearSolver;

namespace ClusteringIP
{
    public class Point
    {
        public int Id { get; set; }
        public double X { get; set; }
        public double Y { get; set; }

    }
    class KMeans
    {
        static void Main(string[] args)
        {
            using (var reader = new StreamReader(@"..\..\..\points.csv"))
            using (var csv = new CsvReader(reader, CultureInfo.InvariantCulture))
            {
                var records = csv.GetRecords<Point>().ToArray();
                var points = new double[records.Count(), 2];
                for(int i=0;i < records.Count();i++)
                {
                    points[i, 0] = records[i].X;
                    points[i, 1] = records[i].Y;

                }

                //var points = new double[,] { { 2.3, 4.4 }, { 5.5, 10 }, { 0.7, 9 }, { 1, 1 } };
                //Number of clusters
                int K = 20;
                var assigement=Run(points, K);
                HashSet<int>[] clusters = new HashSet<int>[K];
                for(int i = 0; i < assigement.GetLength(0); i++)
                {
                    if (clusters[assigement[i]] == null)
                        clusters[assigement[i]] = new HashSet<int>();
                    clusters[assigement[i]].Add(records[i].Id);
                }
            }


        }

        public static Variable[,] Assignement(double[,] Distances, int N, int K)
        {

            //double[,] Distances = new double[N, K] { { 2.1, 1.2 }, { 3, 4 }, { 5, 6 }, { 7, 8 } };
            // Create the linear solver with the GLOP backend.
            Solver solver = Solver.CreateSolver("SimpleLpProgram", "GLOP_LINEAR_PROGRAMMING");

            // Create the variables x and y.
            Variable[,] X = solver.MakeBoolVarMatrix(N, K, "X");


            Console.WriteLine("Number of variables = " + solver.NumVariables());

            //Create a linear constraint each node should belong to only one cluster (Sum j=1..K xij=1 foreach i belong to N)
            for (int i = 0; i < N; i++)
            {
                Constraint nodeInCluster = solver.MakeConstraint(1, 1, "NodeOnlyInOneCluster");
                for (int j = 0; j < K; j++)
                    nodeInCluster.SetCoefficient(X[i, j], 1);
            }
            ////Constraint for balancing
            //for (int j = 0; j < K; j++)
            //{
            //    Constraint balanceCluster = solver.MakeConstraint(0, Math.Ceiling((double)(N / K))+1000);
            //    for (int i = 0; i < N; i++)
            //        balanceCluster.SetCoefficient(X[i, j], 1);
            //}

            Console.WriteLine("Number of constraints = " + solver.NumConstraints());
            //Create objective function, min sum i belong to N, sum j belong K of  Dij*Xij
            Objective obj = solver.Objective();
            for (int i = 0; i < N; i++)
                for (int j = 0; j < K; j++)
                {
                    obj.SetCoefficient(X[i, j], Distances[i, j]);
                }
            obj.Minimization();


            solver.Solve();

            Console.WriteLine("Solution:");
            Console.WriteLine("Objective value = " + solver.Objective().Value());
            //for (int i = 0; i < N; i++)
            //{
            //    for (int j = 0; j < K; j++)
            //    {
            //        Console.Write(X[i, j].SolutionValue() + ", ");
            //    }
            //    Console.Write(System.Environment.NewLine);
            //}

            return X;
        }


        public static List<int> ChooseCentroids(int k, int numberOfPoints)
        {
            var random = new Random();
            var centIndexes = new HashSet<int>();

            while (centIndexes.Count < k)
            {
                var randomIndex = random.Next(0, numberOfPoints);

                centIndexes.Add(randomIndex);
                //Console.WriteLine(randomIndex);
            }

            return centIndexes.ToList();
        }

        public static (int[], double[], double[,]) ComputeNodeCentroidDistances(double[,] centroidsCoord, double[,] pointsCoord)
        {

            //double[,] Distances = new double[pointsCoord.GetLength(0), centroidsCoord.GetLength(0)];
            int[] numberOfMembersPerCluster = new int[centroidsCoord.GetLength(0)];
            double[] minDistances = new double[pointsCoord.GetLength(0)];
            int[] indexes = new int[pointsCoord.GetLength(0)];
            double[,] newCentroidsCoord = new double[centroidsCoord.GetLength(0), centroidsCoord.GetLength(1)];


            for (int i = 0; i < pointsCoord.GetLength(0); i++)
            {
                double minDist = double.MaxValue;
                int minDistIndex = 0;
                for (int k = 0; k < centroidsCoord.GetLength(0); k++)
                {
                    var dist= Math.Sqrt(Math.Pow((pointsCoord[i, 0] - centroidsCoord[k, 0]), 2) + Math.Pow((pointsCoord[i, 1] - centroidsCoord[k, 1]), 2));
                    if (dist< minDist)
                    {
                        minDist = dist;
                        minDistIndex = k;
                    }
                }
                indexes[i] = minDistIndex;
                numberOfMembersPerCluster[minDistIndex] += 1;
                minDistances[i] = minDist;
                newCentroidsCoord[minDistIndex, 0] += pointsCoord[i, 0];
                newCentroidsCoord[minDistIndex, 1] += pointsCoord[i, 1];
            }
            for (int k = 0; k < centroidsCoord.GetLength(0); k++)
            {
                if (numberOfMembersPerCluster[k] > 0)
                {
                    newCentroidsCoord[k, 0] /= numberOfMembersPerCluster[k];
                    newCentroidsCoord[k, 1] /= numberOfMembersPerCluster[k];
                }
                else
                {
                    newCentroidsCoord[k, 0] = centroidsCoord[k, 0];
                    newCentroidsCoord[k, 1] = centroidsCoord[k, 1];
                }

            }

            return (indexes, minDistances, newCentroidsCoord);
        }

        public static double[,] dist2Clusters(double[,] centroidsCoord, double[,] pointsCoord)
        {
            double[,] Distances = new double[pointsCoord.GetLength(0), centroidsCoord.GetLength(0)];
            for (int i = 0; i < pointsCoord.GetLength(0); i++)
            {
                for (int k = 0; k < centroidsCoord.GetLength(0); k++)
                {
                    Distances[i, k] = Math.Sqrt(Math.Pow((pointsCoord[i, 0] - centroidsCoord[k, 0]), 2) + Math.Pow((pointsCoord[i, 1] - centroidsCoord[k, 1]), 2));
                }
            }
            return Distances;
        }

        public static (double[], int[]) GetMinDistancesAndAssignement(double[,] Distances, Variable[,] assignementVectorBinary)
        {
            double[] minDistances = new double[assignementVectorBinary.GetLength(0)];
            int[] assignementVector = new int[assignementVectorBinary.GetLength(0)];
            for (int i = 0; i < assignementVectorBinary.GetLength(0); i++)
            {
                for (int k = 0; k < assignementVectorBinary.GetLength(1); k++)
                {
                    if (assignementVectorBinary[i, k].SolutionValue() == 1)
                    {
                        minDistances[i] = Distances[i, k];
                        assignementVector[i] = k;
                        break;
                    }
                }
             }
            return (minDistances, assignementVector);
        }

        public static double[,] ComputeCentroidsCoord(double[,] pointsCoord, int[] assignementVector, double[,] oldCentroidsCoord)
        {
            double[,] centroidsCoord = new double[oldCentroidsCoord.GetLength(0), oldCentroidsCoord.GetLength(1)];
            int[] numberOfMembersPerCluster = new int[centroidsCoord.GetLength(0)];
            for (int i = 0; i < assignementVector.GetLength(0); i++)
            {
                centroidsCoord[assignementVector[i], 0] += pointsCoord[i, 0];
                centroidsCoord[assignementVector[i], 1] += pointsCoord[i, 1];
                numberOfMembersPerCluster[assignementVector[i]] += 1;
            }
            for (int k = 0; k < oldCentroidsCoord.GetLength(0); k++)
            {
                if (numberOfMembersPerCluster[k] > 0)
                {
                    centroidsCoord[k, 0] /= numberOfMembersPerCluster[k];
                    centroidsCoord[k, 1] /= numberOfMembersPerCluster[k];
                }
                else
                {
                    centroidsCoord[k, 0] = oldCentroidsCoord[k,0];
                    centroidsCoord[k, 1] = oldCentroidsCoord[k,1];
                }

            }
            return centroidsCoord;
        }

        public static int[] Run(double[,] pointsCoord,int K)
        {
            int N = pointsCoord.GetLength(0);
            double[,] centroidsCoord = new double[K, 2];
            double[,] bestCentroids;
            double[,] Distances = new double[N, K];
            int[] assignementVector = new int[N];
            double[] minDistances = new double[N];
            var bestDistances = Double.MaxValue;
            int[] bestAssignement= new int[N];
            int multiStartIterations = 10;
            for(int it = 0; it < multiStartIterations; it++)
            {
                var randomCentroidIndexes = ChooseCentroids(K, N);
                for (int k = 0; k < K; k++)
                {
                    centroidsCoord[k, 0] = pointsCoord[randomCentroidIndexes[k], 0];
                    centroidsCoord[k, 1] = pointsCoord[randomCentroidIndexes[k], 1];
                }

                (assignementVector, minDistances, centroidsCoord) = ComputeNodeCentroidDistances(centroidsCoord, pointsCoord);
                var prevDistances = minDistances.Average();
                var improvement = prevDistances;
                int iterations = 0;
                while (improvement > 0.01 && iterations < 10)
                {
                    (minDistances, assignementVector) = ClusterDistIdxIP(centroidsCoord, pointsCoord);
                    centroidsCoord = ComputeCentroidsCoord(pointsCoord, assignementVector, centroidsCoord);
                    var newDistances = minDistances.Average();
                    improvement = prevDistances - newDistances;
                    prevDistances = newDistances;
                    iterations++;
                }

                if (prevDistances < bestDistances)
                {
                    bestDistances = prevDistances;
                    bestAssignement = assignementVector;
                    bestCentroids = centroidsCoord;
                }
            }
            Console.WriteLine("Best Assignement: {0}",bestDistances);
            foreach(var assign in bestAssignement)
            {
                Console.WriteLine(assign);
            }
            return bestAssignement;

        }

        public static (double[], int[]) ClusterDistIdxIP(double[,] centroidsCoord, double[,] pointsCoord)
        {
            var Distances = dist2Clusters(centroidsCoord, pointsCoord);
            var assignementVectorBinary = Assignement(Distances, pointsCoord.GetLength(0), centroidsCoord.GetLength(0));

            return GetMinDistancesAndAssignement(Distances, assignementVectorBinary);
        }
    }



}
