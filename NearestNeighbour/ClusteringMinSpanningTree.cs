using Basemap.Algorithms.Contracts.Graph;
using Basemap.Algorithms.Contracts.Location;
using Basemap.Algorithms.Engine.Clustering;
using Basemap.Algorithms.Engine.Graph;
using Basemap.Algorithms.Engine.Location;
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
    public class ClusteringMinSpanningTree
    {
        //static void Main(string[] args)
        //{
        //    using (var reader = new StreamReader(@"..\..\..\points.csv"))
        //    using (var csv = new CsvReader(reader, CultureInfo.InvariantCulture))
        //    {
        //        var records = csv.GetRecords<Point>();
                
        //        Console.WriteLine("Carmen do the magic - safe");
                
        //        var graph = BuildGraph(records);

        //        //K is the number of clusters you want in output
        //        var HMSTC = new HMinSpanningTreeClustering<double>(graph, k: 20, new HashSet<IGraphNode>());

        //        var clusters= ComputeClustersConnectedComponents(HMSTC);


        //        Console.ReadKey();
        //    }

        //}


        public static IGraph<IGraphNode> BuildGraph(IEnumerable<Point> origins)
        {
            var directedAdjList = new Dictionary<IGraphNode, HashSet<IGraphEdge<IGraphNode>>>();

            var idNodesMapping = origins.ToDictionary(x => x.Id, x => new GraphNodeWithCoordinates<double>(id: x.Id, demand: 0, position: new RelativeCoordinates(x.X, x.Y)));


            foreach (var source in idNodesMapping)
            {
                directedAdjList.Add(source.Value, new HashSet<IGraphEdge<IGraphNode>>());


                foreach (var destination in idNodesMapping.Where(x => x.Key != source.Key))
                {
                    directedAdjList[source.Value].Add(new GraphEdge<IGraphNode>(source.Value, idNodesMapping[destination.Key], ComputeEuclideanDistance(source.Value, idNodesMapping[destination.Key])));
                }

            }

            return new Graph<IGraphNode>(directedAdjList);
        }


        public static float ComputeEuclideanDistance(GraphNodeWithCoordinates<double> a, GraphNodeWithCoordinates<double> b)
        {
            return (float)Math.Sqrt(Math.Pow((a.Position.Coord1 - b.Position.Coord1), 2) + Math.Pow((a.Position.Coord2 - b.Position.Coord2), 2));
        }

        public static IEnumerable<IEnumerable<IGraphNode>> ComputeClustersConnectedComponents(HMinSpanningTreeClustering<double> HMSTC)
        {


            var mst = HMSTC.Graph.MinSpanningTree();
            var representatives = HMSTC.Graph.Nodes;
            var consideringEdges = mst.SortedUniqueEdges;//The sorting doesn't matter, the uniqueness matters
            int num_clusters = 1;
            while (num_clusters != HMSTC.K)
            {
                (var avg, var stdDev) = ComputeAvgAndStdDev(consideringEdges);
                var edgesToRemove = consideringEdges.Where(x => x.Weight > avg + stdDev).ToHashSet();
                num_clusters += edgesToRemove.Count;
                mst = mst.RemoveEdgesFrom(edgesToRemove);

                if (num_clusters < HMSTC.K)
                {
                    edgesToRemove = mst.SortedUniqueEdges.TakeLast(HMSTC.K - num_clusters).ToHashSet();  //Remove the heaviest (K-num clusters)
                    mst = mst.RemoveEdgesFrom(edgesToRemove);
                    break;
                }
                else if (num_clusters > HMSTC.K)
                {
                    representatives = ComputeCentroids<double>(mst.ConnectedComponents()).ToHashSet();  //Build a new Graph made only by representatives nodes (centroids of prev clustering)
                    consideringEdges = HMSTC.Graph.ComputeEdgesToConsider(representatives); //Compute the edges to add t the graph to connect the centroids together
                    mst = mst.AddEdges(consideringEdges).MinSpanningTree(); //Add the edges and rebuild the min spanning tree again
                    consideringEdges = mst.ComputeEdgesToConsider(representatives); //The new edges to consider are the  ones linking the centroids but inside the mst
                    num_clusters = 1;
                }
            }

            //var con = mst.ConnectedComponents();
            //foreach(var x in con)
            //{
            //    foreach (var n in x)
            //    {
            //        Console.Write(n.Id+",");
            //    }
            //    Console.WriteLine();
            //}

            return mst.ConnectedComponents();
        }

        public static (double avg, double stdDev) ComputeAvgAndStdDev(HashSet<IGraphEdge<IGraphNode>> edges)
        {
            var avg = edges.Sum(x => x.Weight) / edges.Count;
            var stdDev = Math.Sqrt(edges.Select(x => Math.Pow(x.Weight - avg, 2)).Sum() / edges.Count);
            return (avg, stdDev);
        }

        public static IEnumerable<IGraphNode> ComputeCentroids<Type> (IEnumerable<IEnumerable<IGraphNode>> clustersMembers)
        {
            var centroidList = new List<IGraphNode>();
            foreach (var members in clustersMembers)
            {
                var coords = members.Select(x => ((GraphNodeWithCoordinates<Type>)x).Position).ToList();
                var median = coords.First().ComputeMedian(coords); //TO do modify this

                var centroid = FindClosestPoint(median, members);
                centroidList.Add(centroid);

            }

            return centroidList;

        }

        public static IGraphNode FindClosestPoint<Type>(ICoordinates<Type> median, IEnumerable<IGraphNode> members)
        {

            var medianNodeDistances = members.ToDictionary(x => x, x => ((GraphNodeWithCoordinates<Type>)x).Position.Distance(median)).OrderBy(x => x.Value);
            return medianNodeDistances.First().Key;

        }
    }
}
