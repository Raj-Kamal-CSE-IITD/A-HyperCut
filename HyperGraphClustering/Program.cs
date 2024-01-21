using System;
using System.IO;

using MathNet.Numerics.LinearAlgebra;
using System.Collections.Generic;
using System.Linq;
using System.Numerics;
using System.Diagnostics;
using System.Data;
using System.Threading.Tasks.Sources;
using System.Collections;
using System.Runtime.Serialization.Json;

namespace HyperGraphCLustering
{
    class MainClass
    {
        public static double Sqr(double v) { return v * v; }
        const double dt = 0.01;


        public static MathNet.Numerics.LinearAlgebra.Vector<double> Integrate(Func<double, MathNet.Numerics.LinearAlgebra.Vector<double>> f, double T)
        {
            var dummy = f(0.0);
            var res = CreateVector.Dense<double>(dummy.Count);

            var t = 0.0;
            while (t < T)
            {
                double nt = Math.Min(T - t, dt);
                var v = f(t);
                res += v * nt;
                t += nt;
            }
            return res;
        }


        // for Proposed_local
        public static void Proposed_local(string fn, int v_init)
        {
            var H = HyperGraph.Open(fn);

            var time = new System.Diagnostics.Stopwatch();
            time.Start();

            int n = H.n;
            int m = H.m;

            const double eps = 0.9;

            const double dt = 1.0;
            const double T = 30.0;

            var A_cand = new List<double>();
            for (int i = 0; i <= Math.Log(n * m) / Math.Log(1 + eps); i++)
            {
                A_cand.Add(Math.Pow(1 + eps, i) / (n * m));
            }

            var edge_size = new Dictionary<int, int>();
            for (int eid = 0; eid < H.m; eid++)
            {
                edge_size.Add(eid, H.edge_ID_rev[eid].Count());
            }

            double min_conductance = double.MaxValue;

            foreach (double alpha in A_cand)
            {
                var vec = CreateVector.Dense<double>(n);

                vec[v_init] = 1.0;

                vec = HyperGraph.Simulate(H, vec, v_init, dt, T, alpha);

                for (int i = 0; i < n; i++)
                {
                    vec[i] /= H.w_Degree(i);
                }

                int[] index = Enumerable.Range(0, n).ToArray<int>();
                Array.Sort<int>(index, (a, b) => vec[a].CompareTo(vec[b]));

                Array.Reverse(index);

                double vol_V = 0;
                for (int i = 0; i < n; i++) vol_V += H.w_Degree(i);

                var num_contained_nodes = new Dictionary<int, int>();
                for (int eid = 0; eid < H.m; eid++)
                {
                    num_contained_nodes.Add(eid, 0);
                }

                double cut_val = 0;
                double vol_S = 0;
                double conductance = double.MaxValue;
                int best_index = -1;

                foreach (int i in index)
                {
                    vol_S += H.w_Degree(i);
                    if (vol_S <= vol_V / 10.0)
                    {
                        foreach (var e in H.incident_edges[i])
                        {
                            if (num_contained_nodes[e] == 0)
                            {
                                cut_val += H.edge_weights[e];
                            }
                            if (num_contained_nodes[e] == edge_size[e] - 1)
                            {
                                cut_val -= H.edge_weights[e];
                            }
                            num_contained_nodes[e] += 1;
                        }
                        conductance = cut_val / Math.Min(vol_S, vol_V - vol_S);
                        if (conductance < min_conductance)
                        {
                            min_conductance = conductance;
                            best_index = i;
                        }
                    }
                    else
                    {
                        break;
                    }
                }
            }
            time.Stop();
            TimeSpan ts = time.Elapsed;

            Console.WriteLine("conductance: " + min_conductance);
            Console.WriteLine("time(s): " + time.ElapsedMilliseconds / 1000.0);
        }



        // for Proposed_local_round
        public static List<double> Proposed_local_round(HyperGraph H, string out_file, int v_init)
        {
            //var H = HyperGraph.Open(in_file);

            var time = new System.Diagnostics.Stopwatch();
            time.Start();

            int n = H.n;
            int m = H.m;
            double md = m;
            double nd = n;
            const double eps = 0.9;

            const double dt = 1.0;
            const double T = 30.0;
            var A_cand = new List<double>();
            for (int i = 0; i <= Math.Log(nd * md) / Math.Log(1 + eps); i++)
            {
                A_cand.Add(Math.Pow(1 + eps, i) / (nd * md));
            }

            var edge_size = new Dictionary<int, int>();
            for (int eid = 0; eid < m; eid++)
            {
                edge_size.Add(eid, H.edge_ID_rev[eid].Count());
            }

            double min_conductance = double.MaxValue;
            foreach (double alpha in A_cand)
            {

                var vec = CreateVector.Dense<double>(n);

                vec[v_init] = 1.0;

                vec = HyperGraph.Simulate_round(H, vec, v_init, dt, T, alpha);

                for (int i = 0; i < n; i++)
                {
                    vec[i] /= H.w_Degree(i);
                }

                int[] index = Enumerable.Range(0, n).ToArray<int>();
                Array.Sort<int>(index, (a, b) => vec[a].CompareTo(vec[b]));

                Array.Reverse(index);

                double vol_V = 0;
                for (int i = 0; i < n; i++) vol_V += H.w_Degree(i);

                var num_contained_nodes = new Dictionary<int, int>();
                for (int eid = 0; eid < H.m; eid++)
                {
                    num_contained_nodes.Add(eid, 0);
                }

                double cut_val = 0;
                double vol_S = 0;
                double conductance = double.MaxValue;
                int best_index = -1;

                foreach (int i in index)
                {
                    vol_S += H.w_Degree(i);
                    if (vol_S <= vol_V / 10.0)
                    {
                        foreach (var e in H.incident_edges[i])
                        {
                            if (num_contained_nodes[e] == 0)
                            {
                                cut_val += H.edge_weights[e];
                            }
                            if (num_contained_nodes[e] == edge_size[e] - 1)
                            {
                                cut_val -= H.edge_weights[e];
                            }
                            num_contained_nodes[e] += 1;
                        }
                        conductance = cut_val / Math.Min(vol_S, vol_V - vol_S);
                        if (conductance < min_conductance)
                        {
                            min_conductance = conductance;
                            best_index = i;
                        }
                    }
                    else
                    {
                        break;
                    }
                }
            }
            time.Stop();
            TimeSpan ts = time.Elapsed;
            File.AppendAllText(out_file, v_init + "\t" + min_conductance + "\t" + time.ElapsedMilliseconds / 1000.0 + "\n");
            List<double> answer = new()
            {
                min_conductance,
                time.ElapsedMilliseconds / 1000.0
            };
            return answer;
        }


        public static void Alpha_Average_Clustering(HyperGraph H, string out_file, int v_init)
        {
            //var H = HyperGraph.Open(in_file);
            File.AppendAllText(out_file, "Seed vertex: " + v_init + "\n");
            File.AppendAllText(out_file, "alpha\tconductance\t\ttime taken (seconds)\n");
            int n = H.num_nodes;
            int m = H.num_edges;

            var edge_size = new Dictionary<int, int>();
            for (int eid = 0; eid < H.m; eid++)
            {
                edge_size.Add(eid, H.edge_ID_rev[eid].Count);
            }

            var myA = new List<double>() { 1e-5, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95 };
            var M = Matrix<double>.Build;
            var temp_vert = M.Dense(n, 4);
            for (int i = 0; i < n; i++)
            {
                var weight = H.vertex_weights[i];
                if (weight > 0)
                {
                    temp_vert[i, 2] = 1.0 / weight;
                }
            }

            var temp_edge = M.Dense(m, 3);
            for (int i = 0; i < m; i++)
            {
                var weight = H.edge_ID_rev[i].Count;
                if (weight > 0)
                {
                    temp_edge[i, 1] = 1.0 / weight;
                }
            }

            foreach (double alpha in myA)
            {
                var time = new System.Diagnostics.Stopwatch();
                time.Start();
                var vec = HyperGraph.A_HyperCut(H, temp_vert, temp_edge, v_init, alpha);
                for (int i = 0; i < n; i++)
                {
                    vec[i] /= H.vertex_weights[i];
                }
                int[] index = Enumerable.Range(0, n).ToArray<int>();
                Array.Sort<int>(index, (a, b) => vec[a].CompareTo(vec[b]));

                Array.Reverse(index);

                double vol_V = 0;
                for (int i = 0; i < n; i++) vol_V += H.vertex_weights[i];
                var num_contained_nodes = new Dictionary<int, int>();
                for (int eid = 0; eid < H.m; eid++)
                {
                    num_contained_nodes.Add(eid, 0);
                }

                double cut_val = 0;
                double vol_S = 0;
                double conductance = double.MaxValue;
                double min_conductance = double.MaxValue;

                foreach (int i in index)
                {
                    vol_S += H.w_Degree(i);
                    if (vol_S <= vol_V / 2.0)
                    {
                        foreach (var e in H.incident_edges[i])
                        {
                            if (num_contained_nodes[e] == 0)
                            {
                                cut_val += H.edge_weights[e];
                            }
                            if (num_contained_nodes[e] == edge_size[e] - 1)
                            {
                                cut_val -= H.edge_weights[e];
                            }
                            num_contained_nodes[e] += 1;
                        }
                        conductance = cut_val / Math.Min(vol_S, vol_V - vol_S);
                        if (conductance < min_conductance)
                        {
                            min_conductance = conductance;
                        }
                    }
                    else
                    {
                        break;
                    }
                }
                time.Stop();
                TimeSpan ts = time.Elapsed;
                File.AppendAllText(out_file, alpha + "\t" + min_conductance + "\t" + time.ElapsedMilliseconds / 1000.0 + "\n");
            }
            File.AppendAllText(out_file, "\n");
        }

        public static List<double> Average_Clustering(HyperGraph H, string out_file, int v_init, double epsilon)
        {
            //var H = HyperGraph.Open(in_file);

            var time = new System.Diagnostics.Stopwatch();
            time.Start();

            int n = H.num_nodes;
            int m = H.num_edges;

            var edge_size = new Dictionary<int, int>();
            for (int eid = 0; eid < H.m; eid++)
            {
                edge_size.Add(eid, H.edge_ID_rev[eid].Count());
            }

            double min_conductance = double.MaxValue;
            var myA = new List<double>() { 1e-5 };
            var M = Matrix<double>.Build;
            var temp_vert = M.Dense(n, 4);
            for (int i = 0; i < n; i++)
            {
                var weight = H.vertex_weights[i];
                if (weight > 0)
                {
                    temp_vert[i, 2] = 1.0 / weight;
                }
            }

            var temp_edge = M.Dense(m, 3);
            for (int i = 0; i < m; i++)
            {
                var weight = edge_size[i];
                if (weight > 0)
                {
                    temp_edge[i, 1] = 1.0 / weight;
                }
            }

            foreach (double alpha in myA)
            {
                var vec = HyperGraph.A_HyperCut(H, temp_vert, temp_edge, v_init, alpha);

                for (int i = 0; i < n; i++)
                {
                    vec[i] /= H.vertex_weights[i];
                }
                int[] index = Enumerable.Range(0, n).ToArray<int>();
                Array.Sort<int>(index, (a, b) => vec[a].CompareTo(vec[b]));

                Array.Reverse(index);

                double vol_V = 0;
                for (int i = 0; i < n; i++) vol_V += H.vertex_weights[i];
                var num_contained_nodes = new Dictionary<int, int>();
                for (int eid = 0; eid < H.m; eid++)
                {
                    num_contained_nodes.Add(eid, 0);
                }

                double cut_val = 0;
                double vol_S = 0;
                double conductance = double.MaxValue;
                int best_index = -1;

                foreach (int i in index)
                {
                    vol_S += H.w_Degree(i);
                    if (vol_S <= vol_V / 2.0)
                    {
                        foreach (var e in H.incident_edges[i])
                        {
                            if (num_contained_nodes[e] == 0)
                            {
                                cut_val += H.edge_weights[e];
                            }
                            if (num_contained_nodes[e] == edge_size[e] - 1)
                            {
                                cut_val -= H.edge_weights[e];
                            }
                            num_contained_nodes[e] += 1;
                        }
                        conductance = cut_val / Math.Min(vol_S, vol_V - vol_S);
                        if (conductance < min_conductance)
                        {
                            min_conductance = conductance;
                            best_index = i;
                        }
                    }
                    else
                    {
                        break;
                    }
                }
            }
            time.Stop();
            TimeSpan ts = time.Elapsed;
            File.AppendAllText(out_file, v_init + "\t" + min_conductance + "\t" + time.ElapsedMilliseconds / 1000.0 + "\n");
            List<double> answer = new()
            {
                min_conductance,
                time.ElapsedMilliseconds / 1000.0
            };
            return answer;
        }

        public static List<double> Star_Expansion(HyperGraph H, string out_file, int v_init, double epsilon, int m, int n)
        {
            var time = new System.Diagnostics.Stopwatch();
            time.Start();
            var G = Graph.CreateStar(H);
            var M = Matrix<double>.Build;
            var temp_vert = M.Dense(m + n, 4);
            for (int i = 0; i < m + n; i++)
            {
                var weight = G.w_Degree(i);
                if (weight > 0)
                {
                    temp_vert[i, 3] = 1.0 / weight;
                }
            }
            var edge_size = new Dictionary<int, int>();
            for (int eid = 0; eid < H.m; eid++)
            {
                edge_size.Add(eid, H.edge_ID_rev[eid].Count);
            }
            double min_conductance = double.MaxValue;
            var A_cand = new List<double>() { 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95 };
            foreach (double alpha in A_cand)
            {
                var vec = Graph.PowerStar(G, temp_vert, v_init, alpha, m, n, epsilon);
                for (int i = 0; i < n; i++)
                {
                    vec[i] /= H.vertex_weights[i];
                }
                int[] index = Enumerable.Range(0, n).ToArray<int>();
                Array.Sort<int>(index, (a, b) => vec[a].CompareTo(vec[b]));

                Array.Reverse(index);

                double vol_V = 0;
                for (int i = 0; i < n; i++) vol_V += H.vertex_weights[i];
                var num_contained_nodes = new Dictionary<int, int>();
                for (int eid = 0; eid < H.m; eid++)
                {
                    num_contained_nodes.Add(eid, 0);
                }

                double cut_val = 0;
                double vol_S = 0;
                double conductance = double.MaxValue;
                int best_index = -1;

                foreach (int i in index)
                {
                    vol_S += H.w_Degree(i);
                    if (vol_S <= vol_V / 2.0)
                    {
                        foreach (var e in H.incident_edges[i])
                        {
                            if (num_contained_nodes[e] == 0)
                            {
                                cut_val += H.edge_weights[e];
                            }
                            if (num_contained_nodes[e] == edge_size[e] - 1)
                            {
                                cut_val -= H.edge_weights[e];
                            }
                            num_contained_nodes[e] += 1;
                        }
                        conductance = cut_val / Math.Min(vol_S, vol_V - vol_S);
                        if (conductance < min_conductance)
                        {
                            min_conductance = conductance;
                            best_index = i;
                        }
                    }
                    else
                    {
                        break;
                    }
                }
            }
            time.Stop();
            TimeSpan ts = time.Elapsed;
            File.AppendAllText(out_file, v_init + "\t" + min_conductance + "\t" + time.ElapsedMilliseconds / 1000.0 + "\n");
            List<double> answer = new()
            {
                min_conductance,
                time.ElapsedMilliseconds / 1000.0
            };
            return answer;
        }


        public static List<double> Clique_Expansion(HyperGraph H, string out_file, int v_init, double epsilon, int m, int n)
        {
            var time = new System.Diagnostics.Stopwatch();
            time.Start();
            var G = Graph.CreateClique(H);
            var M = Matrix<double>.Build;
            var temp_vert = M.Dense(n, 4);
            for (int i = 0; i < n; i++)
            {
                var weight = G.w_Degree(i);
                if (weight > 0)
                {
                    temp_vert[i, 3] = 1.0 / weight;
                }
            }
            var edge_size = new Dictionary<int, int>();
            for (int eid = 0; eid < H.m; eid++)
            {
                edge_size.Add(eid, H.edge_ID_rev[eid].Count);
            }
            double min_conductance = double.MaxValue;
            var A_cand = new List<double>() { 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95 };
            foreach (double alpha in A_cand)
            {
                var vec = Graph.PowerClique(G, temp_vert, v_init, alpha, n, epsilon);
                for (int i = 0; i < n; i++)
                {
                    vec[i] /= H.vertex_weights[i];
                }
                int[] index = Enumerable.Range(0, n).ToArray<int>();
                Array.Sort<int>(index, (a, b) => vec[a].CompareTo(vec[b]));

                Array.Reverse(index);

                double vol_V = 0;
                for (int i = 0; i < n; i++) vol_V += H.vertex_weights[i];
                var num_contained_nodes = new Dictionary<int, int>();
                for (int eid = 0; eid < H.m; eid++)
                {
                    num_contained_nodes.Add(eid, 0);
                }

                double cut_val = 0;
                double vol_S = 0;
                double conductance = double.MaxValue;
                int best_index = -1;

                foreach (int i in index)
                {
                    vol_S += H.w_Degree(i);
                    if (vol_S <= vol_V / 2.0)
                    {
                        foreach (var e in H.incident_edges[i])
                        {
                            if (num_contained_nodes[e] == 0)
                            {
                                cut_val += H.edge_weights[e];
                            }
                            if (num_contained_nodes[e] == edge_size[e] - 1)
                            {
                                cut_val -= H.edge_weights[e];
                            }
                            num_contained_nodes[e] += 1;
                        }
                        conductance = cut_val / Math.Min(vol_S, vol_V - vol_S);
                        if (conductance < min_conductance)
                        {
                            min_conductance = conductance;
                            best_index = i;
                        }
                    }
                    else
                    {
                        break;
                    }
                }
            }
            time.Stop();
            TimeSpan ts = time.Elapsed;
            File.AppendAllText(out_file, v_init + "\t" + min_conductance + "\t" + time.ElapsedMilliseconds / 1000.0 + "\n");
            List<double> answer = new()
            {
                min_conductance,
                time.ElapsedMilliseconds / 1000.0
            };
            return answer;
        }
        public static void Main(string[] args)
        {
            string in_path = "../../../instance/";
            string out_path = "../../../results/";

            List<string> datasets = new List<string> {"graphprod_LCC.txt"};
            //List<string> datasets = new List<string> {"netscience_LCC.txt"};
            //List<string> datasets = new List<string> {"dblp_kdd_LCC.txt"};
            //List<string> datasets = new List<string> {"opsahl-collaboration_LCC.txt"};
            //List<string> datasets = new List<string> {"dbpedia-writer_LCC.txt"};
            //List<string> datasets = new List<string> {"youtube-groupmemberships_LCC.txt"};
            //List<string> datasets = new List<string> {"dbpedia-recordlabel_LCC.txt"};
            //List<string> datasets = new List<string> {"dbpedia-genre_LCC.txt"};
            //List<string> datasets = new List<string> {"coauth-MAG-History.txt"};
            //List<string> datasets = new List<string> {"coauth-MAG-Geology.txt"};
            //List<string> datasets = new List<string> {"coauth-DBLP-full.txt"};
            //List<string> datasets = new List<string> {"threads-stack-overflow-full.txt"};

            //List<string> datasets = new List<string> {"graphprod_LCC.txt", "netscience_LCC.txt", "dblp_kdd_LCC.txt", "opsahl-collaboration_LCC.txt", "dbpedia-writer_LCC.txt", "youtube-groupmemberships_LCC.txt", "dbpedia-recordlabel_LCC.txt", "dbpedia-genre_LCC.txt", "coauth-MAG-History.txt", "coauth-MAG-Geology.txt", "coauth-DBLP-full.txt", "threads-stack-overflow-full.txt"};
            //List<int> number_of_nodes = new List<int> {86, 330, 5590, 13861, 54909, 88490, 158385, 253968, 1034876, 1261130, 1930379, 3455074};
            //List<string> methods = new List<string> { "Clique" };
            //List<string> methods = new List<string> { "Star" };
            //List<string> methods = new List<string> { "LocalClustering" };
            //List<string> methods = new List<string> { "HyperCut" };
            List<string> methods = new() { "Clique", "Star", "LocalClustering", "HyperCut" };

            const double epsilon = 1e-8;
            int dataset_index = 0;
            List<double> temp_answer = new List<double>();
            List<double> HyperCut_conductance = new List<double>();
            List<double> HyperCut_time = new List<double>();
            List<double> Local_conductance = new List<double>();
            List<double> Local_time = new List<double>();
            List<double> Clique_conductance = new List<double>();
            List<double> Clique_time = new List<double>();
            List<double> Star_conductance = new List<double>();
            List<double> Star_time = new List<double>();
            Random ranndom = new Random(0);
            foreach (var dataset in datasets)
            {
                var H = HyperGraph.Open(in_path + dataset);
                int n = H.num_nodes;
                int m = H.num_edges;
                double avg_edge_size = 0.0;
                foreach (var edge in H.edges)
                {
                    avg_edge_size += edge.Count;
                }
                avg_edge_size /= m;
                Console.WriteLine("dataset: " + dataset);
                List<int> seeds = new List<int>();
                while (seeds.Count < 50)
                {
                    var seed = ranndom.Next(0, n);
                    if (!seeds.Contains(seed))
                    {
                        seeds.Add(seed);
                    }
                }
                string list = string.Empty;
                foreach (var seed in seeds)
                {
                    list += seed + "\t";
                }
                if ((dataset == "dbpedia-recordlabel_LCC.txt") || (dataset == "dbpedia-genre_LCC.txt"))
                {
                    var out_file = File.Create(out_path + "alpha_" + dataset);
                    out_file.Close();
                    foreach (var seed in seeds)
                    {
                        Alpha_Average_Clustering(H, out_path + "alpha_" + dataset, seed);
                    }
                }
                foreach (var method in methods)
                {
                    var out_file = File.Create(out_path + method + "_" + dataset);
                    out_file.Close();
                    File.AppendAllText(out_path + method + "_" + dataset, "# nodes: " + n + "\n# edges: " + m + "\naverage degree: " + H.vertex_weights.Sum() / n + "\naverage edge size: " + avg_edge_size + "\n\nseed vertex\t" + "conductance\t" + "time taken (seconds)" + "\n\n");

                    foreach (var seed in seeds)
                    {

                        if (method == "Star")
                        {
                            temp_answer = Star_Expansion(H, out_path + "Star_" + dataset, seed, epsilon, m, n);
                            Star_conductance.Add(temp_answer[0]);
                            Star_time.Add(temp_answer[1]);
                        }

                        if (method == "Clique")
                        {
                            temp_answer = Clique_Expansion(H, out_path + "Clique_" + dataset, seed, epsilon, m, n);
                            Clique_conductance.Add(temp_answer[0]);
                            Clique_time.Add(temp_answer[1]);
                        }

                        if (method == "LocalClustering")
                        {
                            temp_answer = Proposed_local_round(H, out_path + "LocalClustering_" + dataset, seed);
                            Local_conductance.Add(temp_answer[0]);
                            Local_time.Add(temp_answer[1]);
                        }

                        if (method == "HyperCut")
                        {
                            temp_answer = Average_Clustering(H, out_path + "HyperCut_" + dataset, seed, epsilon);
                            HyperCut_conductance.Add(temp_answer[0]);
                            HyperCut_time.Add(temp_answer[1]);
                        }                        
                    }
                    
                    if (method == "Star")
                    {
                        File.AppendAllText(out_path + "Star_" + dataset, "\n\nConductance\n");
                        for (int i = 0; i < Star_conductance.Count; i++)
                        {
                            File.AppendAllText(out_path + "Star_" + dataset, Star_conductance[i] + ", ");
                        }
                        File.AppendAllText(out_path + "Star_" + dataset, "\n\ntime taken (seconds)\n");
                        for (int i = 0; i < Star_time.Count; i++)
                        {
                            File.AppendAllText(out_path + "Star_" + dataset, Star_time[i] + ", ");
                        }
                    }

                    if (method == "Clique")
                    {
                        File.AppendAllText(out_path + "Clique_" + dataset, "\n\nConductance\n");
                        for (int i = 0; i < Clique_conductance.Count; i++)
                        {
                            File.AppendAllText(out_path + "Clique_" + dataset, Clique_conductance[i] + ", ");
                        }
                        File.AppendAllText(out_path + "Clique_" + dataset, "\n\ntime taken (seconds)\n");
                        for (int i = 0; i < Clique_time.Count; i++)
                        {
                            File.AppendAllText(out_path + "Clique_" + dataset, Clique_time[i] + ", ");
                        }
                    }

                    if (method == "LocalClustering")
                    {
                        File.AppendAllText(out_path + "LocalClustering_" + dataset, "\n\nConductance\n");
                        for (int i = 0; i < Local_conductance.Count; i++)
                        {
                            File.AppendAllText(out_path + "LocalClustering_" + dataset, Local_conductance[i] + ", ");
                        }
                        File.AppendAllText(out_path + "LocalClustering_" + dataset, "\n\ntime taken (seconds)\n");
                        for (int i = 0; i < Local_time.Count; i++)
                        {
                            File.AppendAllText(out_path + "LocalClustering_" + dataset, Local_time[i] + ", ");
                        }
                    }

                    if (method == "HyperCut")
                    {
                        File.AppendAllText(out_path + "HyperCut_" + dataset, "\n\nConductance\n");
                        for (int i = 0; i < HyperCut_conductance.Count; i++)
                        {
                            File.AppendAllText(out_path + "HyperCut_" + dataset, HyperCut_conductance[i] + ", ");
                        }
                        File.AppendAllText(out_path + "HyperCut_" + dataset, "\n\ntime taken (seconds)\n");
                        for (int i = 0; i < HyperCut_time.Count; i++)
                        {
                            File.AppendAllText(out_path + "HyperCut_" + dataset, HyperCut_time[i] + ", ");
                        }

                    }
                    Console.WriteLine("\n");
                }
                dataset_index++;
            }
        }
    }
}
