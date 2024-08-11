using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;

namespace HyperGraphCLustering
{
    public class HyperGraph
    {
        public int num_nodes;
        public int num_edges;
        public List<List<int>> edges = new List<List<int>>();
        public List<double> vertex_weights = new List<double>();
        public List<List<int>> incident_edges = new List<List<int>>();
        public List<double> edge_weights = new List<double>();
        public Dictionary<List<int>, int> edge_ID = new Dictionary<List<int>, int>();
        public Dictionary<int, List<int>> edge_ID_rev = new Dictionary<int, List<int>>();

        public double w_Degree(int v)
        {
            double sum = 0;
            foreach (int e in incident_edges[v])
            {
                sum += edge_weights[e];
            }
            return sum;
        }

        public int n
        {
            get
            {
                return incident_edges.Count;
            }
        }

        public int m
        {
            get
            {
                return edges.Count;
            }
        }

        public static HyperGraph Open(string fn)
        {
            var fs = new FileStream(fn, FileMode.Open);
            var sr = new StreamReader(fs);

            var edges = new List<List<int>>();
            var edge_weights = new List<double>();
            int vertex_num = 0;
            for (string line; (line = sr.ReadLine()) != null;)
            {
                var words = line.Split();
                var edge = new List<int>();
                int i = 0;
                foreach (var word in words)
                {
                    if (i < words.Length - 1)
                    {
                        int v = int.Parse(word);
                        edge.Add(v);
                        if (v >= vertex_num) vertex_num = v + 1;
                        i++;
                    }
                }
                edges.Add(edge);
                edge_weights.Add(double.Parse(words.Last()));
            }

            var H = new HyperGraph();
            H.num_nodes = vertex_num;
            H.num_edges = edges.Count;
            H.edges = edges;
            H.edge_weights = edge_weights;
            for (int v = 0; v < vertex_num; v++)
            {
                H.incident_edges.Add(new List<int>());
                H.vertex_weights.Add(0.0);
            }
            for (int i = 0; i < edges.Count; i++)
            {
                var edge = edges[i];
                double weight = edge_weights[i];
                foreach (var v in edge)
                {
                    H.incident_edges[v].Add(i);
                    H.vertex_weights[v] += weight;
                }
                H.edge_ID.Add(edge, i);
                H.edge_ID_rev.Add(i, edge);
            }
            fs.Close();
            return H;
        }

        public Vector<double> T(Vector<double> vec, int v_init, double alpha)
        {
            var res = CreateVector.Dense<double>(n);

            const double eps = 1e-8;
            foreach (var edge in edges)
            {
                var argmaxs = new List<int>();
                var argmins = new List<int>();
                double maxval = double.MinValue, minval = double.MaxValue;
                foreach (var v in edge)
                {
                    var val = vec[v] / w_Degree(v);
                    if (val > maxval + eps)
                    {
                        maxval = val;
                        argmaxs.Clear();
                        argmaxs.Add(v);
                    }
                    else if (val > maxval - eps)
                    {
                        argmaxs.Add(v);
                    }

                    if (val < minval - eps)
                    {
                        minval = val;
                        argmins.Clear();
                        argmins.Add(v);
                    }
                    else if (val < minval + eps)
                    {
                        argmins.Add(v);
                    }
                }
                foreach (var v in argmaxs)
                {
                    res[v] += edge_weights[edge_ID[edge]] * (maxval - minval) / argmaxs.Count;
                }
                foreach (var v in argmins)
                {
                    res[v] -= edge_weights[edge_ID[edge]] * (maxval - minval) / argmins.Count;
                }
            }

            var res_init = vec;
            res_init[v_init] -= 1;

            var mix = (1 - alpha) * res + alpha * res_init;

            return mix;
        }

        public static Vector<double> Iterate(HyperGraph H, Vector<double> vec, int v_init, double dt, double alpha)
        {
            var dv = H.T(vec, v_init, alpha);
            var res = vec;
            res -= dv * dt;
            return res;
        }

        public static Vector<double> Simulate(HyperGraph H, Vector<double> vec, int v_init, double dt, double T, double alpha)
        {
            var cur_time = 0.0;
            while (cur_time < T)
            {
                var next_time = Math.Min(cur_time + dt, T);
                vec = Iterate(H, vec, v_init, next_time - cur_time, alpha);
                cur_time = next_time;
            }
            return vec;
        }

        public Vector<double> T_round(Vector<double> vec, int v_init, double alpha, List<List<int>> active_edges)
        {
            var res = CreateVector.Dense<double>(n);

            const double eps = 1e-8;

            foreach (var edge in active_edges)
            {
                var argmaxs = new List<int>();
                var argmins = new List<int>();
                double maxval = double.MinValue, minval = double.MaxValue;
                foreach (var v in edge)
                {
                    var val = vec[v] / w_Degree(v);
                    if (val > maxval + eps)
                    {
                        maxval = val;
                        argmaxs.Clear();
                        argmaxs.Add(v);
                    }
                    else if (val > maxval - eps)
                    {
                        argmaxs.Add(v);
                    }

                    if (val < minval - eps)
                    {
                        minval = val;
                        argmins.Clear();
                        argmins.Add(v);
                    }
                    else if (val < minval + eps)
                    {
                        argmins.Add(v);
                    }
                }
                foreach (var v in argmaxs)
                {
                    res[v] += edge_weights[edge_ID[edge]] * (maxval - minval) / argmaxs.Count;
                }
                foreach (var v in argmins)
                {
                    res[v] -= edge_weights[edge_ID[edge]] * (maxval - minval) / argmins.Count;
                }
            }

            var res_init = CreateVector.Dense<double>(n);
            vec.CopyTo(res_init);
            res_init[v_init] -= 1;

            var mix = (1 - alpha) * res + alpha * res_init;

            return mix;
        }

        public static Vector<double> Iterate_round(HyperGraph H, Vector<double> vec, int v_init, double dt, double alpha, List<List<int>> active_edges, List<int> active)
        {
            var dv = H.T_round(vec, v_init, alpha, active_edges);
            var res = vec;
            res -= dv * dt;

            for (int i = 0; i < H.n; i++)
            {
                if (res[i] < 1e-5)
                {
                    res[i] = 0;
                }
            }

            var new_active_edges = new List<int>();

            for (int i = 0; i < H.n; i++)
            {
                if (vec[i] == 0 && res[i] != 0)
                {
                    foreach (var f in H.incident_edges[i])
                    {
                        if (active[f] == 0)
                        {
                            new_active_edges.Add(f);
                            active[f] = 1;
                        }
                    }
                }
            }

            foreach (var e in new_active_edges)
            {
                active_edges.Add(H.edge_ID_rev[e]);
            }

            return res;
        }

        public static Vector<double> Simulate_round(HyperGraph H, Vector<double> vec, int v_init, double dt, double T, double alpha)
        {
            var active_edges = new List<List<int>>();
            var active = new List<int>(new int[H.m]);

            foreach (var e in H.incident_edges[v_init])
            {
                active_edges.Add(H.edge_ID_rev[e]);
                active[e] = 1;
            }

            var cur_time = 0.0;
            while (cur_time < T)
            {
                var next_time = Math.Min(cur_time + dt, T);
                vec = Iterate_round(H, vec, v_init, next_time - cur_time, alpha, active_edges, active);
                cur_time = next_time;
            }
            return vec;
        }

        public static Vector<double> A_HyperCut(HyperGraph H, Matrix<double> vec_vertices, Matrix<double> vec_edges, int v_init, double alpha, double epsilon)
        {
            var pre_vec = vec_vertices.Column(0);
            pre_vec.Clear();
            var next_vec = pre_vec;
            next_vec[v_init] = 1.0;
            double beta = 1.0 - alpha;
            var init_row = vec_vertices.Row(v_init);
            init_row[0] = alpha;
            init_row[1] = beta;
            if (init_row[2] > 0)
            {
                init_row[3] = init_row[1] * init_row[2];
            }
            else
            {
                init_row[0] = 1.0;
                vec_vertices.SetRow(v_init, init_row);
                return vec_vertices.Column(0);
            }
            vec_vertices.SetRow(v_init, init_row);
            HashSet<int> active_vertices = new()
            {
                v_init
            };
            HashSet<int> active_edges = new();
            foreach (var active_vertex in active_vertices)
            {
                foreach (var edge_index in H.incident_edges[active_vertex])
                {
                    active_edges.Add(edge_index);
                }
            }
            var err_vec = next_vec;
            while (err_vec.L2Norm() > epsilon)
            {
                pre_vec = next_vec;
                foreach (var active_edge_index in active_edges)
                {
                    var active_edge = H.edge_ID_rev[active_edge_index];
                    foreach (var v in active_edge)
                    {
                        vec_edges[active_edge_index, 0] += vec_vertices[v, 3];
                        active_vertices.Add(v);
                    }
                }
                vec_edges.SetColumn(2, vec_edges.Column(0).PointwiseMultiply(vec_edges.Column(1)));
                vec_vertices.ClearColumn(1);
                foreach (var active_vertex in active_vertices)
                {
                    foreach (var incident_edge in H.incident_edges[active_vertex])
                    {
                        vec_vertices[active_vertex, 1] += vec_edges[incident_edge, 2];
                        active_edges.Add(incident_edge);
                    }
                }
                vec_edges.ClearColumns(0, 2);
                next_vec = vec_vertices.Column(0) + vec_vertices.Column(1);
                vec_vertices.SetColumn(0, vec_vertices.Column(0) + alpha * vec_vertices.Column(1));
                vec_vertices.SetColumn(1, beta * vec_vertices.Column(1));
                vec_vertices.SetColumn(3, vec_vertices.Column(1).PointwiseMultiply(vec_vertices.Column(2)));
                err_vec = next_vec - pre_vec;
            }
            return vec_vertices.Column(0) + vec_vertices.Column(1);
        }
    }





    public class Graph
    {
        public List<Dictionary<int, double>> adj_list = new List<Dictionary<int, double>>();

        public double w_Degree(int v)
        {
            double sum = 0;
            foreach (var neighbor_val in adj_list[v].Values)
            {
                sum += neighbor_val;
            }
            return sum;
        }

        public static Graph CreateStar(HyperGraph H)
        {
            int n = H.num_nodes;
            var G = new Graph();
            for (int i = 0; i < H.num_nodes + H.num_edges; i++)
            {
                G.adj_list.Add(new Dictionary<int, double>());
            }
            foreach (var edge in H.edges)
            {
                var edge_ID = H.edge_ID[edge];
                var edge_size = edge.Count;
                var edge_weight = H.edge_weights[edge_ID] / edge_size;
                for (int i = 0; i < edge_size; i++)
                {
                    G.adj_list[n + edge_ID][edge[i]] = edge_weight;
                    G.adj_list[edge[i]][n + edge_ID] = edge_weight;
                }
            }
            return G;
        }

        public static Vector<double> PowerStar(Graph G, Matrix<double> vec_vertices, int v_init, double alpha, int m, int n, double epsilon)
        {
            var pre_vec = CreateVector.Dense<double>(m + n);
            pre_vec[v_init] = 1.0;
            double beta = 1.0 - alpha;
            var init_row = vec_vertices.Row(v_init);
            init_row[0] = alpha + beta * 0.5;
            init_row[1] = beta * 0.5;
            if (init_row[3] > 0)
            {
                init_row[2] = init_row[1] * init_row[3];
            }
            else
            {
                return pre_vec.SubVector(0, n);
            }
            vec_vertices.SetRow(v_init, init_row);
            HashSet<int> active_vertices = new()
            {
                v_init
            };
            HashSet<int> temp_active_vertices = new();
            var nhbrs = G.adj_list[v_init].Keys;
            var measure = vec_vertices[v_init, 2];
            foreach (var nhbr in nhbrs)
            {
                temp_active_vertices.Add(nhbr);
                vec_vertices[nhbr, 0] += measure * G.adj_list[v_init][nhbr];
            }
            active_vertices.UnionWith(temp_active_vertices);
            temp_active_vertices.Clear();
            var err_vec = vec_vertices.Column(0).Subtract(pre_vec);
            while (err_vec.L1Norm() > epsilon)
            {
                pre_vec = vec_vertices.Column(0);
                vec_vertices.SetColumn(1, pre_vec * beta * 0.5);
                vec_vertices.SetColumn(0, vec_vertices.Column(1));
                vec_vertices[v_init, 0] += alpha;
                vec_vertices.SetColumn(2, vec_vertices.Column(1).PointwiseMultiply(vec_vertices.Column(3)));
                foreach (var vertex in active_vertices)
                {
                    measure = vec_vertices[vertex, 2];
                    nhbrs = G.adj_list[vertex].Keys;
                    foreach (var nhbr in nhbrs)
                    {
                        temp_active_vertices.Add(nhbr);
                        vec_vertices[nhbr, 0] += measure * G.adj_list[vertex][nhbr];
                    }
                }
                active_vertices.UnionWith(temp_active_vertices);
                temp_active_vertices.Clear();
                err_vec = vec_vertices.Column(0).Subtract(pre_vec);
            }
            return vec_vertices.Column(0).SubVector(0, n);
        }

        public static Graph CreateClique(HyperGraph H)
        {
            var G = new Graph();
            for (int i = 0; i < H.num_nodes; i++)
            {
                G.adj_list.Add(new Dictionary<int, double>());
            }
            foreach (var edge in H.edges)
            {
                var edge_weight = H.edge_weights[H.edge_ID[edge]];
                var edge_size = edge.Count;
                for (int i = 0; i < edge_size - 1; i++)
                {
                    var adj_list = G.adj_list[edge[i]];
                    for (int j = i + 1; j < edge_size; j++)
                    {
                        if (adj_list.ContainsKey(edge[j]))
                        {
                            G.adj_list[edge[i]][edge[j]] += edge_weight;
                            G.adj_list[edge[j]][edge[i]] += edge_weight;
                        }
                        else
                        {
                            G.adj_list[edge[i]][edge[j]] = edge_weight;
                            G.adj_list[edge[j]][edge[i]] = edge_weight;
                        }
                    }
                }
            }
            return G;
        }
        public static Vector<double> PowerClique(Graph G, Matrix<double> vec_vertices, int v_init, double alpha, int n, double epsilon)
        {
            var pre_vec = CreateVector.Dense<double>(n);
            pre_vec[v_init] = 1.0;
            double beta = 1.0 - alpha;
            var init_row = vec_vertices.Row(v_init);
            init_row[0] = alpha + beta * 0.5;
            init_row[1] = beta * 0.5;
            if (init_row[3] > 0)
            {
                init_row[2] = init_row[1] * init_row[3];
            }
            else
            {
                return pre_vec;
            }
            vec_vertices.SetRow(v_init, init_row);
            HashSet<int> active_vertices = new()
            {
                v_init
            };
            HashSet<int> temp_active_vertices = new();
            var nhbrs = G.adj_list[v_init].Keys;
            var measure = vec_vertices[v_init, 2];
            foreach (var nhbr in nhbrs)
            {
                temp_active_vertices.Add(nhbr);
                vec_vertices[nhbr, 0] += measure * G.adj_list[v_init][nhbr];
            }
            active_vertices.UnionWith(temp_active_vertices);
            temp_active_vertices.Clear();
            var err_vec = vec_vertices.Column(0).Subtract(pre_vec);
            while (err_vec.L1Norm() > epsilon)
            {
                pre_vec = vec_vertices.Column(0);
                vec_vertices.SetColumn(1, pre_vec * beta * 0.5);
                vec_vertices.SetColumn(0, vec_vertices.Column(1));
                vec_vertices[v_init, 0] += alpha;
                vec_vertices.SetColumn(2, vec_vertices.Column(1).PointwiseMultiply(vec_vertices.Column(3)));
                foreach (var vertex in active_vertices)
                {
                    measure = vec_vertices[vertex, 2];
                    nhbrs = G.adj_list[vertex].Keys;
                    foreach (var nhbr in nhbrs)
                    {
                        temp_active_vertices.Add(nhbr);
                        vec_vertices[nhbr, 0] += measure * G.adj_list[vertex][nhbr];
                    }
                }
                active_vertices.UnionWith(temp_active_vertices);
                temp_active_vertices.Clear();
                err_vec = vec_vertices.Column(0).Subtract(pre_vec);
            }
            return vec_vertices.Column(0);
        }

    }

}
