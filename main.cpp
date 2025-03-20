// 1次ホモロジー基底の計算

#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <Eigen/Core>

struct HasherEigenVector2i
{
    std::size_t operator()(Eigen::Vector2i const& key) const
    {
        std::hash<int> hasher;
        std::size_t seed = 0;
        seed ^= hasher(key(0)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hasher(key(1)) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

struct Halfedge
{
    int idx;
    int h_opp;
    std::vector<Eigen::Vector3i> &F;

    Halfedge(std::vector<Eigen::Vector3i> &f) : F(f) {}

    int face()
    {
        return int(idx/3);
    }
    int h_next()
    {
        if(idx % 3 == 2) return idx-2;
        else return idx+1;
    }
    int h_prev()
    {
        if(idx % 3 == 0) return idx+2;
        else return idx-1;
    }
    int v_src()
    {
        return F[face()][(idx+1)%3];
    }

    int v_tgt()
    {
        return F[face()][(idx+2)%3];
    }
};

struct Mesh
{
    std::vector<Eigen::Vector3d> V;
    std::vector<Eigen::Vector3i> F;
    std::vector<Eigen::Vector2i> E;
    std::map<std::pair<int, int>, int> edge_map;

    std::vector<Halfedge> HE;
    std::map<std::pair<int, int>, int> he_map;
    std::vector<int> h_out;

    void make_halfedge_list()
    {
        for (int i = 0; i < F.size(); i++)
        {
            for (int j = 0; j < 3; j++)
            {
                Halfedge h(F);
                h.idx = 3*i+j;

                auto key = std::make_pair(h.v_src(), h.v_tgt());
                auto keyswap = std::make_pair(h.v_tgt(), h.v_src());
                if (he_map.contains(keyswap))
                {
                    h.h_opp = he_map.at(keyswap);
                    HE[he_map.at(keyswap)].h_opp = 3*i+j;
                }
                else
                {
                    h.h_opp = -1;
                    he_map.emplace(key, 3*i+j);
                }
                HE.push_back(h);
            }
        }
        // h_out を計算
        h_out.resize(V.size());
        for (int i = 0; i < HE.size(); i++)
        {
            //  h_out に境界半辺が保存されていない場合のみ更新
            if (HE[ h_out[HE[i].v_src()] ].h_opp != -1)
            {
                h_out[HE[i].v_src()] = i;
            }
        }
        // 境界半辺の h_opp に次の境界半辺を保存する
        /*for (int i=0; i<h_out.size(); i++) {
            if(HE[h_out[i]].h_opp == -1){
                HE[h_out[i]].h_opp = -h_out[ HE[h_out[i]].v_tgt(F) ] - 1;
                // h_i1.h_opp = -i2-1;
            }
        }*/
    }

    void make_edge_map()
    {
        int edge_idx = 0;
        for (auto& f : F)
        {
            for (int j = 0; j < 3; j++)
            {
                auto edge = std::make_pair(f(j), f((j+1)%3));
                auto edge_opp = std::make_pair(f((j+1)%3), f(j));
                if (! edge_map.contains(edge) && ! edge_map.contains(edge_opp))
                {
                    edge_map.emplace(edge, edge_idx);
                    edge_idx ++;
                }
            }
        }
        E.resize(edge_map.size());
        for (const auto& [edge, idx] : edge_map)
        {
            E[idx](0) = edge.first;
            E[idx](1) = edge.second;
        }
    }
};

struct UnionFind
{
    std::vector<int> parents;

    UnionFind(int n) : parents(n, -1) {}

    int find(int x)
    {
        int rank = 0;
        while (parents[x] >= 0)
        {
            x = parents[x];
            rank ++;
        }
        parents[x] = -(rank + 1);

        return x;
    }

    void unite(int x, int y)
    {
        int rootx = find(x);
        int rooty = find(y);

        if (rootx > rooty) parents[rootx] = rooty;
        else if (rootx < rooty) parents[rooty] = rootx;
    }
};

class TreeCotree
{

public:
    std::vector<Eigen::Vector2i> T;     // 全域木
    std::vector<Eigen::Vector2i> C;     // Cotree
    std::vector<Eigen::Vector2i> X;     // 残りの辺集合
//    std::vector<Eigen::Vector2i> non_separating_cycle;
    std::vector<std::vector<Eigen::Vector2i>> non_separating_cycle;     // 2g個の非分離サイクル

    // 全域木の計算
    std::vector<Eigen::Vector2i> kruskal(int num_vertices,
                                         std::vector<Eigen::Vector2i> const& E)
    {
        UnionFind uf(num_vertices);
        std::vector<Eigen::Vector2i> T;     // 全域木

        for (const auto& edge : E)
        {
            if ( uf.find(edge(0)) != uf.find(edge(1)) )
            {
                uf.unite(edge(0), edge(1));
                T.push_back(edge);
            }
        }
        return T;
    }

    // 双対グラフの辺から主グラフの辺へのマップを作る
    std::map<std::pair<int, int>, Eigen::Vector2i> make_dual_edge_map(std::vector<Eigen::Vector2i> const& E,
                            std::vector<Halfedge>& HE,
                            std::map<std::pair<int, int>, int> const& he_map)
    {
        std::map<std::pair<int, int>, Eigen::Vector2i> dual_edge_map;
        for (int i = 0; i < HE.size(); i++)
        {
            auto key = std::make_pair(HE[i].face(), HE[ HE[i].h_opp ].face());
            auto key_swap = std::make_pair(HE[ HE[i].h_opp ].face(), HE[i].face());
            if ( ! dual_edge_map.contains(key_swap) )
            {
                Eigen::Vector2i edge;
                edge(0) = HE[i].v_src();
                edge(1) = HE[i].v_tgt();
                dual_edge_map.emplace(key, edge);
            }
        }
        return dual_edge_map;
    }

    // Cotree を計算
    void compute_C(int f_num,
                   std::vector<Eigen::Vector2i> const& E,
                   std::vector<Halfedge>& HE,
                   std::map<std::pair<int, int>, int> const& he_map)
    {
        // E\T を計算
        std::vector<Eigen::Vector2i> E_T;       // E\T
        std::unordered_set<Eigen::Vector2i, HasherEigenVector2i> T_set;
        for (const auto& tree : T)
        {
            T_set.emplace(tree);
        }
        for (const auto& edge : E)
        {
            if ( ! T_set.contains(edge) )
            {
                E_T.push_back(edge);
            }
        }

        // (E\T)*を計算
        std::vector<Eigen::Vector2i> dualE_T;   // (E\T)*
        for (const auto& et : E_T)
        {
            Eigen::Vector2i co_edge;
            auto key = std::make_pair(et(0), et(1));
            auto key_swap = std::make_pair(et(1), et(0));
            if (he_map.contains(key))
            {
                co_edge(0) = HE[ he_map.at(key) ].face();
                co_edge(1) = HE[ HE[ he_map.at(key) ].h_opp ].face();
            }
            else
            {
                co_edge(0) = HE[ he_map.at(key_swap) ].face();
                co_edge(1) = HE[ HE[ he_map.at(key_swap) ].h_opp ].face();
            }
            dualE_T.push_back(co_edge);
        }

        // C*を計算
        std::vector<Eigen::Vector2i> dualC = kruskal(f_num, dualE_T);

        // Cを計算
        std::map<std::pair<int, int>, Eigen::Vector2i> dual_edge_map = make_dual_edge_map(E, HE, he_map);
        for (const auto& dc : dualC)
        {
            auto key = std::make_pair(dc(0), dc(1));
            auto key_swap = std::make_pair(dc(1), dc(0));
            if (dual_edge_map.contains(key))
            {
                C.push_back(dual_edge_map.at(key));
            }
            else
            {
                C.push_back(dual_edge_map.at(key_swap));
            }
        }
    }

    // 残りの辺集合Xを計算
    void compute_X(std::vector<Eigen::Vector2i> const& E)
    {
        // T ∪ C を計算
        std::unordered_set<Eigen::Vector2i, HasherEigenVector2i> TUC_set;
        for (const auto& tree : T)
        {
            TUC_set.emplace(tree);
        }
        for (const auto& cotree : C)
        {
            TUC_set.emplace(cotree);
        }
        // 残りの辺集合Xを計算
        for (const auto& edge : E)
        {
            if ( ! TUC_set.contains(edge) )
            {
                X.push_back(edge);
            }
        }
    }

    // 深さ優先探索でcurrentからtargetへのパスを計算
    bool find_path_DFS(int current,
                       int target,
                       std::unordered_map<int, std::vector<int>> const& adj_tree_list,
                       std::vector<int>& path,
                       std::vector<bool>& visited)
    {
        visited[current] = true;
        path.push_back(current);

        if (current == target) return true;

        for (auto neighbor : adj_tree_list.at(current))
        {
            if (!visited[neighbor])
            {
                if (find_path_DFS(neighbor, target, adj_tree_list, path, visited))
                {
                    return true;
                }
            }
        }
        path.pop_back();
        return false;
    }

    // x_edgeの基本サイクルを計算
    std::vector<Eigen::Vector2i> compute_basis_sycle(std::unordered_map<int, std::vector<int>> adj_tree_list,
                                                     Eigen::Vector2i const& x_edge)
    {
        // 木内の u から v へのパスを計算
        std::vector<int> path;                          // u から v へのパス
        std::vector<bool> visited(adj_tree_list.size(), false); // 訪問フラグ
        find_path_DFS(x_edge(0), x_edge(1), adj_tree_list, path, visited);

        // サイクルを構成
        std::vector<Eigen::Vector2i> cycle;
        for (int i = 0; i < path.size() - 1; i++)
        {
            cycle.emplace_back(path[i], path[i + 1]);
        }
        cycle.push_back(x_edge);

        return cycle;
    }

    // 2g個の非分離サイクルを計算
    void compute_non_separating_cycle(int v_num,
                                      int f_num,
                                      std::vector<Eigen::Vector2i> const& E,
                                      std::vector<Halfedge>& HE,
                                      std::map<std::pair<int, int>, int> const& he_map)
    {
        // 全域木Tを計算
        T = kruskal(v_num, E);

        // G*の全域木Cを計算
        compute_C(f_num, E, HE, he_map);

        // 残りの辺集合Xを計算
        compute_X(E);

        // 非分離サイクル（Xの基本サイクル）を計算
        // 全域木を隣接リストに変換
        std::unordered_map<int, std::vector<int>> adj_tree_list;
        for (const auto& edge : T)
        {
            adj_tree_list[edge(0)].push_back(edge(1));
            adj_tree_list[edge(1)].push_back(edge(0));
        }
        for (auto& x_edge : X)
        {
            std::vector<Eigen::Vector2i> nsc = compute_basis_sycle(adj_tree_list, x_edge);
//            non_separating_cycle.insert(non_separating_cycle.end(), nsc.begin(), nsc.end());
            non_separating_cycle.push_back(nsc);
        }


    }
};

void export_obj(std::string name, std::string type, std::vector<Eigen::Vector3d> const& vert, std::vector<Eigen::Vector2i> const& tree)
{
    std::ofstream of;
    name.erase(name.length()-4);
    std::string filename = name + "_" + type + ".obj";
    of.open(filename, std::ios::out);
    for(auto& v : vert)
    {
        of << "v " << v(0) << " " << v(1) << " " << v(2) << std::endl;
    }
    for(auto& t : tree)
    {
        of << "l " << t(0)+1 << " " << t(1)+1 << std::endl;
    }
    of.close();
}

void read_obj(std::string const& filename, Mesh& mesh)
{
    std::ifstream ifs(filename);
    if (ifs.fail())
    {
        std::cerr << "Failed to open file." << "\n";
        std::exit(1);
    }

    std::string line;
    while (std::getline(ifs, line))
    {
        std::istringstream string_in(line);
        std::string type;
        string_in >> type;

        if (type == "v")
        {
            Eigen::Vector3d v;
            string_in >> v(0) >> v(1) >> v(2);
            mesh.V.push_back(v);
        }
        else if (type == "f")
        {
            Eigen::Vector3i f;
            string_in >> f(0) >> f(1) >> f(2);
            f -= Eigen::Vector3i{1, 1, 1};
            mesh.F.push_back(f);
        }
    }
}


int main(int argc, const char * argv[])
{
    // Read mesh file...........................................
    std::string filename;
    if (argc != 2)
    {
        std::cout << "wrong command line argument" << std::endl;
        std::exit(1);
    }
    filename = std::string(argv[1]);

    auto start = std::chrono::system_clock::now();

    Mesh mesh;
    read_obj(filename, mesh);
    mesh.make_halfedge_list();
    mesh.make_edge_map();

    // Tree-Cotree 分解....................................................................................
    TreeCotree tc;
    tc.compute_non_separating_cycle((int)mesh.V.size(), (int)mesh.F.size(), mesh.E, mesh.HE, mesh.he_map);

    auto end = std::chrono::system_clock::now();

    // 出力....................................................................................
//    export_obj(filename, "ST", mesh.V, tc.T);
//    export_obj(filename, "C", mesh.V, tc.C);
//    export_obj(filename, "X", mesh.V, tc.X);
//    export_obj(filename, "nsc", mesh.V, tc.non_separating_cycle);
    for (int i = 0; i < tc.non_separating_cycle.size(); i++)
    {
        export_obj(filename, "nsc" + std::to_string(i+1), mesh.V, tc.non_separating_cycle[i]);
    }

    using namespace std::chrono_literals;
    std::cout << "Execute time\n";
    std::cout << (end - start) / 1.0s << " s\n";

    return 0;
}
