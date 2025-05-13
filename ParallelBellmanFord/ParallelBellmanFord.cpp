#include <iostream>
#include <vector>
#include <limits>
#include <random>
#include <chrono>
#include <fstream>
#include <algorithm>
#include <omp.h>
#include <unordered_map>

using namespace std;
using namespace std::chrono;

const int INF = numeric_limits<int>::max();

struct Edge {
    int u, v;
    int weight;
};

struct Graph {
    int V;
    vector<Edge> edges;
    vector<vector<pair<int, int>>> adj; // adjacency list (neighbor, weight)

    Graph() : V(0) {}
    Graph(int vertices) : V(vertices), adj(vertices) {}

    void addEdge(int u, int v, int weight) {
        edges.push_back({ u, v, weight });
        adj[u].emplace_back(v, weight);
    }

    // Generate a complete graph with random weights more efficiently
    static Graph generateCompleteGraph(int V, int minWeight = -100, int maxWeight = 100) {
        Graph g(V);
        random_device rd;
        mt19937 gen(rd());
        uniform_int_distribution<> dist(minWeight, maxWeight);

        // Pre-allocate memory
        g.edges.reserve(V * (V - 1));

        for (int u = 0; u < V; ++u) {
            g.adj[u].reserve(V - 1);
            for (int v = 0; v < V; ++v) {
                if (u != v) {
                    int weight = dist(gen);
                    g.addEdge(u, v, weight);
                }
            }
        }
        return g;
    }
};

// Cache for generated graphs with LRU policy
const size_t MAX_CACHE_SIZE = 5;
unordered_map<int, Graph> graphCache;
vector<int> cacheAccessOrder;

struct BFResult {
    vector<int> distances;
    vector<int> predecessors;
    bool hasNegativeCycle;
    double executionTime;
};

// Optimized sequential Bellman-Ford using adjacency list
BFResult bellmanFordSequential(const Graph& graph, int source) {
    auto start = high_resolution_clock::now();

    const int V = graph.V;
    vector<int> dist(V, INF);
    vector<int> pred(V, -1);
    dist[source] = 0;

    bool changed;
    for (int i = 1; i <= V - 1; ++i) {
        changed = false;
        for (const Edge& edge : graph.edges) {
            if (dist[edge.u] != INF && dist[edge.u] + edge.weight < dist[edge.v]) {
                dist[edge.v] = dist[edge.u] + edge.weight;
                pred[edge.v] = edge.u;
                changed = true;
            }
        }
        if (!changed) break;
    }

    // Check for negative-weight cycles
    bool hasNegativeCycle = false;
    for (const Edge& edge : graph.edges) {
        if (dist[edge.u] != INF && dist[edge.u] + edge.weight < dist[edge.v]) {
            hasNegativeCycle = true;
            break;
        }
    }

    auto end = high_resolution_clock::now();
    double duration = duration_cast<milliseconds>(end - start).count() / 1000.0;

    return { dist, pred, hasNegativeCycle, duration };
}

// Optimized Parallel Bellman-Ford Version 1 - Edge-based parallelization
BFResult bellmanFordParallelV1(const Graph& graph, int source) {
    auto start = high_resolution_clock::now();

    const int V = graph.V;
    vector<int> dist(V, INF);
    vector<int> pred(V, -1);
    dist[source] = 0;

    const int numEdges = graph.edges.size();
    const int chunkSize = max(1000, numEdges / (omp_get_max_threads() * 10));

    bool changed;
    for (int i = 1; i <= V - 1; ++i) {
        changed = false;

#pragma omp parallel for schedule(dynamic, chunkSize) reduction(||:changed)
        for (int j = 0; j < numEdges; ++j) {
            const Edge& edge = graph.edges[j];
            if (dist[edge.u] != INF) {
                int newDist = dist[edge.u] + edge.weight;
                if (newDist < dist[edge.v]) {
#pragma omp critical
                    {
                        if (newDist < dist[edge.v]) {
                            dist[edge.v] = newDist;
                            pred[edge.v] = edge.u;
                            changed = true;
                        }
                    }
                }
            }
        }

        if (!changed) break;
    }

    // Check for negative-weight cycles in parallel
    bool hasNegativeCycle = false;
#pragma omp parallel for
    for (int j = 0; j < numEdges; ++j) {
        const Edge& edge = graph.edges[j];
        if (dist[edge.u] != INF && dist[edge.u] + edge.weight < dist[edge.v]) {
            hasNegativeCycle = true;
#pragma omp cancel for
        }
    }

    auto end = high_resolution_clock::now();
    double duration = duration_cast<milliseconds>(end - start).count() / 1000.0;

    return { dist, pred, hasNegativeCycle, duration };
}

// Optimized Parallel Bellman-Ford Version 2 - Vertex-based with local updates
BFResult bellmanFordParallelV2(const Graph& graph, int source) {
    auto start = high_resolution_clock::now();

    const int V = graph.V;
    vector<int> dist(V, INF);
    vector<int> pred(V, -1);
    dist[source] = 0;

    // Use separate arrays for updates to avoid race conditions
    vector<int> newDist(V, INF);
    vector<int> newPred(V, -1);

    bool changed;
    for (int i = 1; i <= V - 1; ++i) {
        changed = false;
#pragma omp parallel for
        for (int v = 0; v < V; ++v) {
            newDist[v] = dist[v];
            newPred[v] = pred[v];
        }

#pragma omp parallel for reduction(||:changed) schedule(dynamic)
        for (int u = 0; u < V; ++u) {
            if (dist[u] == INF) continue;

            for (const auto& edge : graph.adj[u]) {
                int v = edge.first;
                int weight = edge.second;
                if (dist[u] + weight < newDist[v]) {
                    newDist[v] = dist[u] + weight;
                    newPred[v] = u;
                    changed = true;
                }
            }
        }

        if (!changed) break;

#pragma omp parallel for
        for (int v = 0; v < V; ++v) {
            if (newDist[v] < dist[v]) {
                dist[v] = newDist[v];
                pred[v] = newPred[v];
            }
        }
    }

    // Check for negative-weight cycles
    bool hasNegativeCycle = false;
#pragma omp parallel for
    for (int u = 0; u < V; ++u) {
        if (dist[u] == INF) continue;

        for (const auto& edge : graph.adj[u]) {
            if (dist[u] + edge.second < dist[edge.first]) {
                hasNegativeCycle = true;
#pragma omp cancel for
            }
        }
    }

    auto end = high_resolution_clock::now();
    double duration = duration_cast<milliseconds>(end - start).count() / 1000.0;

    return { dist, pred, hasNegativeCycle, duration };
}

// Optimized Parallel Bellman-Ford Version 3 - Hybrid approach
BFResult bellmanFordParallelV3(const Graph& graph, int source) {
    auto start = high_resolution_clock::now();

    const int V = graph.V;
    vector<int> dist(V, INF);
    vector<int> pred(V, -1);
    dist[source] = 0;

    // Partition the vertices for better load balancing
    const int numThreads = omp_get_max_threads();
    const int chunkSize = max(1, V / (numThreads * 4));

    bool changed;
    for (int i = 1; i <= V - 1; ++i) {
        changed = false;

#pragma omp parallel for schedule(dynamic, chunkSize) reduction(||:changed)
        for (int u = 0; u < V; ++u) {
            if (dist[u] == INF) continue;

            for (const auto& edge : graph.adj[u]) {
                int v = edge.first;
                int weight = edge.second;
                if (dist[u] + weight < dist[v]) {
#pragma omp critical
                    {
                        if (dist[u] + weight < dist[v]) {
                            dist[v] = dist[u] + weight;
                            pred[v] = u;
                            changed = true;
                        }
                    }
                }
            }
        }

        if (!changed) break;
    }

    // Check for negative-weight cycles
    bool hasNegativeCycle = false;
#pragma omp parallel for schedule(dynamic, chunkSize)
    for (int u = 0; u < V; ++u) {
        if (dist[u] == INF) continue;

        for (const auto& edge : graph.adj[u]) {
            if (dist[u] + edge.second < dist[edge.first]) {
                hasNegativeCycle = true;
#pragma omp cancel for
            }
        }
    }

    auto end = high_resolution_clock::now();
    double duration = duration_cast<milliseconds>(end - start).count() / 1000.0;

    return { dist, pred, hasNegativeCycle, duration };
}

bool verifyResults(const BFResult& seqResult, const BFResult& parResult) {
    if (seqResult.hasNegativeCycle != parResult.hasNegativeCycle) {
        return false;
    }

    if (!seqResult.hasNegativeCycle) {
        return seqResult.distances == parResult.distances &&
            seqResult.predecessors == parResult.predecessors;
    }

    return true;
}

void runTests(const vector<int>& sizes, const string& outputFile) {
    ofstream out(outputFile);
    out << "Algorithm,GraphSize,ExecutionTime,Correctness\n";

    // Limit cache size and manage memory
    for (int size : sizes) {
        cout << "\nTesting with graph size: " << size << endl;

        // Generate or retrieve from cache
        Graph graph;
        auto it = graphCache.find(size);
        if (it == graphCache.end()) {
            cout << "Generating graph... ";
            graph = Graph::generateCompleteGraph(size);

            // Manage cache size
            if (graphCache.size() >= MAX_CACHE_SIZE) {
                int toRemove = cacheAccessOrder.front();
                graphCache.erase(toRemove);
                cacheAccessOrder.erase(cacheAccessOrder.begin());
            }
            graphCache[size] = graph;
            cacheAccessOrder.push_back(size);
            cout << "Done" << endl;
        }
        else {
            graph = it->second;
            // Update access order
            cacheAccessOrder.erase(remove(cacheAccessOrder.begin(), cacheAccessOrder.end(), size), cacheAccessOrder.end());
            cacheAccessOrder.push_back(size);
        }

        cout << "Running sequential version... ";
        auto seqResult = bellmanFordSequential(graph, 0);
        cout << "Done (" << seqResult.executionTime << "s)" << endl;

        cout << "Running parallel version 1... ";
        auto v1Result = bellmanFordParallelV1(graph, 0);
        cout << "Done (" << v1Result.executionTime << "s)" << endl;

        cout << "Running parallel version 2... ";
        auto v2Result = bellmanFordParallelV2(graph, 0);
        cout << "Done (" << v2Result.executionTime << "s)" << endl;

        cout << "Running parallel version 3... ";
        auto v3Result = bellmanFordParallelV3(graph, 0);
        cout << "Done (" << v3Result.executionTime << "s)" << endl;

        // Verify correctness
        bool v1Correct = verifyResults(seqResult, v1Result);
        bool v2Correct = verifyResults(seqResult, v2Result);
        bool v3Correct = verifyResults(seqResult, v3Result);

        // Write results to CSV
        out << "Sequential," << size << "," << seqResult.executionTime << "," << "1\n";
        out << "ParallelV1," << size << "," << v1Result.executionTime << "," << v1Correct << "\n";
        out << "ParallelV2," << size << "," << v2Result.executionTime << "," << v2Correct << "\n";
        out << "ParallelV3," << size << "," << v3Result.executionTime << "," << v3Correct << "\n";
    }

    out.close();
    cout << "\nResults saved to " << outputFile << endl;
}

int main() {
    // number of threads
    omp_set_num_threads(omp_get_max_threads());
    cout << "Using " << omp_get_max_threads() << " threads\n";

    vector<int> graphSizes = { 50, 100, 200, 500, 1000, 2000 }; 
    string outputFile = "bellman_ford_results.csv";

    runTests(graphSizes, outputFile);

    return 0;
}