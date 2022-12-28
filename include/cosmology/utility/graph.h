#pragma once

#include <vector>

namespace graph {
    enum VertexColor {
        Gray, Black
    };

    template<class Vertex, class Graph, class Visitor>
    void DepthFirstSearch(const Graph &graph, Vertex origin, Visitor &visitor) {
        std::vector<Vertex> vertices;
        std::unordered_map<Vertex, VertexColor> color;
        vertices.push_back(origin);

        while (!vertices.empty()) {
            auto vertex = vertices.back();

            if (!color.count(vertex)) {
                color[vertex] = VertexColor::Gray;
                visitor.discoverVertex(vertex);
                for (Vertex neighbor: getVertexNeighbors(graph, vertex)) {
                    vertices.push_back(neighbor);
                }
            } else {
                vertices.pop_back();
                if (color[vertex] == VertexColor::Gray) {
                    color[vertex] = VertexColor::Black;
                    visitor.leaveVertex(vertex);
                }
            }
        }
    }

// graph's vertices are 0,1,...n-1, where n is number of vertices
    struct Graph {
        std::vector<std::vector<int>> edges;
    };

    const std::vector<int> &getVertexNeighbors(const Graph &graph, int vertex) {
        return graph.edges[vertex];
    }

    Graph transposeGraph(const Graph &graph) {
        Graph transpose;
        transpose.edges.resize(graph.edges.size());

        for (int from = 0; from < graph.edges.size(); ++from) {
            for (int to: graph.edges[from]) {
                transpose.edges[to].push_back(from);
            }
        }

        return transpose;
    }

    std::pair<Graph, std::unordered_map<int, int>> buildSubgraph(const Graph &graph, const std::vector<int> &vertices) {
        Graph subgraph;
        subgraph.edges.resize(vertices.size());
        auto vertexToIndex = createVectorToIndex(vertices);

        for (int from = 0; from < graph.edges.size(); ++from) {
            if (vertexToIndex.count(from)) {
                int fromIndex = vertexToIndex.at(from);

                for (int to: graph.edges[from]) {
                    if (vertexToIndex.count(to)) {
                        subgraph.edges[fromIndex].push_back(vertexToIndex.at(to));
                    }
                }
            }
        }

        return {subgraph, vertexToIndex};
    }

    class TopologicalVisitor {
    public:
        void discoverVertex(int vertex) {
            if (visited_.count(vertex)) {
                throw std::logic_error("graph contains cycles - topological order is impossible to define");
            }

            visited_.insert(vertex);
        }

        void leaveVertex(int vertex) {
            order_.push_back(vertex);
        }

        bool isVisited(int vertex) const {
            return visited_.count(vertex);
        }

        std::vector<int> getOrder() {
            return std::move(order_);
        }

    private:
        std::unordered_set<int> visited_;
        std::vector<int> order_;
    };

    std::vector<int> getTopologicalOrder(const Graph &graph) {
        TopologicalVisitor visitor;

        for (int vertex = 0; vertex < graph.edges.size(); ++vertex) {
            if (!visitor.isVisited(vertex)) {
                DepthFirstSearch(graph, vertex, visitor);
            }
        }

        auto order = visitor.getOrder();
        std::reverse(order.begin(), order.end());

        return order;
    }
} // namespace graph