// Digraph.hpp
//
// ICS 46 Spring 2022
// Project #5: Rock and Roll Stops the Traffic
//
// This header file declares a class template called Digraph, which is
// intended to implement a generic directed graph.  The implementation
// uses the adjacency lists technique, so each vertex stores a linked
// list of its outgoing edges.
//
// Along with the Digraph class template is a class DigraphException
// and a couple of utility structs that aren't generally useful outside
// of this header file.
//
// In general, directed graphs are all the same, except in the sense
// that they store different kinds of information about each vertex and
// about each edge; these two types are the type parameters to the
// Digraph class template.

#ifndef DIGRAPH_HPP
#define DIGRAPH_HPP

#include <exception>
#include <functional>
#include <list>
#include <map>
#include <utility>
#include <vector>
#include <iostream>
#include <queue>



// DigraphExceptions are thrown from some of the member functions in the
// Digraph class template, so that exception is declared here, so it
// will be available to any code that includes this header file.

class DigraphException : public std::runtime_error
{
public:
    DigraphException(const std::string& reason);
};


inline DigraphException::DigraphException(const std::string& reason)
    : std::runtime_error{reason}
{
}



// A DigraphEdge lists a "from vertex" (the number of the vertex from which
// the edge points), a "to vertex" (the number of the vertex to which the
// edge points), and an EdgeInfo object.  Because different kinds of Digraphs
// store different kinds of edge information, DigraphEdge is a struct template.

template <typename EdgeInfo>
struct DigraphEdge
{
    int fromVertex;
    int toVertex;
    EdgeInfo einfo;
};



// A DigraphVertex includes two things: a VertexInfo object and a list of
// its outgoing edges.  Because different kinds of Digraphs store different
// kinds of vertex and edge information, DigraphVertex is a struct template.

template <typename VertexInfo, typename EdgeInfo>
struct DigraphVertex
{
    VertexInfo vinfo;
    std::list<DigraphEdge<EdgeInfo>> edges;
};



// Digraph is a class template that represents a directed graph implemented
// using adjacency lists.  It takes two type parameters:
//
// * VertexInfo, which specifies the kind of object stored for each vertex
// * EdgeInfo, which specifies the kind of object stored for each edge
//
// You'll need to implement the member functions declared here; each has a
// comment detailing how it is intended to work.
//
// Each vertex in a Digraph is identified uniquely by a "vertex number".
// Vertex numbers are not necessarily sequential and they are not necessarily
// zero- or one-based.

template <typename VertexInfo, typename EdgeInfo>
class Digraph
{
public:
    // The default constructor initializes a new, empty Digraph so that
    // contains no vertices and no edges.
    Digraph();

    // The copy constructor initializes a new Digraph to be a deep copy
    // of another one (i.e., any change to the copy will not affect the
    // original).
    Digraph(const Digraph& d);

    // The move constructor initializes a new Digraph from an expiring one.
    Digraph(Digraph&& d) noexcept;

    // The destructor deallocates any memory associated with the Digraph.
    ~Digraph() noexcept;

    // The assignment operator assigns the contents of the given Digraph
    // into "this" Digraph, with "this" Digraph becoming a separate, deep
    // copy of the contents of the given one (i.e., any change made to
    // "this" Digraph afterward will not affect the other).
    Digraph& operator=(const Digraph& d);

    // The move assignment operator assigns the contents of an expiring
    // Digraph into "this" Digraph.
    Digraph& operator=(Digraph&& d) noexcept;

    // vertices() returns a std::vector containing the vertex numbers of
    // every vertex in this Digraph.
    std::vector<int> vertices() const;

    // edges() returns a std::vector of std::pairs, in which each pair
    // contains the "from" and "to" vertex numbers of an edge in this
    // Digraph.  All edges are included in the std::vector.
    std::vector<std::pair<int, int>> edges() const;

    // This overload of edges() returns a std::vector of std::pairs, in
    // which each pair contains the "from" and "to" vertex numbers of an
    // edge in this Digraph.  Only edges outgoing from the given vertex
    // number are included in the std::vector.  If the given vertex does
    // not exist, a DigraphException is thrown instead.
    std::vector<std::pair<int, int>> edges(int vertex) const;

    // vertexInfo() returns the VertexInfo object belonging to the vertex
    // with the given vertex number.  If that vertex does not exist, a
    // DigraphException is thrown instead.
    VertexInfo vertexInfo(int vertex) const;

    // edgeInfo() returns the EdgeInfo object belonging to the edge
    // with the given "from" and "to" vertex numbers.  If either of those
    // vertices does not exist *or* if the edge does not exist, a
    // DigraphException is thrown instead.
    EdgeInfo edgeInfo(int fromVertex, int toVertex) const;

    // addVertex() adds a vertex to the Digraph with the given vertex
    // number and VertexInfo object.  If there is already a vertex in
    // the graph with the given vertex number, a DigraphException is
    // thrown instead.
    void addVertex(int vertex, const VertexInfo& vinfo);

    // addEdge() adds an edge to the Digraph pointing from the given
    // "from" vertex number to the given "to" vertex number, and
    // associates with the given EdgeInfo object with it.  If one
    // of the vertices does not exist *or* if the same edge is already
    // present in the graph, a DigraphException is thrown instead.
    void addEdge(int fromVertex, int toVertex, const EdgeInfo& einfo);

    // removeVertex() removes the vertex (and all of its incoming
    // and outgoing edges) with the given vertex number from the
    // Digraph.  If the vertex does not exist already, a DigraphException
    // is thrown instead.
    void removeVertex(int vertex);

    // removeEdge() removes the edge pointing from the given "from"
    // vertex number to the given "to" vertex number from the Digraph.
    // If either of these vertices does not exist *or* if the edge
    // is not already present in the graph, a DigraphException is
    // thrown instead.
    void removeEdge(int fromVertex, int toVertex);

    // vertexCount() returns the number of vertices in the graph.
    int vertexCount() const noexcept;

    // edgeCount() returns the total number of edges in the graph,
    // counting edges outgoing from all vertices.
    int edgeCount() const noexcept;

    // This overload of edgeCount() returns the number of edges in
    // the graph that are outgoing from the given vertex number.
    // If the given vertex does not exist, a DigraphException is
    // thrown instead.
    int edgeCount(int vertex) const;

    // isStronglyConnected() returns true if the Digraph is strongly
    // connected (i.e., every vertex is reachable from every other),
    // false otherwise.
    bool isStronglyConnected() const;

    // findShortestPaths() takes a start vertex number and a function
    // that takes an EdgeInfo object and determines an edge weight.
    // It uses Dijkstra's Shortest Path Algorithm to determine the
    // shortest paths from the start vertex to every other vertex
    // in the graph.  The result is returned as a std::map<int, int>
    // where the keys are vertex numbers and the value associated
    // with each key k is the predecessor of that vertex chosen by
    // the algorithm.  For any vertex without a predecessor (e.g.,
    // a vertex that was never reached, or the start vertex itself),
    // the value is simply a copy of the key.
    std::map<int, int> findShortestPaths(
        int startVertex,
        std::function<double(const EdgeInfo&)> edgeWeightFunc) const;


private:
    // Add whatever member variables you think you need here.  One
    // possibility is a std::map where the keys are vertex numbers
    // and the values are DigraphVertex<VertexInfo, EdgeInfo> objects.
    std::map<int, DigraphVertex<VertexInfo, EdgeInfo>> _vertices;
    int _vCount;
    int _eCount;
    


    // You can also feel free to add any additional member functions
    // you'd like (public or private), so long as you don't remove or
    // change the signatures of the ones that already exist.

    void getEdges(int vertex, std::vector<std::pair<int, int>>& edges) const;

    bool containsVertex(int vertex) const;

    int DFTr(int vertex, std::map<int, bool>& visited) const;

};




// You'll need to implement the member functions below.  There's enough
// code in place to make them compile, but they'll all need to do the
// correct thing instead.

template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>::Digraph()
    : _vCount{0}, _eCount{0}
{
}


template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>::Digraph(const Digraph& d)
    : _vertices{d._vertices}, _vCount{d._vCount}, _eCount{d._eCount}
{
}


template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>::Digraph(Digraph&& d) noexcept
{
    std::swap(_vertices, d._vertices);
    std::swap(_vCount, d._vCount);
    std::swap(_eCount, d._eCount);
}


template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>::~Digraph() noexcept
{
}


template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>& Digraph<VertexInfo, EdgeInfo>::operator=(const Digraph& d)
{
    _vertices = d._vertices;
    _vCount = d._vCount;
    _eCount = d._eCount;
    return *this;
}


template <typename VertexInfo, typename EdgeInfo>
Digraph<VertexInfo, EdgeInfo>& Digraph<VertexInfo, EdgeInfo>::operator=(Digraph&& d) noexcept
{
    std::swap(_vertices, d._vertices);
    std::swap(_vCount, d._vCount);
    std::swap(_eCount, d._eCount);
    return *this;
}


template <typename VertexInfo, typename EdgeInfo>
std::vector<int> Digraph<VertexInfo, EdgeInfo>::vertices() const
{
    std::vector<int> v;
    for (auto it = _vertices.begin(); it != _vertices.end(); ++it)
    {
        v.push_back(it->first);
    }
    return v;
}


template <typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo, EdgeInfo>::getEdges(int vertex, 
    std::vector<std::pair<int, int>>& edges) const
{
    for (DigraphEdge e : _vertices.at(vertex).edges)
    {
        edges.push_back(std::pair<int, int>(e.fromVertex, e.toVertex));
    }
}


template <typename VertexInfo, typename EdgeInfo>
std::vector<std::pair<int, int>> Digraph<VertexInfo, EdgeInfo>::edges() const
{
    std::vector<std::pair<int, int>> edges;
    for (auto it = _vertices.begin(); it != _vertices.end(); ++it)
    {
        getEdges(it->first, edges); 
    }
    return edges;
}


template <typename VertexInfo, typename EdgeInfo>
bool Digraph<VertexInfo, EdgeInfo>::containsVertex(int vertex) const
{
    if (_vertices.count(vertex) == 0)
    {
        return false;
    }
    return true;
}


template <typename VertexInfo, typename EdgeInfo>
std::vector<std::pair<int, int>> Digraph<VertexInfo, EdgeInfo>::edges(int vertex) const
{
    if (!containsVertex(vertex))
    {
        throw DigraphException("Does not contain vertex");
    }
    std::vector<std::pair<int, int>> edges;
    getEdges(vertex, edges); 
    return edges;
}


template <typename VertexInfo, typename EdgeInfo>
VertexInfo Digraph<VertexInfo, EdgeInfo>::vertexInfo(int vertex) const
{
    if (!containsVertex(vertex))
    {
        throw DigraphException("Does not contain vertex");
    }
    return _vertices.at(vertex).vinfo;
}


template <typename VertexInfo, typename EdgeInfo>
EdgeInfo Digraph<VertexInfo, EdgeInfo>::edgeInfo(int fromVertex, int toVertex) const
{
    if (!containsVertex(fromVertex) || !containsVertex(toVertex))
    {
        throw DigraphException("Does not contain vertex");
    }
    for (DigraphEdge e : _vertices.at(fromVertex).edges)
    {
        if (e.toVertex == toVertex)
        {
            return e.einfo;
        }
    }
    throw DigraphException("Graph doesn't contain edge");
}


template <typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo, EdgeInfo>::addVertex(int vertex, const VertexInfo& vinfo)
{
    if (containsVertex(vertex))
    {
        throw DigraphException("already contains vertex");
    }
    _vertices[vertex] = DigraphVertex<VertexInfo, EdgeInfo>{vinfo};
    ++_vCount;
}


template <typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo, EdgeInfo>::addEdge(int fromVertex, int toVertex, const EdgeInfo& einfo)
{
    if (!containsVertex(fromVertex) || !containsVertex(toVertex))
    {
        throw DigraphException("Does not contain vertex");
    }
    for (DigraphEdge e :  _vertices[fromVertex].edges)
    {
        if (e.toVertex == toVertex)
        {
            throw DigraphException("Graph already contains edge");
        }
    }
    _vertices[fromVertex].edges.push_back(DigraphEdge<EdgeInfo>{fromVertex, toVertex, einfo});
    ++_eCount;
}


template <typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo, EdgeInfo>::removeVertex(int vertex)
{
    if (!containsVertex(vertex))
    {
        throw DigraphException("Graph does not contain vertex");
    }
    _vertices.erase(vertex);

    for (auto it = _vertices.begin(); it != _vertices.end(); ++it)
    {
        it->second.edges.remove_if([vertex](DigraphEdge<EdgeInfo> e){return e.toVertex == vertex;});
    }
    --_vCount;
}


template <typename VertexInfo, typename EdgeInfo>
void Digraph<VertexInfo, EdgeInfo>::removeEdge(int fromVertex, int toVertex)
{
    if (!containsVertex(fromVertex) || !containsVertex(toVertex))
    {
        throw DigraphException("Does not contain vertex");
    }
    auto& edges = _vertices[fromVertex].edges;  
    int beforeSize = edges.size();

    edges.remove_if([toVertex](DigraphEdge<EdgeInfo> e){return e.toVertex == toVertex;});
    if (beforeSize == edges.size())
    {
        throw DigraphException("Does not contain edge");
    }
    --_eCount;
}


template <typename VertexInfo, typename EdgeInfo>
int Digraph<VertexInfo, EdgeInfo>::vertexCount() const noexcept
{
    return _vCount;
}


template <typename VertexInfo, typename EdgeInfo>
int Digraph<VertexInfo, EdgeInfo>::edgeCount() const noexcept
{
    return _eCount;  
}


template <typename VertexInfo, typename EdgeInfo>
int Digraph<VertexInfo, EdgeInfo>::edgeCount(int vertex) const
{
    return _vertices.at(vertex).edges.size(); 
}


template <typename VertexInfo, typename EdgeInfo>
int Digraph<VertexInfo, EdgeInfo>::DFTr(int vertex, std::map<int, bool>& visited) const
{
    int numVisited = 1;
    visited[vertex] = true; 

    for (DigraphEdge e : _vertices.at(vertex).edges)
    {
        if (!visited[e.toVertex]) 
        {
            numVisited += DFTr(e.toVertex, visited);
        }
    }
    return numVisited;
}


template <typename VertexInfo, typename EdgeInfo>
bool Digraph<VertexInfo, EdgeInfo>::isStronglyConnected() const
{
    std::map<int, bool> visited;

    for (int vertex : vertices())
    {
        for (int vertex : vertices())
        {
            visited[vertex] = false; 
        }
        int numVisited = DFTr(vertex, visited);
        if (numVisited != _vCount)
        {
            return false;
        }
    }
    return true;
}


template <typename VertexInfo, typename EdgeInfo>
std::map<int, int> Digraph<VertexInfo, EdgeInfo>::findShortestPaths(
    int startVertex,
    std::function<double(const EdgeInfo&)> edgeWeightFunc) const
{
    std::map<int, bool> k;
    std::map<int, int> p;
    std::map<int, double> d;

    for (int v :  vertices())
    {
        k[v] = false; 
        p[v] = v; 
        d[v] = std::numeric_limits<double>::infinity();
    }
    d[startVertex] = 0;

    std::map<int, int> pr;
    pr[startVertex] = 0;

    auto cmp = [&](int left, int right){return pr[left] > pr[right];};
    std::priority_queue<int, std::vector<int>, decltype(cmp)> pq(cmp); 
    pq.push(startVertex);

    while (!pq.empty())
    {
        int v = pq.top(); 
        pq.pop();

        if (k[v] == false)
        {
            k[v] = true;

            for (DigraphEdge e : _vertices.at(v).edges)
            {
                int w = e.toVertex; 

                if (d[w] > d[v] + edgeWeightFunc(e.einfo))
                {
                    d[w] = d[v] + edgeWeightFunc(e.einfo); 
                    p[w] = v;

                    pr[w] = d[w];
                    pq.push(w);
                }
            }
        }
    }
    return p;
}



#endif

