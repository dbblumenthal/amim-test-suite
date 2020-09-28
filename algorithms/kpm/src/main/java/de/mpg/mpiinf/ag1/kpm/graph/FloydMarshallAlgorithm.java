/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package de.mpg.mpiinf.ag1.kpm.graph;

import de.mpg.mpiinf.ag1.kpm.shortestpath.ShortestPathways;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.util.Pair;

import java.util.ArrayList;
import java.util.Stack;

/**
 *
 * @author nalcaraz
 */
public class FloydMarshallAlgorithm {

    public Graph<Node, Edge> g;
 
    private double[][] dist;
 
    private Edge[][] edgeTo;
 
    private int n;
 
    private int m;
 
    ArrayList<Node> nodeList;
 
    ShortestPathways shortestPaths;
 
    private boolean isDirected;

    public FloydMarshallAlgorithm(Graph<Node, Edge> g, boolean isDirected) {
        this.g = g;
        this.isDirected = isDirected;
        n = g.getVertexCount();
        m = g.getEdgeCount();
        dist = new double[n][n];
        edgeTo = new Edge[n][n];
        nodeList = new ArrayList<Node>(g.getVertices());
        shortestPaths = new ShortestPathways();
        
        // initialize distances to infinity
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                dist[i][j] = Double.POSITIVE_INFINITY;
                edgeTo[i][j] = null;
            }
        }
        // initialize distances using edge-weighted digraph's
        for (int i = 0; i < n; i++) {
            Node v = nodeList.get(i);

            for (Edge edge : g.getIncidentEdges(v)) {
                int s = nodeList.indexOf(edge.getSource());
                int t = nodeList.indexOf(edge.getTarget());
                dist[s][t] = edge.getWeight();
                edgeTo[s][t] = edge;
                if (!isDirected) {
                    dist[t][s] = edge.getWeight();
                    edgeTo[t][s] = edge;
                }
            }
            if (dist[i][i] > 0.0) {
                dist[i][i] = 0.0;
                edgeTo[i][i] = null;
            }

            
   
        }
        
         // Floyd-Warshall updates
        for (int i = 0; i < n; i++) {
            // compute shortest paths using only 0, 1, ..., i as intermediate vertices
            for (int v = 0; v < n; v++) {
                if (edgeTo[v][i] == null) {
                    continue;
                }  // optimization
                for (int w = 0; w < n; w++) {
                    if (dist[v][w] > dist[v][i] + dist[i][w]) {
                        dist[v][w] = dist[v][i] + dist[i][w];
                        edgeTo[v][w] = edgeTo[i][w];
                    }
                }
            }
        }
        
//        for (int i = 0; i < n; i++) {
//            for (int j = 0; j < n; j++) {
//                if (edgeTo[i][j] != null) {
//                    System.out.println(edgeTo[i][j]);
//                }
//            }
//        }
 
        fillShortestMap();

    }
    
    public boolean hasPath(int s, int t) {
        return dist[s][t] < Double.POSITIVE_INFINITY;
    }
    
    public boolean hasPath(Node s, Node t) {
        int indexS = nodeList.indexOf(s);
        int indexT = nodeList.indexOf(t);
        return hasPath(indexS, indexT);
    }
    
    public double dist(int s, int t) {
        return dist[s][t];
    }
    
    public double dist(Node s, Node t) {
        int indexS = nodeList.indexOf(s);
        int indexT = nodeList.indexOf(t);
        return dist[indexS][indexT];
    }

    

    private Stack<Edge> getShortestPath(int s, int t) {
        if (!hasPath(s,t)) {
            return null;
        }
        Stack<Edge> path = new Stack<Edge>();
        Edge e;
        int start = s;
        int next = t;
        Node v = nodeList.get(s);
        Node w = nodeList.get(t);
        System.out.println("START NODE: " + v.getId());
        System.out.println("END NODE: " + w.getId());

        while ((e = edgeTo[s][next]) != null) {           
            if (!path.contains(e)) {
                path.push(e);
                System.out.println("FROM: " + e.getSource().getId());
                System.out.println("TO: " + e.getTarget().getId());
 
            } else {
                break;
            }
            next = nodeList.indexOf(e.getSource());
        }
        return path;
        
    }

    private void fillShortestMap() {
        for (int i = 0; i < nodeList.size(); i++) {
            for (int j = 0; j < nodeList.size(); j++) {
                if ((i != j) && hasPath(i,j)) {
                    Pair<Node> ends = new Pair(nodeList.get(i),nodeList.get(j));
                    Stack<Edge> path = getShortestPath(i, j);
                    System.out.println("START SHORTEST PATH");
                    shortestPaths.put(ends, new LinearPath(path));
                }
            }
        }
    }
    
    
//    private void fillShortestMapUndirected() {
//        for (int i = 0; i < nodeList.size() - 1; i++) {
//            for (int j = i + 1; j < nodeList.size(); j++) {
//                if (dist[i][j] < Integer.MAX_VALUE) {
//                    Pair<Integer> endsI= new Pair(i,j);
//                    Pair<Node> ends = new Pair(nodeList.get(i),nodeList.get(j));
//                    Stack<Edge> path = getShortestPath(i, j);
//                    shortestPaths.put(ends, new LinearPath(path));
//                }
//            }
//        }
//    }
    
    public ShortestPathways getShortestPaths() {
        return shortestPaths;
    }
    
    
}
