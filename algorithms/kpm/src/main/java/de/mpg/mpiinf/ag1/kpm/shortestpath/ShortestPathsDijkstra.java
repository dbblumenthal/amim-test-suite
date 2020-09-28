package de.mpg.mpiinf.ag1.kpm.shortestpath;

import de.mpg.mpiinf.ag1.kpm.graph.Edge;
import de.mpg.mpiinf.ag1.kpm.graph.LinearPath;
import de.mpg.mpiinf.ag1.kpm.graph.Node;
import edu.uci.ics.jung.algorithms.shortestpath.DijkstraShortestPath;
import edu.uci.ics.jung.graph.Graph;

import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author nalcaraz
 */
public class ShortestPathsDijkstra {
    
    Graph<Node, Edge> g;
    
    boolean isDirected;
    
    List<Node> nodeList;
    
    ShortestPathways shortestPathways;
    
    public ShortestPathsDijkstra(Graph<Node, Edge> g, boolean isDirected) {
        this.g = g;        
        this.isDirected = isDirected;
        nodeList = new ArrayList<Node>(g.getVertices());
        shortestPathways = new ShortestPathways();
        int n = nodeList.size();
        DijkstraShortestPath<Node,Edge> dsp = new DijkstraShortestPath<Node, Edge>(g);
        if (!isDirected) {
            for (int i = 0; i < n - 1; i++) {
                for (int j = i + 1; j < n; j++) {
                    Node source = nodeList.get(i);

                    Node target = nodeList.get(j);

                    List<Edge> path = dsp.getPath(source, target);
                    if (path.isEmpty()) {
                        //System.out.println("EMPTY !");
                    } else {
                        shortestPathways.addPathway(new LinearPath(path));

                    }

                }
            }
        } else {
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        Node source = nodeList.get(i);

                        Node target = nodeList.get(j);

                        List<Edge> path = dsp.getPath(source, target);
                        if (path.isEmpty()) {
                            //System.out.println("EMPTY !");
                        } else {
                            shortestPathways.addPathway(new LinearPath(path));
                        }
                    }

                }
            }
        }
    }
    
    public ShortestPathways getShortestPathways() {
        return shortestPathways;
    }
    
}
