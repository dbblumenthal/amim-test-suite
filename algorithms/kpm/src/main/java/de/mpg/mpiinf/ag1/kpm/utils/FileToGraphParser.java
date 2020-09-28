/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package de.mpg.mpiinf.ag1.kpm.utils;

import de.mpg.mpiinf.ag1.kpm.graph.Edge;
import de.mpg.mpiinf.ag1.kpm.graph.Node;
import edu.uci.ics.jung.graph.DirectedSparseGraph;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author nalcaraz
 */
public class FileToGraphParser {
    
    public static Graph<Node, Edge> parseSif(String file, boolean isDirected) {
        
        Graph g = null;

        if (isDirected) {
            g = new DirectedSparseGraph<Node, Edge>();            
        } else {
            g = new UndirectedSparseGraph<Node, Edge>();
        }
        BufferedReader br = null;        
        try {
            br = new BufferedReader(new FileReader(file));
            String line = "";
            while ((line = br.readLine()) != null) {
                String[] fields = line.split("\t");
                String id1 = fields[0].trim();
                String id2 = fields[2].trim();
                
                Node node1 = new Node(id1);
                Node node2 = new Node(id2);
                Edge edge = new Edge(id1 + "-" + id2, node1, node2);
                g.addEdge(edge, node1, node2);
                
            }
            br.close();
        } catch (FileNotFoundException ex) {
            Logger.getLogger(FileToGraphParser.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ioe) {
            Logger.getLogger(FileToGraphParser.class.getName()).log(Level.SEVERE, null, ioe);
        } 
        
        return g;
    }
    
    public static Graph<Node, Edge> parseTable(String file, boolean isDirected, 
            boolean hasNodeWeights, boolean hasEdgeWeights) {
        Graph g = null;
        Map<String, Node> addedNodes = new HashMap<String, Node>();
        Map<String, Edge> addedEdges = new HashMap<String, Edge>();
        
        
        if (isDirected) {
            g = new DirectedSparseGraph<Node, Edge>();
        } else {
            g = new UndirectedSparseGraph<Node, Edge>();
        }
        
                BufferedReader br = null;        
        try {
            br = new BufferedReader(new FileReader(file));
            String line = "";
            while ((line = br.readLine()) != null) {
                String[] fields = line.split("\t");
                String id1 = fields[0].trim();
                String id2 = fields[1].trim();
                double nodew1 = 1.0;
                double nodew2 = 1.0;
                double edgew = 1.0;
                int eIndex = 2;
                if (hasNodeWeights) {
                    nodew1 = Double.parseDouble(fields[2].trim());
                    nodew2 = Double.parseDouble(fields[3].trim());
                    eIndex = 4;
                }
                if (hasEdgeWeights) {
                    edgew = Double.parseDouble(fields[eIndex].trim());
                }
                Node node1 = new Node(id1, nodew1);
                if (addedNodes.containsKey(id1)) {
                    node1 = addedNodes.get(id1);
                } else {
                    addedNodes.put(id1, node1);
                }
                Node node2 = new Node(id2, nodew2);
                if (addedNodes.containsKey(id2)) {
                    node2 = addedNodes.get(id2);
                } else {
                    addedNodes.put(id2, node2);
                }                
                String edgeId = id1 + "-" + id2;
                Edge edge = new Edge(edgeId, edgew, node1, node2);
                if (addedEdges.containsKey(edgeId)) {
                    edge = addedEdges.get(edgeId);
                } else {
                    addedEdges.put(edgeId, edge);
                }
                g.addEdge(edge, node1, node2);
                
            }
            br.close();
        } catch (FileNotFoundException ex) {
            Logger.getLogger(FileToGraphParser.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ioe) {
            Logger.getLogger(FileToGraphParser.class.getName()).log(Level.SEVERE, null, ioe);
        } 

        
        return g;
        
    }
    
    public static Graph<Node, Edge> parseTable(String file, boolean isDirected) {

        boolean hasNodeWeights = false;
        boolean hasEdgeWeights = false;

        BufferedReader br = null;
        try {
            br = new BufferedReader(new FileReader(file));
            String line = br.readLine();

            String[] fields = line.split("\t");
            int cols = fields.length;
            if (cols == 5) {
                hasNodeWeights = true;
                hasEdgeWeights = true;
            } else if (cols == 3) {
                hasEdgeWeights = true;
            } else if (cols == 4) {
                hasNodeWeights = true;
            }

            
            br.close();
        } catch (FileNotFoundException ex) {
            Logger.getLogger(FileToGraphParser.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ioe) {
            Logger.getLogger(FileToGraphParser.class.getName()).log(Level.SEVERE, null, ioe);
        }
        return parseTable(file, isDirected, hasNodeWeights, hasEdgeWeights);
    }
    
}
