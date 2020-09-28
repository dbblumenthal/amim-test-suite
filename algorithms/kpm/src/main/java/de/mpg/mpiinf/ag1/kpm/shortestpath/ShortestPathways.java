/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package de.mpg.mpiinf.ag1.kpm.shortestpath;

import de.mpg.mpiinf.ag1.kpm.Parameters;
import de.mpg.mpiinf.ag1.kpm.graph.Edge;
import de.mpg.mpiinf.ag1.kpm.graph.LinearPath;
import de.mpg.mpiinf.ag1.kpm.graph.Node;
import de.mpg.mpiinf.ag1.kpm.utils.StatisticsUtility;
import edu.uci.ics.jung.graph.util.Pair;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author nalcaraz
 */
public class ShortestPathways extends HashMap<Pair<Node>, LinearPath> {
    
    Map<Node, Set<Pair<Node>>> nodeCount;
    
    Map<Edge, Set<Pair<Node>>> edgeCount;
    
    public ShortestPathways() {
        super();
        nodeCount = new HashMap<Node, Set<Pair<Node>>>();
        edgeCount = new HashMap<Edge, Set<Pair<Node>>>();        
    }

    public Map<Node, Set<Pair<Node>>> getNodeCount() {
        return nodeCount;
    }

    public Map<Edge, Set<Pair<Node>>> getEdgeCount() {
        return edgeCount;
    }
        
    public void addPathway(LinearPath path) {
        Pair<Node> endPoints = path.getEndNodes();
        if (containsKey(endPoints)) {
            return;
        }
        
        put(endPoints, path);
        
        for (Node node: path.getNodeList()) {            
            if (!nodeCount.containsKey(node)) {
                Set<Pair<Node>> pathSet = new HashSet<Pair<Node>>();
                pathSet.add(endPoints);
                nodeCount.put(node, pathSet);
            } else {
                nodeCount.get(node).add(endPoints);
            }
        }
        
        for (Edge edge: path.getEdgeList()) {
            if (!edgeCount.containsKey(edge)) {
                Set<Pair<Node>> pathSet = new HashSet<Pair<Node>>();
                pathSet.add(endPoints);
                edgeCount.put(edge, pathSet);
            } else {
                edgeCount.get(edge).add(endPoints);
            }

        }
    }
    
    public Collection<LinearPath> getContainingPaths(Node node) {
        Set<Pair<Node>> paths = nodeCount.get(node);
        ArrayList<LinearPath> ret = new ArrayList<LinearPath>(paths.size());
        for (Pair<Node> endPoints: paths) {
            ret.add(get(endPoints));
        }
        return ret;
    }
    
    public Collection<LinearPath> getContainingPaths(Edge edge) {
        Set<Pair<Node>> paths = edgeCount.get(edge);
        ArrayList<LinearPath> ret = new ArrayList<LinearPath>(paths.size());
        for (Pair<Node> endPoints: paths) {
            ret.add(get(endPoints));
        }
        return ret;
    }
    
    public List<Node> getNodesByHits() {        
        List<Node> nodes = new ArrayList<Node>(nodeCount.keySet());
        Collections.sort(nodes, NodeHitComparator);
        return nodes;
    }
    
    public List<Edge> getEdgesByHits() {
        List<Edge> edges = new ArrayList<Edge>(edgeCount.keySet());
        Collections.sort(edges, EdgeHitComparator);
        return edges;
    }
    
    public List<LinearPath> getPathsSortedByLength() {
        ArrayList<LinearPath> paths = new ArrayList<LinearPath>(values());
        Collections.sort(paths, LinearPath.LengthComparator);
        return paths;
    }
    
    public List<LinearPath> getPathsSortedByEdgeWeight() {
        ArrayList<LinearPath> paths = new ArrayList<LinearPath>(values());
        Collections.sort(paths, LinearPath.EdgeWeightComparator);
        return paths;
    }
    
    public List<LinearPath> getPathsSortedByNodeWeight() {
        ArrayList<LinearPath> paths = new ArrayList<LinearPath>(values());
        Collections.sort(paths, LinearPath.NodeWeightComparator);
        return paths;
    }
    
    public List<LinearPath> getPathsSortedBy(Comparator cmp) {
        ArrayList<LinearPath> paths = new ArrayList<LinearPath>(values());
        Collections.sort(paths, cmp);
        return paths;
    }
    
    public void writeShortestPathsFile(String file, Comparator cmp) {
        DecimalFormat df1 = StatisticsUtility.getIntegerFormat(size());
        DecimalFormat df2 = StatisticsUtility.getDoubleFormat(2);
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(checkFileName(file)));
            bw.write("ID" + "\t");
            bw.write("PATH" + "\t");
            bw.write("LENGTH" + "\t");
            bw.write("NODE WEIGHT" + "\t");
            bw.write("EDGE WEIGHT" + "\n");
            int i = 0;
            for (LinearPath path : getPathsSortedBy(cmp)) {
                i++;
                String pathID = df1.format(i);
                String pathNodes = path.getId();
                String pathLength = String.valueOf(path.getLength());
                String nodeWeight = df2.format(path.getNodeWeight());
                String edgeWeight = df2.format(path.getEdgeWeight());
                bw.write(pathID + "\t");
                bw.write(pathNodes + "\t");
                bw.write(pathLength + "\t");
                bw.write(nodeWeight + "\t");
                bw.write(edgeWeight + "\n");
                
            }
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(ShortestPathways.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    public void writeIndividualShortestPaths(String filePrefix, Comparator cmp, int minlength, int maxpaths, Parameters params) {
        DecimalFormat df = StatisticsUtility.getIntegerFormat(size());
        int i = 0;
        Iterator<LinearPath> paths = getPathsSortedBy(cmp).iterator();
        LinearPath path;
        while ((path = paths.next()) != null) {
            i++;
            if (i > maxpaths || path.length < minlength) {
                break;
            }
            BufferedWriter bw = null;            
            try {                
                String pathID = df.format(i);            
                String file = filePrefix + "-" + pathID + "-" + params.SUFFIX + 
                		params.FILE_EXTENSION;
                bw = new BufferedWriter(new FileWriter(checkFileName(file)));
                bw.write("SOURCE" + "\t");
                bw.write("TARGET" + "\n");
                for (Edge edge: path.getEdgeList()) {
                    bw.write(edge.getSource().getId());
                    bw.write("\t");
                    bw.write(edge.getTarget().getId());
                    bw.write("\n");                    
                }
                bw.close();
            } catch (IOException ex) {
                Logger.getLogger(ShortestPathways.class.getName()).log(Level.SEVERE, null, ex);
            } finally {
                try {
                    bw.close();
                } catch (IOException ex) {
                    Logger.getLogger(ShortestPathways.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }

        
    }
    
    public void writeShortestPathsNodeStats(String file, Parameters params) {
        BufferedWriter bw = null;
        double den = (double)size();
        DecimalFormat df = StatisticsUtility.getDoubleFormat(3);
        if (!params.IS_GRAPH_DIRECTED) {
            den *= 2.0;
        }
        try {
            bw = new BufferedWriter(new FileWriter(checkFileName(file)));
            bw.write("NODE" + "\t");
            bw.write("HITS" + "\t");
            bw.write("BETWEENESS" + "\n");
            for (Node node: getNodesByHits()) {
                int hits = nodeCount.get(node).size();
                double between = (double) hits / den;
                String bet = df.format(between);
                bw.write(node.getId() + "\t");
                bw.write(hits + "\t");
                bw.write(bet + "\n");
            }
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(ShortestPathways.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                bw.close();
            } catch (IOException ex) {
                Logger.getLogger(ShortestPathways.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        
    }
    
   public void writeShortestPathsEdgeStats(String file, Parameters params) {
        BufferedWriter bw = null;
        double den = (double)size();
        if (!params.IS_GRAPH_DIRECTED) {
            den *= 2.0;
        }
        DecimalFormat df = StatisticsUtility.getDoubleFormat(3);
        try {
            bw = new BufferedWriter(new FileWriter(checkFileName(file)));
            bw.write("SOURCE" + "\t");
            bw.write("TARGET" + "\t");            
            bw.write("HITS" + "\t");
            bw.write("BETWEENESS" + "\n");

            for (Edge edge: getEdgesByHits()) {
                int hits = edgeCount.get(edge).size();
                double between = (double) hits / den;
                String bet = df.format(between);
                bw.write(edge.getSource().getId() + "\t");
                bw.write(edge.getTarget().getId() + "\t");                
                bw.write(hits + "\t");
                bw.write(bet + "\n");
            }
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(ShortestPathways.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            try {
                bw.close();
            } catch (IOException ex) {
                Logger.getLogger(ShortestPathways.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        
    }
   
   public static String checkFileName(String file) {
       int count = 0;
       while ((new File(file)).exists()) {
           count++;
           
           if (count == 1) {
               int index = file.lastIndexOf(".");
               file = file.substring(0, index) + 
                       "(1)" + file.substring(index, file.length());
           } else {
               int aux = count - 1;
               file = file.replace("(" + aux + ")", "(" + count + ")");
           }
           
       }
       return file;
   }

    public Comparator<Node> getNodeHitComparator() {
        return NodeHitComparator;
    }

    public Comparator<Edge> getEdgeHitComparator() {
        return EdgeHitComparator;
    }
   
   
    
   public Comparator<Node> NodeHitComparator = new Comparator<Node>() {

        @Override
        public int compare(Node o1, Node o2) {
            //Descending
            return (new Integer(nodeCount.get(o2).size()))
                    .compareTo(new Integer(nodeCount.get(o1).size()));
        }
        
    };
    
   public Comparator<Edge> EdgeHitComparator = new Comparator<Edge>() {

        @Override
        public int compare(Edge o1, Edge o2) {
            //Descending
            return (new Integer(edgeCount.get(o2).size()))
                    .compareTo(new Integer(edgeCount.get(o1).size()));
        }
        
    };
    
    
}
