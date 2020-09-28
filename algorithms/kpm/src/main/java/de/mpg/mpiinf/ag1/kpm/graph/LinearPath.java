/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package de.mpg.mpiinf.ag1.kpm.graph;

import edu.uci.ics.jung.graph.util.Pair;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Stack;

/**
 *
 * @author nalcaraz
 */
public class LinearPath {
    
    List<Node> nodeList;
    
    List<Edge> edgeList;
    
    String id;
    
    double nodeWeight;

    double edgeWeight;
    
    public int length;
        
    public LinearPath() {
        id = "";
        nodeList = new ArrayList<Node>();        
        edgeList = new ArrayList<Edge>();
        nodeWeight = 0.0;
        edgeWeight = 0.0;
        length = 0;
    }
    
    public LinearPath(int length) {
        id = "";
        nodeList = new ArrayList<Node>(length + 1);        
        edgeList = new ArrayList<Edge>(length);
        nodeWeight = 0.0;
        edgeWeight = 0.0;
        this.length = length;

    }
    
    public LinearPath(List<Edge> path) {
        Edge first = path.get(0);
        Node firstS = first.getSource();
        Node firstT = first.getTarget();
        id = first.getSource() + "-" + first.getTarget();
        length = path.size();
        nodeList = new ArrayList<Node>(length + 1);
        nodeList.add(firstS);
        nodeList.add(firstT);
        edgeList = new ArrayList<Edge>(length);      
        edgeList.add(first);
        nodeWeight = firstS.getWeight() + firstT.getWeight();
        edgeWeight = first.getWeight();
        for (int i = 1; i < length; i++) {
            Edge edge = path.get(i);
            Node source = edge.getSource();
            Node target = edge.getTarget();
            if (!nodeList.contains(source)) {
                nodeList.add(source);
                id += "-" + source.getId();
                nodeWeight += source.getWeight();
            } else if (!nodeList.contains(target)){
                nodeList.add(target);
                id += "-" + target.getId();
                nodeWeight += target.getWeight();                
            }
            edgeList.add(edge);
            edgeWeight += edge.getWeight();
        }
    }
    
    public String getId() {
        return id;
    }

    public List<Node> getNodeList() {
        return nodeList;
    }

    public List<Edge> getEdgeList() {
        return edgeList;
    }
    
    
    public double getNodeWeight() {
        return nodeWeight;
    }
    
    public double getEdgeWeight() {
        return edgeWeight;
    }
        
    public int getLength() {
        return length;
    }
    
    public boolean addNode(Node node) {
        if (nodeList.contains(node)) {
            return false;
        }
        nodeList.add(node);
        nodeWeight += node.getWeight();
        if (id.equals("")) {
            id = node.getId();
        } else {
            id += "-" + node.getId();
        }
        return true;
    }
    
    public boolean addEdge(Edge edge) {
        if (edgeList.contains(edge)) {
            return false;
        }
        edgeList.add(edge);
        edgeWeight += edge.getWeight();
        length++;
        return true;
    }
    
    public Node getSourceNode() {
        return nodeList.get(0);
    }
    
    public Node getTargetNode() {
        return nodeList.get(nodeList.size() - 1);
    }
    
    public Pair<Node> getEndNodes() {
        return new Pair<Node>(getSourceNode(), getTargetNode());
    }
    
    public Edge getSourceEdge() {
        return edgeList.get(0);
    }
    
    public Edge getTargetEdge() {
        return edgeList.get(edgeList.size() - 1);
    }
    
    public Pair<Edge> getEndEdges() {
        return new Pair<Edge>(getSourceEdge(), getTargetEdge());
    }
    
    
    @Override
    public boolean equals(Object obj) {
        if (obj == null) {
            return false;
        }
        if (! (obj instanceof LinearPath)) {
            return false;        
        }
        return id.equals(((LinearPath)obj).getId());
    }

    @Override
    public int hashCode() {
        int hash = 7;
        hash = 29 * hash + (this.nodeList != null ? this.nodeList.hashCode() : 0);
        return hash;
    }
    
    public static Comparator<LinearPath> LengthComparator 
            = new Comparator<LinearPath>() {

        @Override
        public int compare(LinearPath o1, LinearPath o2) {
            Integer length1 = new Integer(o1.getLength());
            Integer length2 = new Integer(o2.getLength());
            
            // Descending
            return length2.compareTo(length1);
        }
        
    };
    
    public static Comparator<LinearPath> NodeWeightComparator 
            = new Comparator<LinearPath>() {

        @Override
        public int compare(LinearPath o1, LinearPath o2) {
            Double weight1 = new Double(o1.getNodeWeight());
            Double weight2 = new Double(o2.getNodeWeight());
            
            // Descending
            return weight2.compareTo(weight1);
        }
        
    };
    
    public static Comparator<LinearPath> EdgeWeightComparator 
            = new Comparator<LinearPath>() {

        @Override
        public int compare(LinearPath o1, LinearPath o2) {
            Double weight1 = new Double(o1.getEdgeWeight());
            Double weight2 = new Double(o2.getEdgeWeight());
            
            // Descending
            return weight2.compareTo(weight1);
        }
        
    };
    
    public boolean isSubPathway(LinearPath other) {
        if (other.getLength() > this.getLength()) {
            return false;
        }
        
        int last = this.getEdgeList().indexOf(other.getEdgeList().get(0));
        if (last == -1) {
            return false;
        }
        for (int i = 1; i < other.getEdgeList().size(); i++) {
            Edge edge = other.getEdgeList().get(i);
            int index = this.getEdgeList().indexOf(edge);
            if (index == -1 || index != last + 1) {
                return false;
            } else {
                last = index;
            }
        } 
        return true;        
    }
    
    public boolean isSuperPathway(LinearPath other) {
        return other.isSubPathway(this);
    }
    
}
