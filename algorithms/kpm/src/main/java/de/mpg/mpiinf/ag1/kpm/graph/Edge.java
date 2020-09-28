package de.mpg.mpiinf.ag1.kpm.graph;

/**
 *
 * @author nalcaraz
 */
public class Edge implements Comparable {
    
    private String id;
    
    private double weight;
    
    private Node source;
    
    private Node target;
    
    public Edge(String id, Node source, Node target) {
        this.id = id;
        weight = 1.0;
        this.source = source;
        this.target = target;
    }
    
    public Edge(String id, double weight, Node source, Node target) {
        this.id = id;
        this.weight = weight;
        this.source = source;
        this.target = target;
    }
    
    public String getId() {
        return id;
    }
    
    public double getWeight() {
        return weight;
    }
    
    public void setWeight() {
        this.weight = weight;
    }
    
    public Node getSource() {
        return source;
    }
    
    public Node getTarget() {
        return target;
    }

    @Override
    public int compareTo(Object o) {
        return (new Double(weight)).compareTo(new Double(((Edge)o).getWeight()));
    }
    
    @Override
    public boolean equals(Object o) {
        return id.equals(((Edge)o).getId());
    }
    
    @Override
    public String toString() {
        return id;
    }
     
}
