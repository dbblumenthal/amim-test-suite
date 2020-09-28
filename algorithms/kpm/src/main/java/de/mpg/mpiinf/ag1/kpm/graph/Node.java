package de.mpg.mpiinf.ag1.kpm.graph;

/**
 *
 * @author nalcaraz
 */
public class Node implements Comparable {
    
    private String id;
    
    private double weight;
    
    
    public Node(String id) {
        this.id = id;
        weight = 1.0;
    }
    
    public Node(String id, double weight) {
        this.id = id;
        this.weight = weight;
    }
    
    public String getId() {
        return id;
    }
    
    public double getWeight() {
        return weight;
    }
    
    public void setWeight(double weight) {
        this.weight = weight;
    }
    
    @Override
    public int compareTo(Object o) {
        return (new Double(weight)).compareTo(new Double(((Node)o).getWeight()));
    }
    
    @Override
    public boolean equals(Object o) {
        return id.equals(((Node)o).getId());
    }
    
    @Override
    public String toString() {
        return id;
    }
    
}
