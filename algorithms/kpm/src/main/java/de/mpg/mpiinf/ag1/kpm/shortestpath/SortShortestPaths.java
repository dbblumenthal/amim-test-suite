package de.mpg.mpiinf.ag1.kpm.shortestpath;

import java.util.Comparator;

import de.mpg.mpiinf.ag1.kpm.graph.LinearPath;

/**
 *
 * @author nalcaraz
 */
public enum SortShortestPaths {
    NODE_WEIGHT, EDGE_WEIGHT, LENGTH;
    
    public Comparator<LinearPath> getComparator() {
        switch(this) {
            case LENGTH:
                return LinearPath.LengthComparator;
            case NODE_WEIGHT:
                return LinearPath.NodeWeightComparator;
            case EDGE_WEIGHT:
                return LinearPath.EdgeWeightComparator;
            default:
                return LinearPath.LengthComparator;                
        }
    }
}
