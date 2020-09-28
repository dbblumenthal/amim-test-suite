/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package de.mpg.mpiinf.ag1.kpm.graph;

import dk.sdu.kpm.graph.GeneEdge;
import dk.sdu.kpm.graph.GeneNode;
import dk.sdu.kpm.graph.KPMGraph;
import edu.uci.ics.jung.algorithms.shortestpath.UnweightedShortestPath;

/**
 *
 * @author nalcaraz
 */
public class LinearPathwayExtractor {
    
        
    UnweightedShortestPath<GeneNode, GeneEdge> uwsp;
    
    public LinearPathwayExtractor(KPMGraph g) {
        uwsp = new UnweightedShortestPath<GeneNode, GeneEdge>(g);        
    }
    
    
            
    
}
