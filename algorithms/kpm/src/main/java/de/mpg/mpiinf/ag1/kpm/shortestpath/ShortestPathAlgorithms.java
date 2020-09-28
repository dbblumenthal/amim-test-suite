package de.mpg.mpiinf.ag1.kpm.shortestpath;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import de.mpg.mpiinf.ag1.kpm.Parameters;
import de.mpg.mpiinf.ag1.kpm.graph.Edge;
import de.mpg.mpiinf.ag1.kpm.graph.LinearPath;
import de.mpg.mpiinf.ag1.kpm.graph.Node;
import de.mpg.mpiinf.ag1.kpm.utils.FileToGraphParser;
import de.mpg.mpiinf.ag1.kpm.utils.StatisticsUtility;
import dk.sdu.kpm.KPMSettings;
import dk.sdu.kpm.graph.GeneEdge;
import dk.sdu.kpm.graph.GeneNode;
import dk.sdu.kpm.graph.KPMGraph;
import dk.sdu.kpm.graph.Result;
import edu.uci.ics.jung.graph.DirectedSparseGraph;
import edu.uci.ics.jung.graph.Graph;
import edu.uci.ics.jung.graph.UndirectedSparseGraph;
import edu.uci.ics.jung.graph.util.Pair;

public class ShortestPathAlgorithms {
	public static Parameters shortestPathways(List<Result> results, KPMGraph mainGraph, Parameters params, KPMSettings kpmSettings) {
		Graph<Node, Edge> g;
		if (params.IS_GRAPH_DIRECTED) {
			g = new DirectedSparseGraph<Node, Edge>();
		} else {
			g = new UndirectedSparseGraph<Node, Edge>();
		}

		for (int i = 0; i < params.PATHWAYS_SHORTEST_PATHS; i++) {
			Set<String> addedEdges = new HashSet<String>();
			Map<String, Node> addedNodes = new HashMap<String, Node>();

			Result result = results.get(i);

			List<GeneNode> gNodes = 
					new ArrayList<GeneNode>(result.getVisitedNodes().values());
			List<GeneEdge> gEdges = 
					new ArrayList<GeneEdge>(mainGraph.getConnectingEdges(gNodes));

			for (GeneEdge gEdge: gEdges) {
				Pair<GeneNode> pair = mainGraph.getEndpoints(gEdge);
				GeneNode gNode1 = pair.getFirst();
				String gNode1Id = gNode1.getNodeId();
				Node node1;
				if (addedNodes.containsKey(gNode1Id)) {
					node1 = addedNodes.get(gNode1Id);
				} else {
					node1 = new Node(gNode1Id, 
							gNode1.getAverageExpressedCasesNormalized());
					addedNodes.put(gNode1Id, node1);
				}

				GeneNode gNode2 = pair.getSecond();
				String gNode2Id = gNode2.getNodeId();
				Node node2;
				if (addedNodes.containsKey(gNode2Id)) {
					node2 = addedNodes.get(gNode2Id);
				} else {
					node2 = new Node(gNode2Id, 
							gNode2.getAverageExpressedCasesNormalized());
					addedNodes.put(gNode2Id, node2);
				}
				double total = 0.0;
				for (String expId: kpmSettings.NUM_CASES_MAP.keySet()) {
					int[] exp1 = gNode1.getDifferenceIntMap().get(expId);
					int[] exp2 = gNode2.getDifferenceIntMap().get(expId);
					int matches = 0;
					int ncases = kpmSettings.NUM_CASES_MAP.get(expId);
					for (int k = 0; k < ncases; k++) {
						if (Math.abs(exp1[k]) == Math.abs(exp2[k])) {
							matches++;
						}
					}
					double avg = (double)matches / (double)ncases;
					total += avg;
				}
				total = total / (double)kpmSettings.NUM_CASES_MAP.size();
				String edgeId = node1.getId() + "-" + node2.getId();
				Edge edge = new Edge(edgeId, total, node1, node2);
				if (!addedEdges.contains(edge.getId())) {
					g.addEdge(edge, node1, node2);
					addedEdges.add(edgeId);
				} 
			}
			DecimalFormat df = 
					StatisticsUtility.getIntegerFormat(kpmSettings.NUM_SOLUTIONS);
			String pathId = "-PATHWAY_" + df.format(i + 1);
			params = shortestPathways(g, pathId, params);
		}
		return params;
	}

	public static Parameters shortestPathways(Graph<Node, Edge> g, String pathwayId, Parameters params) {
		System.out.println("Vertices: " + g.getVertexCount());
		System.out.println("Edges: " + g.getEdgeCount());


		System.out.println("\n*********** EXTRACTING SHORTEST PATHS ***************\n");
		ShortestPathsDijkstra spd = new ShortestPathsDijkstra(g, params.IS_GRAPH_DIRECTED);
		ShortestPathways paths = spd.getShortestPathways();

		System.out.println("\n*********** WRITING RESULTS TO FILES ***************\n");


		if (params.GENERATE_SHORTEST_PATHS_STATS_FILE) {
			System.out.println("Generating shortests paths stats file...");

			String file = params.RESULTS_FOLDER + params.FILE_SEPARATOR + params.SHORTEST_PATHS_STATS_FILE
					+ pathwayId + "-" + params.SUFFIX  + params.FILE_EXTENSION;
			paths.writeShortestPathsFile(file, params.SORT_SHORTEST_PATHS_BY.getComparator());
		}

		if (params.GENERATE_SHORTEST_PATH_FILES) {
			System.out.println("Generating shortests paths individual files...");
			String file = params.RESULTS_FOLDER + params.FILE_SEPARATOR +
					params.SHORTEST_PATH_FILE + pathwayId;
			paths.writeIndividualShortestPaths(file,
					LinearPath.LengthComparator, 
					params.MINIMUM_LENGTH_SHORTEST_PATHS, params.SHORTEST_PATHS_REPORTED, params);
		}

		if (params.GENERATE_SHORTEST_PATHS_NODE_STATS_FILE) {
			System.out.println("Generating shortests paths node stats file...");
			String file = params.RESULTS_FOLDER + params.FILE_SEPARATOR + 
					params.SHORTEST_PATHS_NODE_STATS_FILE 
					+ pathwayId + "-" + params.SUFFIX + params.FILE_EXTENSION;
			paths.writeShortestPathsNodeStats(file, params);
		}

		if (params.GENERATE_SHORTEST_PATHS_EDGE_STATS_FILE) {
			System.out.println("Generating shortests paths edge stats file...");
			String file = params.RESULTS_FOLDER + params.FILE_SEPARATOR + 
					params.SHORTEST_PATHS_EDGE_STATS_FILE 
					+ pathwayId + "-" + params.SUFFIX  + params.FILE_EXTENSION;
			paths.writeShortestPathsEdgeStats(file, params);
		}

		System.out.println("\n*********** FINISHED SHORTEST PATHS "+" ***************\n");
		return params;
	}

	public static Parameters shortestPathways(String graphFile, Parameters params) {
		Graph<Node, Edge> g = null;
		System.out.println("\n*********** CREATING GRAPH ***************\n");
		if (graphFile.endsWith("sif")) {
			g = FileToGraphParser.parseSif(graphFile, params.IS_GRAPH_DIRECTED);
		} else {
			g = FileToGraphParser.parseTable(graphFile, params.IS_GRAPH_DIRECTED, false, false);
		}

		return shortestPathways(g, "", params);
	}

}
