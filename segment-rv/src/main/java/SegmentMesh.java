import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Scanner;
import java.util.Set;
import java.io.File;
import java.io.FileNotFoundException;

public class SegmentMesh {

    static Set<Edge> makeEdges(Scanner sc) throws FileNotFoundException {
        Set<Edge> edges = new HashSet<>();
        int edgeIndex = 0;

        while (sc.hasNextLine()) {
            String line = sc.nextLine();
            String[] ln = line.split(" ");
            int n = Integer.parseInt(ln[0]);
            for (int i = 1; i <= n; i++) {
                int m = Integer.parseInt(ln[i]);
                edges.add(new Edge(edgeIndex, m, 1));
            }
            edgeIndex++;
        }

        return edges;
    }


    /**
     * @param n the number of nodes
     * @param m the number of edges
     * @param s the starting vertex (0 <= s < n)
     * @param edges the set of edges of the graph, with endpoints labelled between 1 and n inclusive.
     * @return The dictionary containing all the distances from s
     */
    public static Map<Integer, Integer> computeDistances(int n, int m, int s, Set<Edge> edges) {
        // Set up the network
        Map<Integer, LinkedList<Vertex>> network = new HashMap<>();
        Map<Integer, Integer> path = new HashMap<>();
        for (int i = 0; i <= n; i++) {
            network.put(i, new LinkedList<>());
            path.put(i, Integer.MAX_VALUE);
        }

        // Populate the network
        for (Edge e : edges) {
            network.get(e.from).add(new Vertex(e.to, e.weight));
        }

        // Meta data for the traversal
        PriorityQueue<Vertex> pq = new PriorityQueue<>();
        Set<Integer> cloud = new HashSet<>();

        // Initialise single source
        Vertex source = new Vertex(s, 0);
        path.put(s, 0);
        pq.add(source);

        while(cloud.size() != n) {
            Vertex u = pq.poll();

            if (!cloud.contains(u.index)) {
                cloud.add(u.index);

                for (Vertex v : network.get(u.index)) {
                    int tempDistance = path.get(u.index) + v.weight;
                    if (tempDistance < path.get(v.index)) {
                        path.put(v.index, tempDistance);
                    }

                    pq.add(new Vertex(v.index, 1, tempDistance));
                }
            }
        }

        return path;
    }

}
