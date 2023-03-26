import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Set;
import java.util.stream.Collectors;
import org.junit.jupiter.api.*;
import java.io.File;
import java.io.FileNotFoundException;

public class TestDijkstra {
    @Test
    public void segmentApex() throws FileNotFoundException {
        // Load in the data edges
        File f = new File("beutel.txt");
        Scanner sc = new Scanner(f);
        Set<Edge> edges = SegmentMesh.makeEdges(sc);

        // Define meta-parameters -> starting point
        int n = 938;
        int m = edges.size();
        int srcApex = 41;
        int neighboursApex = 228;

        // Obtain the dictionary
        var pathApex = SegmentMesh.computeDistances(n, m, srcApex, edges);
        // Only get the vertices that are close enough to the apex
        List<Integer> sortedApex = pathApex
            .entrySet()
            .stream()
            .sorted(Map.Entry.comparingByValue(Comparator.naturalOrder()))
            .filter(x -> x.getValue() <= neighboursApex)
            .map(Map.Entry::getKey)
            .collect(Collectors.toList());

        // Septal or Free Wall
        int srcSide = 605;
        int neighbourSide = 4;
        var pathSide = SegmentMesh.computeDistances(n, m, srcSide, edges);

        List<Integer> septalBody = pathSide
            .entrySet()
            .stream()
            .sorted(Map.Entry.comparingByValue(Comparator.naturalOrder()))
            .filter(x -> x.getValue() <= neighbourSide)
            .map(Map.Entry::getKey)
            .collect(Collectors.toList());

        // Take only the septal body values not contained by apex part
        var septal = septalBody
            .stream()
            .filter(x -> !sortedApex.contains(x))
            .collect(Collectors.toList());

        String region = "fb  ";
        for (var i : septal) {
            System.out.println(region + i);
        }
    }

    @Test
    public void inflowOutflow() throws FileNotFoundException {
        // Load in the data edges
        File f = new File("beutel.txt");
        Scanner sc = new Scanner(f);
        Set<Edge> edges = SegmentMesh.makeEdges(sc);

        // Define meta-parameters -> starting point
        int n = 938; // 162 938
        int m = edges.size();
        int s = 661; // 102 63

        // Obtain the dictionary
        var path = SegmentMesh.computeDistances(n, m, s, edges);
        int neighbourDistance = 4;

        // Sort the dictionary by values, and take only where the distance < R
        // where R is a specified radius
        List<Map.Entry<Integer, Integer>> sorted = path
            .entrySet()
            .stream()
            .sorted(Map.Entry.comparingByValue(Comparator.naturalOrder()))
            .filter(x -> x.getValue() <= neighbourDistance)
            .collect(Collectors.toList());

        // Print out the desired number of neighbours (only the index)
        String region = "ifw  ";
        for (var entry : sorted) {
            System.out.println(region + entry.getKey());
        }

    }

    @Test
    public void test1() throws FileNotFoundException {
        File f = new File("region.txt");
        Scanner sc = new Scanner(f);
        SegmentMesh.uniqueRegion(sc);
    }

}
