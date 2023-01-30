public class Vertex implements Comparable<Vertex> {
    int from;
    int index;
    int weight;
    int distance;

    public Vertex(int index, int w) {
        this.index = index;
        this.weight = w;
        this.distance = Integer.MAX_VALUE;
    }

    public Vertex(int index, int w, int d) {
        this.index = index;
        this.weight = w;
        this.distance = d;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Vertex vertex = (Vertex) o;
        return index == vertex.index && weight == vertex.weight && distance == vertex.distance;
    }

    @Override
    public int compareTo(Vertex o) {
        return this.distance - o.distance;
    }
}
