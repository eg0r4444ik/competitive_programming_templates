import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;

public class UnorderedGraph {
    // Неориентированный граф

    // Реализация графа
    public List<Node> nodes = new ArrayList<>();

    public class Node {
        public int value;
        public int distance = 0;
        public boolean visited;
        public List<Node> neighbors = new ArrayList<>();
        public Node(int value) {
            this.value = value;
        }
        public Node() {
        }
    }

    public Node addNode(int value) {
        Node newNode = new Node(value);
        nodes.add(newNode);
        return newNode;
    }

    public void addEdge(Node a, Node b) {
        a.neighbors.add(b);
        b.neighbors.add(a);
    }

    public void readGraph(int nodesCount, int edgesCount, Main.InputReader in, PrintWriter out){
        ArrayList<Node> nodes = new ArrayList<>();
        for (int i = 0; i < nodesCount; i++)
            nodes.add(new Node());
        for (int i = 0; i < edgesCount; i++) {
            int u = in.nextInt() - 1, v = in.nextInt() - 1;
            nodes.get(u).neighbors.add(nodes.get(v));
            nodes.get(v).neighbors.add(nodes.get(u));
        }
    }


    // Алгоритмы на графах

    // Обход в глубину dfs O(|V|+|E|)
    public void dfs(Node node) {
        node.visited = true;
        /* тут код при заходе в вершину */
        for (Node next: node.neighbors) {
            if (!next.visited){
                dfs(next);
            }
        }
        /* тут код при выходе из вершины */
    }

    // Обход в ширину bfs O(|V|+|E|)
    public void bfs(Node node) { // node - вершина начала обхода
        Queue<Node> queue = new LinkedList<>();
        node.visited = true;
        queue.add(node);
        while (queue.size() > 0) {
            Node current = queue.poll();
            /* Посетили вершину current */
            for (Node other: current.neighbors) {
                if (other.visited)
                    continue;
                other.visited = true;
                queue.add(other);
            }
        }
    }

    // Поиск кратчайшего пути
    public int distance(int S, int F) { // S - вершина из которой идут. F - вершина в которую идут
        Queue<Node> queue = new LinkedList<>();
        queue.add(nodes.get(S));
        nodes.get(S).visited = true;
        while (!queue.isEmpty()) {
            Node node = queue.poll();
            for (Node other : node.neighbors) {
                if (other.visited)
                    continue;
                other.visited = true;
                other.distance = node.distance + 1;
                queue.add(other);
            }
        }
        int answer = nodes.get(F).distance;
        return answer;
    }

    // Поиск цикла в графе
    public static boolean cycleDfs(Node node, Node prevNode) {
        if (node.visited)
            return true;
        node.visited = true;
        for (Node other: node.neighbors) {
            if (other != prevNode && cycleDfs(other, node))
                return true;
        }
        return false;
    }
    public static boolean hasCycle(List<Node> graph) {
        for (int i = 0; i < graph.size(); i++) {
            Node node = graph.get(i);
            if (!node.visited && cycleDfs(node, null))
                return true;
        }
        return false;
    }


}
