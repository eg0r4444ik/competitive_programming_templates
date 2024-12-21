
import java.util.Scanner;

public class DisjointSetSystem {

    // Система непересекающихся множеств

    static Scanner scanner = new Scanner(System.in);
    public static int n = scanner.nextInt();
    public static int[] p = new int[n+1];
    public static int[] r = new int[n+1];

    public static int find(int x){
        if(x == p[x]){
            return x;
        }
        return find(p[x]);
    }

    public static void union(int x, int y) {
        x = find(x);
        y = find(y);
        if (r[x] > r[y]) {
            int t = x;
            x = y;
            y = t;
        }
        if(x != y) {
            r[y] += r[x];
            p[x] = y;
        }
    }

    public static boolean get(int x, int y){
        return find(x) == find(y);
    }

}