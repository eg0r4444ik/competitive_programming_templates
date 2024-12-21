import java.util.ArrayList;
import java.util.List;
import java.io.*;
import java.util.*;


public class Matching {
        public static void main(String[] args) throws IOException {
            InputStream inputStream = System.in;
            OutputStream outputStream = System.out;
//        inputStream = new FileInputStream(new File("bridges.in"));
//        outputStream = new FileOutputStream(new File("bridges.out"));
            InputReader in = new InputReader(inputStream);
            PrintWriter out = new PrintWriter(outputStream);
            int t = 1;
//        t = in.nextInt();
            for (int i = 0; i < t; i++) {
                solve(in, out);
            }
            out.close();
        }

        static int[] mt;
        static int n, k;
        static boolean[] used;
        static Node[] graph;

        static class Node{
            int val;
            List<Node> neighbours = new ArrayList<>();

            public Node(int val) {
                this.val = val;
            }
        }

        static void solve(InputReader in, PrintWriter out) {
            n = in.nextInt();
            k = in.nextInt();
            graph = new Node[n+k];
            used = new boolean[n+k];
            mt = new int[k];
            for(int i = 0; i < n+k; i++){
                graph[i] = new Node(i);
            }
            for(int i = 0; i < n; i++){
                int num = in.nextInt();
                while(num != 0){
                    graph[i].neighbours.add(graph[n+num-1]);
                    graph[n+num-1].neighbours.add(graph[i]);
                    num = in.nextInt();
                }
            }

            Arrays.fill(mt, -1);
            Arrays.fill(used, false);


            for (int i = 0; i < n; i++) {
                if(tryKuhn(graph[i])){
                    Arrays.fill(used, false);
                }
            }

            List<String> res = new ArrayList<>();
            for(int i = 0; i < mt.length; i++){
                if(mt[i] != -1){
                    res.add((mt[i]+1) + " " + (i+1));
                }
            }

            out.println(res.size());
            for(String s : res){
                out.println(s);
            }
        }

        static boolean tryKuhn(Node curr){
            if(used[curr.val]){
                return false;
            }
            used[curr.val] = true;
            for(Node node : curr.neighbours){
                int idx = node.val - n;
                if(mt[idx] == -1){
                    mt[idx] = curr.val;
                    return true;
                }else if(tryKuhn(graph[mt[idx]])){
                    mt[idx] = curr.val;
                    return true;
                }
            }

            return false;
        }

        static class InputReader {
            public BufferedReader reader;
            public StringTokenizer tokenizer;

            public InputReader(InputStream stream) {
                reader = new BufferedReader(new InputStreamReader(stream), 32768);
                tokenizer = null;
            }

            public String next() {
                while (tokenizer == null || !tokenizer.hasMoreTokens()) {
                    try {
                        tokenizer = new StringTokenizer(reader.readLine());
                    } catch (IOException e) {
                        throw new RuntimeException(e);
                    }
                }
                return tokenizer.nextToken();
            }

            int nextInt() {
                return Integer.parseInt(next());
            }

            long nextLong() {
                return Long.parseLong(next());
            }
        }
}
