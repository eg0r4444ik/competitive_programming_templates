import java.io.*;
import java.util.*;


public class Run {

    public static void main(String[] args) throws IOException, InterruptedException {
        InputStream inputStream = System.in;
        OutputStream outputStream = System.out;
//        inputStream = new FileInputStream(new File("palindrome.in"));
//        outputStream = new FileOutputStream(new File("out1.txt"));
        InputReader in = new InputReader(inputStream);
        PrintWriter out = new PrintWriter(outputStream);
        int t = 1;
        t = in.nextInt();
        for(int i = 0; i < t; i++){
            solve(in, out);
            out.flush();
        }
        out.close();
    }

    static class SegmentTree{

        static class Node{
            long modify;
            long function;

            public Node(long modify, long function) {
                this.modify = modify;
                this.function = function;
            }
        }

        static long NEUTRAL_ELEMENT = Long.MIN_VALUE;
        static long NO_OPERATION = 0;
        static long MOD = 1000000007;

        static Node[] tree;
        static int size;

        public static long op_modify(long a, long b, int len){
            return (a+b);
        }

        public static long op_func(long a, long b){
            return Math.max(a, b);
        }

        public static void propagate(int x, int lx, int rx){
            if(tree[x].modify == NO_OPERATION || rx-lx == 1) return;
            int m = (lx+rx)/2;
            tree[2*x+1].modify = op_modify(tree[2*x+1].modify, tree[x].modify, 1);
            tree[2*x+1].function = op_modify(tree[2*x+1].function, tree[x].modify, m-lx);
            tree[2*x+2].modify = op_modify(tree[2*x+2].modify, tree[x].modify, 1);
            tree[2*x+2].function = op_modify(tree[2*x+2].function, tree[x].modify, rx-m);
            tree[x].modify = NO_OPERATION;
        }

        public static void init(int n){
            size = 1;
            while(size < n){
                size *= 2;
            }

            tree = new Node[2*size-1];
        }

        public static void build(long[] a, int x, int lx, int rx){
            if(rx-lx == 1){
                if(lx < a.length) {
                    tree[x] = new Node(NO_OPERATION, a[lx]);
                }else{
                    tree[x] = new Node(NO_OPERATION, 0);
                }
            }else{
                int m = (lx+rx)/2;
                build(a, 2*x+1, lx, m);
                build(a, 2*x+2, m, rx);
                tree[x] = new Node(NO_OPERATION, op_func(tree[2*x+1].function, tree[2*x+2].function));
            }
        }

        public static void build(long[] a){
            init(a.length);
            build(a, 0, 0, size);
        }


        public static void modify(int l, int r, long v, int x, int lx, int rx){
            propagate(x, lx, rx);
            if(l >= rx || r <= lx){
                return;
            }
            if(lx >= l && rx <= r){
                tree[x].modify = op_modify(tree[x].modify, v, 1);
                tree[x].function = op_modify(tree[x].function, v, rx-lx);
                return;
            }
            int m = (lx+rx)/2;
            modify(l, r, v, 2*x+1, lx, m);
            modify(l, r, v, 2*x+2, m, rx);
            tree[x].function = op_func(tree[2*x+1].function, tree[2*x+2].function);
            tree[x].function = op_modify(tree[x].function, tree[x].modify, rx-lx);
        }

        public static void modify(int l, int r, long v){
            modify(l, r, v, 0, 0, size);
        }

        public static long func(int l, int r, int x, int lx, int rx){
            propagate(x, lx, rx);
            if(l >= rx || r <= lx){
                return NEUTRAL_ELEMENT;
            }
            if(lx >= l && rx <= r){
                return tree[x].function;
            }
            int m = (lx+rx)/2;
            long res1 = func(l, r, 2*x+1, lx, m);
            long res2 = func(l, r, 2*x+2, m, rx);
            return op_func(res1, res2);
        }

        public static long func(int l, int r){
            return func(l, r, 0, 0, size);
        }

        public static long find(int l, int r, long k, long b, int x, int lx, int rx){
            propagate(x, lx, rx);
            if(l >= rx || r <= lx){
                return NEUTRAL_ELEMENT;
            }
            if(lx >= l && rx <= r){
                if(rx-lx == 1){
                    return Math.min(tree[x].function, k*(lx+1) + b);
                }
                int m = (lx+rx)/2;
                long midVal = k*m + b;
                if (tree[2*x+2].function >= midVal) {
                    long right = find(l, r, k, b, 2*x+2, m, rx);
                    return right;
                } else{
                    long left = find(l, r, k, b, 2*x+1, lx, m);
                    return Math.max(left, tree[2*x+2].function);
                }
            }
            int m = (lx+rx)/2;
            return Math.max(find(l, r, k, b, 2*x+1, lx, m), find(l, r, k, b, 2*x+2, m, rx));
        }

        public static long find(int l, int r, long k, long b){
            return find(l, r, k, b, 0, 0, size);
        }
    }

    static void solve(InputReader in, PrintWriter out) throws InterruptedException {
        int n = in.nextInt();
        int q = in.nextInt();
        long[] a = new long[n];
        for(int i = 0; i < n; i++){
            a[i] = in.nextInt();
        }

        SegmentTree.build(a);
        for(int i = 0; i < q; i++){
            String com = in.next();
            if(com.equals("+")){
                int l = in.nextInt()-1;
                int r = in.nextInt();
                long x = in.nextLong();
                SegmentTree.modify(l, r, x);
                for(int j = l; j < r; j++){
                    a[j] += x;
                }
            }else{
                int l = in.nextInt()-1;
                int r = in.nextInt();
                long k = in.nextLong();
                long b = in.nextLong();

                long res = Long.MIN_VALUE;
                for(int j = l; j < r; j++){
                    res = Math.max(res, Math.min(k*(j+1)+b, a[j]));
                }
                if(res != SegmentTree.find(l, r, k, b)){
                    for(int j = 0; j < 10; j++){
                        out.println();
                    }
                    out.println("NO!!!");
                    out.println(res + " != " + SegmentTree.find(l, r, k, b));
                    for(int j = 0; j < 10; j++){
                        out.println();
                    }
                }else{
                    out.println("YES");
                }
            }
        }
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