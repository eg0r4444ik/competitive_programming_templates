class DSU{
        public static int[] p;
        public static int[] r;

        public DSU(int n) {
            p = new int[n+1];
            r = new int[n+1];
        }

        public static int find(int x){
            if(x == p[x]){
                return x;
            }
            p[x] = find(p[x]);
            return p[x];
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
