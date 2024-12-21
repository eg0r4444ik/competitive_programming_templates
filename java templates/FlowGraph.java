import java.util.*;

public class FlowGraph{
        static class Edge{
            int to;
            long flow, cost, cap;

            public Edge(int to, long cap, long flow, long cost) {
                this.to = to;
                this.cap = cap;
                this.flow = flow;
                this.cost = cost;
            }
        }

        static class Pair{
            long flow, cost;

            public Pair(long flow, long cost) {
                this.flow = flow;
                this.cost = cost;
            }
        }

        public static List<Edge> edges;
        public static List<Integer>[] graph;
        static boolean[] used;
        static int[] lst;

        static int[] p;
        static int[] pe;
        static long[] dist;
        static int s, t, n;

        public static void init(int k) {
            n = k;
            edges = new ArrayList<>();
            graph = new List[n];
            used = new boolean[n];
            lst = new int[n];
            p = new int[n];
            pe = new int[n];
            dist = new long[n];
            s = 0;
            t = n-1;
            for(int i = 0; i < n; i++){
                graph[i] = new ArrayList<>();
            }
        }

        public static void addEdge(int from, int to, long cap, int flow, int cost){
            graph[from].add(edges.size());
            edges.add(new Edge(to, cap, flow, cost));
            graph[to].add(edges.size());
            edges.add(new Edge(from, 0, -flow, -cost));
        }

        public static void addBoundedEdge(int from, int to, long l, long r, int flow, int cost){
            graph[from].add(edges.size());
            edges.add(new Edge(to, r-l, flow, cost));
            graph[to].add(edges.size());
            edges.add(new Edge(from, 0, -flow, -cost));

            graph[from].add(edges.size());
            edges.add(new Edge(n-1, l, flow, cost));
            graph[n-1].add(edges.size());
            edges.add(new Edge(from, 0, -flow, -cost));

            graph[0].add(edges.size());
            edges.add(new Edge(to, l, flow, cost));
            graph[to].add(edges.size());
            edges.add(new Edge(0, 0, -flow, -cost));
        }

        public static long res(int x){
            return edges.get(x).cap-edges.get(x).flow;
        }

        public static long dfs(int x, long f, int k){
            if(used[x]){return 0;}
            if(x == t){return f;}
            used[x] = true;

            for(int e : graph[x]){
                if(res(e) < k){continue;}
                long pushed = dfs(edges.get(e).to, Math.min(f, res(e)), k);
                if(pushed != 0){
                    edges.get(e).flow += pushed;
                    edges.get(e^1).flow -= pushed;
                    return pushed;
                }
            }
            return 0;
        }

        static Pair augment(){
            spfa();
            if(p[t] == -1) return new Pair(0, 0);
            long mf = Integer.MAX_VALUE;
            int mc = 0;
            int curr = t;
            while(curr != s){
                int e = pe[curr];
                mf = Math.min(mf, res(pe[curr]));
                mc += edges.get(e).cost;
                curr = p[curr];
            }
            curr = t;
            while(curr != s){
                int e = pe[curr];
                edges.get(e).flow += mf;
                edges.get(e^1).flow -= mf;
                curr = p[curr];
            }
            return new Pair(mf, mc*mf);
        }

        static void spfa(){
            Arrays.fill(p, -1);
            Arrays.fill(pe, -1);
            Arrays.fill(dist, Integer.MAX_VALUE);
            for(int i = 0; i < n; i++){
                p[i] = -1;
                pe[i] = -1;
                dist[i] = Integer.MAX_VALUE;
            }
            dist[s] = 0;
            Queue<Integer> queue = new ArrayDeque<>();
            queue.add(s);
            boolean[] inq = new boolean[n];
            inq[s] = true;
            while(!queue.isEmpty()){
                int k = queue.peek();
                queue.poll();
                inq[k] = false;
                for(int e : graph[k]){
                    if(res(e) == 0) continue;
                    int to = edges.get(e).to;
                    long w = edges.get(e).cost;
                    if(dist[to] > dist[k] + w){
                        dist[to] = dist[k] + w;
                        p[to] = k;
                        pe[to] = e;
                        if(!inq[to]){
                            queue.add(to);
                            inq[to] = true;
                        }
                    }
                }
            }
        }

        static boolean bfs(){
            for(int i = 0; i < n; i++){
                dist[i] = Integer.MAX_VALUE;
            }
            dist[s] = 0;
            Queue<Integer> queue = new ArrayDeque<>();
            queue.add(s);
            while(!queue.isEmpty()){
                int k = queue.peek();
                queue.poll();
                for(int e : graph[k]){
                    if(res(e) == 0) continue;
                    int to = edges.get(e).to;
                    if(dist[to] > dist[k] + 1){
                        dist[to] = dist[k]+1;
                        queue.add(to);
                    }
                }
            }

            if(dist[t] != Integer.MAX_VALUE){
                return true;
            }
            return false;
        }

        static long dinicDfs(int x, long mx){
            if(x == t) return mx;
            int sum = 0;
            for(;lst[x] < graph[x].size(); lst[x]++){
                int e = graph[x].get(lst[x]);
                int to = edges.get(e).to;
                long r = res(e);
                if(r == 0 || dist[to] != dist[x]+1) continue;
                long pushed = dinicDfs(to, Math.min(mx-sum, r));
                sum += pushed;
                edges.get(e).flow += pushed;
                edges.get(e^1).flow -= pushed;
                if(sum == mx) break;
            }
            return sum;
        }

        static int dinic(){
            int flow = 0;
            while(true){
                if(!bfs()) break;
                Arrays.fill(lst, 0);
                while(true){
                    long f = dinicDfs(s, Integer.MAX_VALUE);
                    if(f == 0) break;
                    flow += f;
                }
            }
            return flow;
        }

        public static long maxFlow(){
            long ans = 0;
            int k = 1 << 30;
            while (k >= 1) {
                while (true) {
                    Arrays.fill(used, false);
                    long f = dfs(s, Integer.MAX_VALUE, k);
                    if (f == 0) {
                        break;
                    }
                    ans += f;
                }
                k /= 2;
            }
            return ans;
        }

        public static Pair minCostMaxFlow(){
            long flow = 0;
            long cost = 0;
            while(true){
                Pair pair = augment();
                if (pair.flow == 0){
                    break;
                }
                flow += pair.flow;
                cost += pair.cost;
            }
            return new Pair(flow, cost);
        }
    }
