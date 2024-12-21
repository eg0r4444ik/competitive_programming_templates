import java.util.ArrayList;
import java.util.List;

class LCA{

        static class Node{
            int val;
            boolean isVisited;
            List<Node> neighbours;

            public Node(int val) {
                this.val = val;
                isVisited = false;
                neighbours = new ArrayList<>();
            }
        }

        int[][] up;
        int[] dist;
        Node[] tree;

        public LCA(int n) {
            up = new int[n+1][31];

            dist = new int[n+1];
            dist[1] = 0;
            for(int i = 0; i <= 30; i++){
                up[1][i] = 0;
            }

            tree = new Node[n+1];
            for(int i = 0; i < n+1; i++){
                tree[i] = new Node(i);
            }
        }

        public void add(int a, int b){
            tree[a].neighbours.add(tree[b]);
            tree[b].neighbours.add(tree[a]);

            dist[b] = dist[a]+1;
            up[b][0] = a;
            for(int j = 1; j <= 30; j++){
                up[b][j] = up[up[b][j-1]][j-1];
            }
        }

        public int lca(int a, int b){
            if (dist[a] < dist[b]){
                int tmp = a;
                a = b;
                b = tmp;
            }

            for (int j = 30; j >= 0; j--){
                if (dist[up[a][j]] >= dist[b] && up[a][j] != 0){
                    a = up[a][j];
                }
            }
            if(a == b){
                return b;
            }else {
                for (int j = 30; j >= 0; j--) {
                    if (up[a][j] == up[b][j]) {
                        continue;
                    } else {
                        a = up[a][j];
                        b = up[b][j];
                    }
                }
                return up[a][0];
            }
        }
    }
