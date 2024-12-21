class SegmentTree{

        static class Node{
            long modify;
            long function;

            public Node(long modify, long function) {
                this.modify = modify;
                this.function = function;
            }
        }

        long NEUTRAL_ELEMENT = 0;
        long NO_OPERATION = Long.MIN_VALUE;
        long MOD = 1000000007;

        Node[] tree;
        int size;

        public long op_modify(long a, long b, int len){
            if(b == NO_OPERATION){
                return a;
            }
            return b*len;
        }

        public long op_func(long a, long b){
            return a+b;
        }

        public void propagate(int x, int lx, int rx){
            if(tree[x].modify == NO_OPERATION || rx-lx == 1) return;
            int m = (lx+rx)/2;
            tree[2*x+1].modify = op_modify(tree[2*x+1].modify, tree[x].modify, 1);
            tree[2*x+1].function = op_modify(tree[2*x+1].function, tree[x].modify, m-lx);
            tree[2*x+2].modify = op_modify(tree[2*x+2].modify, tree[x].modify, 1);
            tree[2*x+2].function = op_modify(tree[2*x+2].function, tree[x].modify, rx-m);
            tree[x].modify = NO_OPERATION;
        }

        public void init(int n){
            size = 1;
            while(size < n){
                size *= 2;
            }

            tree = new Node[2*size-1];
        }

        public void build(long[] a, int x, int lx, int rx){
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

        public void build(long[] a){
            init(a.length);
            build(a, 0, 0, size);
        }


        public void modify(int l, int r, long v, int x, int lx, int rx){
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

        public void modify(int l, int r, long v){
            modify(l, r, v, 0, 0, size);
        }

        public long func(int l, int r, int x, int lx, int rx){
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

        public long func(int l, int r){
            return func(l, r, 0, 0, size);
        }
    }
