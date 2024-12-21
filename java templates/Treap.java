import java.util.Random;

public class Treap {
    // По неявным ключам
    static Random rnd = new Random(104360023);


    static class Node {
        int cnt, y;
        int val;
        Node l, r;

        public Node(int val) {
            y = rnd.nextInt();
            this.val = val;
            this.cnt = 1;
        }
    }

    static class Pair {
        Node l, r;

        public Pair(Node l, Node r) {
            this.l = l;
            this.r = r;
        }
    }

    static Node root = null;

    static int cnt(Node t) {
        return t != null ? t.cnt : 0;
    }

    static void upd_cnt(Node t) {
        if (t != null) {
            t.cnt = cnt(t.l) + cnt(t.r) + 1;
        }
    }

    static public Node merge(Node l, Node r) {
        if (l == null || r == null) {
            return l == null ? r : l;
        }
        if (l.y > r.y) {
            l.r = merge(l.r, r);
            upd_cnt(l);
            return l;
        } else {
            r.l = merge(l, r.l);
            upd_cnt(r);
            return r;
        }
    }

    static public Pair split(Node t, int x) {
        if (t == null) {
            return new Pair(null, null);
        }
        if (cnt(t.l) + 1 <= x) {
            Pair p = split(t.r, x - cnt(t.l) - 1);
            t.r = p.l;
            upd_cnt(t);
            return new Pair(t, p.r);
        } else {
            Pair p = split(t.l, x);
            t.l = p.r;
            upd_cnt(t);
            return new Pair(p.l, t);
        }
    }

    static void insert(int pos) {
        Pair p = split(root, pos);
        Node node = new Node(pos);
        root = merge(p.l, merge(node, p.r));
    }
}