class Trie{

    public static class Node{
        int[] to = new int[26];
        boolean term;
    }

    Node[] data = new Node[10000000];
    int root = 0;
    int free = 1;

    public void Tree(){
        data[0] = new Node();
    }

    public void add(String s){
        int cur = root;
        char[] ch = s.toCharArray();
        for(char c : ch) {
            if(data[cur].to[c-97] == 0) {
                data[cur].to[c-97] = free++;
                cur = data[cur].to[c-97];
                data[cur] = new Node();
            } else {
                cur = data[cur].to[c - 97];
            }
        }
        data[cur].term = true;
    }

    public boolean find(String s){
        int cur = root;
        char[] ch = s.toCharArray();
        for(char c : ch) {
            if(data[cur].to[c-97] == 0) {
                return false;
            }
            cur = data[cur].to[c-97];
        }
        if(data[cur].term){
            return true;
        } else{
            return false;
        }
    }
}
