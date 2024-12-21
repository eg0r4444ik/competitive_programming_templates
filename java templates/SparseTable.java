class SparseTable{
        long[][] table;
        int[] pows;
 
        public SparseTable(long[] a) {
            int n = a.length;
            this.table = new long[n][20];
            for(int i = 0; i < n; i++){
                table[i][0] = a[i];
            }
 
            for(int k = 1; k < 20; k++){
                for(int i = 0; i + (1 << k) - 1 < n; i++){
                    table[i][k] = Math.min(table[i][k-1], table[i + (1 << (k-1))][k-1]);
                }
            }
 
            pows = new int[n+1];
            int curr = 0;
            for(int i = 0; i <= n; i++){
                while(1 << (curr+1) <= i){
                    curr++;
                }
                pows[i] = curr;
            }
        }
 
        public long query(int l, int r){
            int k = pows[r-l+1];
            return Math.min(table[l][k], table[r - (1 << k) + 1][k]);
        }
    }
