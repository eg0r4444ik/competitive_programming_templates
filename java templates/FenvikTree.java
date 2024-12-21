import java.io.*;

public class FenvikTree {

    static long[] fenv;
    
    public static long fenvSum(int idx){
        int start = idx&(idx+1);
        long res = fenv[idx];
        if(start != 0){
            res += fenvSum(start-1);
        }
        return res;
    }

    static void init(long[] a){
        int n = a.length;
        long[] pref = new long[n];
        pref[0] = a[0];
        for(int i = 1; i < n; i++){
            pref[i] = pref[i-1] + a[i];
        }

        fenv = new long[n];
        fenv[0] = pref[0];
        for(int i = 1; i < n; i++){
            fenv[i] = (i & (i + 1)) == 0 ? pref[i] : pref[i] - pref[(i & (i + 1)) - 1];
        }
    }
}
