public class Z_Function {

    public static long[] zFunc(String str){
        char[] arr = str.toCharArray();
        long[] z = new long[arr.length];
        z[0] = 0;
        int l = 0, r = 0;
        for(int i = 1; i < z.length; i++){
            if(i <= r){
                z[i] = Math.min(z[i-l], r-i+1);
            }
            while(z[i] + i < z.length && arr[(int)z[i]] == arr[(int)z[i] + i]){
                z[i]++;
            }
            if(r < z[i] + i - 1){
                l = i;
                r = i + (int)z[i] - 1;
            }
        }
        return z;
    }

}
