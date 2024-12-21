public class PrefixFunction {

    public static long[] prefixFunc(String str){
        char[] arr = str.toCharArray();
        long[] pref = new long[arr.length];
        pref[0] = 0;
        for(int i = 1; i < pref.length; i++){
            int j = (int)pref[i-1];
            while(j > 0 && arr[j] != arr[i]){
                j = (int)pref[j-1];
            }
            pref[i] = j + (arr[j] == arr[i] ? 1 : 0);
        }
        return pref;
    }

}
