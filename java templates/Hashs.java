import java.io.PrintWriter;
import java.util.ArrayList;

public class Hashs {
    // Дано две строки a и b. Вам требуется узнать, где в строке a можно найти строку b как подстроку и выписать все такие позиции.
    static void solve(Main.InputReader in, PrintWriter out) {
        String a = in.next();
        String b = in.next();
        if(a.length() < b.length()){
            out.println(0);
            return;
        }
        long p = 31;
        long m1 = 1000000007;
        long m2 = 1000000009;

        long[] pows1 = new long[a.length()];
        long pow1 = 1;
        for(int i = 0; i < pows1.length; i++){
            pows1[i] = pow1;
            pow1 = (pow1*p)%m1;
        }

        long[] pows2 = new long[a.length()];
        long pow2 = 1;
        for(int i = 0; i < pows2.length; i++){
            pows2[i] = pow2;
            pow2 = (pow2*p)%m2;
        }

        long[] pref1 = new long[a.length()+1];
        pref1[0] = 0;
        for(int i = 1; i < pref1.length; i++){
            pref1[i] = (pref1[i-1] + ((a.charAt(i-1)-'a'+1)*pows1[i-1])%m1)%m1;
        }

        long hash1 = 0;
        for(int i = 0; i < b.length(); i++){
            hash1 = (hash1 + ((b.charAt(i)-'a'+1)*pows1[i])%m1)%m1;
        }


        long[] pref2 = new long[a.length()+1];
        pref2[0] = 0;
        for(int i = 1; i < pref2.length; i++){
            pref2[i] = (pref2[i-1] + ((a.charAt(i-1)-'a'+1)*pows2[i-1])%m2)%m2;
        }

        long hash2 = 0;
        for(int i = 0; i < b.length(); i++){
            hash2 = (hash2 + ((b.charAt(i)-'a'+1)*pows2[i])%m2)%m2;
        }

        ArrayList<Integer> res = new ArrayList<>();
        for(int l = 1; l < a.length() - b.length() + 2; l++){
            long hash11 = (pref1[l+b.length()-1] - pref1[l-1] + m1)%m1;
            long hash12 = (pref2[l+b.length()-1] - pref2[l-1] + m2)%m2;
            long hash = (hash1 * pows1[l - 1] % m1) * m1 + (hash2 * pows2[l - 1] % m2);
            if(hash == hash11 * m1 + hash12){
                res.add(l);
            }
        }

        out.println(res.size());
        for(int  i : res){
            out.print(i + " ");
        }
    }
}
