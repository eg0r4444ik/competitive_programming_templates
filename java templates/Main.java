import com.sun.jdi.IntegerType;

import java.io.*;
import java.util.*;


public class Main {

    public static void main(String[] args) throws IOException, InterruptedException {
        InputStream inputStream = System.in;
        OutputStream outputStream = System.out;
//        inputStream = new FileInputStream(new File("invtrans.in"));
//        outputStream = new FileOutputStream(new File("invtrans.out"));
        InputReader in = new InputReader(inputStream);
        PrintWriter out = new PrintWriter(outputStream);
        int t = 1;
//        t = in.nextInt();
        for(int i = 0; i < t; i++){
            solve(in, out);
            out.flush();
        }
        out.close();
    }

    static void solve(InputReader in, PrintWriter out) throws InterruptedException {

    }

    static class InputReader {
        public BufferedReader reader;
        public StringTokenizer tokenizer;

        public InputReader(InputStream stream) {
            reader = new BufferedReader(new InputStreamReader(stream), 32768);
            tokenizer = null;
        }

        public String next() {
            while (tokenizer == null || !tokenizer.hasMoreTokens()) {
                try {
                    tokenizer = new StringTokenizer(reader.readLine());
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }
            return tokenizer.nextToken();
        }

        int nextInt() {
            return Integer.parseInt(next());
        }

        long nextLong() {
            return Long.parseLong(next());
        }
    }
}