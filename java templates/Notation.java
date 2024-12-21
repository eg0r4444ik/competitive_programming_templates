//import javax.xml.stream.events.StartDocument;
//import java.io.*;
//import java.util.*;
//public class Notation {
//    public static void main(String[] args) throws IOException {
//        InputStream inputStream = System.in;
//        OutputStream outputStream = System.out;
//        //inputStream = new FileInputStream(new File("input.txt"));
//        //outputStream = new FileOutputStream(new File("output.txt"));
//        InputReader in = new InputReader(inputStream);
//        PrintWriter out = new PrintWriter(outputStream);
//        solve(in, out);
//        out.close();
//    }
//
//    public static void solve(InputReader in, PrintWriter out) throws IOException {
//        Stack<String> stack = new Stack<>();
//        String inf = in.reader.readLine();
//        StringBuilder postf = new StringBuilder();
//        int[] skob = new int[inf.length()];
//        int count = 0;
//        int sum = 0;
//        for (int i = 0; i < inf.length(); i++) {
//            char cc = inf.charAt(i);
//            String s = String.valueOf(cc);
//            if (s.equals("(")) {
//                sum += 1;
//                skob[count] = sum;
//                count += 1;
//            }
//            if (s.equals(")")) {
//                sum -= 1;
//                skob[count] = sum;
//                count += 1;
//            }
//        }
//
//        for (int i = 0; i < count; i++) {
//            if (skob[i] < 0) {
//                out.println("WRONG");
//                return;
//            }
//        }
//
//        if (count > 0) {
//            if (skob[count - 1] != 0) {
//                out.println("WRONG");
//                return;
//            }
//        }
//
//        for (int i = 0; i < inf.length(); i++) {
//            char cc = inf.charAt(i);
//            String s = String.valueOf(cc);
//            if (!Character.isDigit(cc) && !s.equals("+") && !s.equals("-") && !s.equals("*") && !s.equals("(") && !s.equals(")") && !s.equals(" ")) {
//                out.println("WRONG");
//                return;
//            }
//        }
//
//        String[] sss = inf.split("\\s+");
//
//        for (int i = 0; i < sss.length; i++) {
//            if (i == 0) {
//                if (sss[i].equals("-")) {
//                    out.println("WRONG");
//                    return;
//                }
//            }
//            if (sss.length > 1) {
//                if (sss[i].matches("[0-9]") && sss[i + 1].matches("[0-9]")) {
//                    out.println("WRONG");
//                    return;
//                }
//            }
//
//            String[] ssss = inf.split("");
//            if (ssss[ssss.length - 1].matches("[*+-]")) {
//                out.println("WRONG");
//                return;
//            }
//            if (ssss.length == 1) {
//                break;
//            }
//            for (int j = 0; j < ssss.length - 1; j++) {
//                if (ssss[j].matches("[*+-]") && ssss[j + 1].matches("[*+-]")) {
//                    out.println("WRONG");
//                    return;
//                }
//                if (ssss[j].matches("[*+-]") && ssss[j + 1].matches("[)]")) {
//                    out.println("WRONG");
//                    return;
//                }
//                if (ssss[j].matches("[(]") && ssss[j + 1].matches("[*+-]")) {
//                    out.println("WRONG");
//                    return;
//                }
//            }
//        }
//
//        for (int i = 0; i < inf.length(); i++) {
//            Character c = inf.charAt(i);
//            if (Character.isDigit(c)) {
//                postf.append(c);
//            } else {
//                String s = String.valueOf(c);
//                if (!Objects.equals(s, " ")) {
//                    postf.append(" ");
//                }
//                switch (s) {
//                    case "+" -> {
//                        while (!stack.isEmpty()) {
//                            if (!Objects.equals(stack.peek(), "(")) {
//                                postf.append(" ").append(stack.pop()).append(" ");
//                            } else break;
//                        }
//                        stack.push("+");
//                    }
//                    case "-" -> {
//                        while (!stack.isEmpty()) {
//                            if (!Objects.equals(stack.peek(), "(")) {
//                                postf.append(" ").append(stack.pop()).append(" ");
//                            } else break;
//                        }
//                        stack.push("-");
//                    }
//                    case "*" -> {
//                        while (!stack.isEmpty()) {
//                            if (!Objects.equals(stack.peek(), "(") && !Objects.equals(stack.peek(), "+") && !Objects.equals(stack.peek(), "-"))
//                                postf.append(" ").append(stack.pop()).append(" ");
//                            else break;
//                        }
//                        stack.push("*");
//                    }
//                    case "(" -> {
//                        stack.push("(");
//                    }
//                    case ")" -> {
//                        while (!Objects.equals(stack.peek(), "(")) {
//                            postf.append(stack.pop()).append(" ");
//                        }
//                        stack.pop();
//                    }
//                }
//            }
//        }
//        while (!stack.isEmpty()) {
//            postf.append(" ").append(stack.pop()).append(" ");
//        }
//
//        String[] ss = postf.toString().trim().split("\\s+");
//
//        Stack<Integer> stack1 = new Stack<>();
//
//        for (String string : ss) {
//            if (!Objects.equals(string, "+") && !Objects.equals(string, "-") && !Objects.equals(string, "*")) {
//                stack1.push(Integer.valueOf(string));
//            } else {
//                int a = stack1.pop();
//                int b = stack1.pop();
//                switch (string) {
//                    case "+" -> stack1.push(a + b);
//                    case "-" -> stack1.push(b - a);
//                    case "*" -> stack1.push(a * b);
//                }
//            }
//        }
//        out.println(stack1.pop());
//
//    }
//
//    static class InputReader {
//        public BufferedReader reader;
//        public StringTokenizer tokenizer;
//
//        public InputReader(InputStream stream) {
//            reader = new BufferedReader(new InputStreamReader(stream), 32768);
//            tokenizer = null;
//        }
//
//        public String next() {
//            while (tokenizer == null || !tokenizer.hasMoreTokens()) {
//                try {
//                    tokenizer = new StringTokenizer(reader.readLine());
//                } catch (IOException e) {
//                    throw new RuntimeException(e);
//                }
//            }
//            return tokenizer.nextToken();
//        }
//
//        int nextInt() {
//            return Integer.parseInt(next());
//        }
//
//        long nextLong() {
//            return Long.parseLong(next());
//        }
//    }
//}
