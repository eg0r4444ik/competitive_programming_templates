import java.util.ArrayList;
import java.util.List;

public class Сombinatoric {
    // Комбинаторика

    // Проверка числа на простоту O(√N)
    public static boolean isPrime(int N) {
        for (int x = 2; x * x <= N; x++) {
            if (N % x == 0)
                return false;
        }
        return N >= 2;
    }


    static final int N = 10000000;
    static boolean[] isComposite = new boolean[N + 1];

    // Решето Эратосфена(если false, то число простое) O(N*logN)
    public static void resh() {
        for (int i = 2; i <= N; i++) {
            if (!isComposite[i]) {
                for (int j = i * 2; j <= N; j += i)
                    isComposite[j] = true;
            }
        }
    }

    // Все простые числа от 1 до N
    public static List<Integer> getPrimesList() {
        List<Integer> numbers = new ArrayList<>();
        for (int i = 2; i <= N; i++)
            if (!isComposite[i])
                numbers.add(i);
        return numbers;
    }

    static int[] divArray = new int[N + 1];

    // Решето Эратосфена с факторизацией чисел
    public static void reshFactorize() {
        for (int i = 2; i <= N; i++) {
            if (divArray[i] == 0) {
                for (int j = i; j <= N; j += i)
                    divArray[j] = i;
            }
        }
    }

    //Факторизация числа O(logN)
    public static List<Integer> factorize(int n) {
        List<Integer> numbers = new ArrayList<>();
        while (n > 1) {
            numbers.add(divArray[n]);
            n = n / divArray[n];
        }
        return numbers;
    }

    // Возведение в степень O(pow)
    public static long simplePow(int value, int pow, int M) {
        long result = 1;
        for (int i = 1; i <= pow; i++)
            result = (result * value) % M;
        return result;
    }

    // Бинарное возведение в степень O(log(pow))
    public static long fastPow(long num, long pow, long M) {
        if (pow == 0)
            return 1;
        if (pow % 2 == 0) {
            long result = fastPow(num, pow / 2, M);
            return (result * result) % M;
        }
        return (num * fastPow(num, pow - 1, M)) % M;
    }

    // Нахождение НОД (алгоритм Евклида) O(log(min(a,b))
    static long gcd(long a, long b) {
        if (b == 0)
            return a;
        else
            return gcd(b, a % b);
    }

    // Деление по модулю
    static long div(long a, long b, long M) {
        // (a/b) (mod M) == a*(b^(-1)) (mod M)
        // b^(-1) == b^(M-2) (mod M)
        return a * fastPow(b, M - 2, M);
    }

    // Расширенный алгоритм Евклида
    static class Gcd {
        int x, y, g;

        public Gcd(int x, int y, int g) {
            this.x = x;
            this.y = y;
            this.g = g;
        }
    }

    public static Gcd gcd(int a, int b) {
        if (b == 0) {
            return new Gcd(1, 0, a);
        }
        Gcd g = gcd(b, a % b);
        int x = g.y;
        int y = g.x - g.y * (a / b);
        return new Gcd(x, y, g.g);
    }

    // X0 = Xg*(c/g)
    // Y0 = Yg*(c/g)
    // X = X0 + k*(b/g)
    // Y = Y0 - k*(a/g)
}