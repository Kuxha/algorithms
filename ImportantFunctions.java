
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Queue;
import java.util.Random;
import java.util.Stack;
import java.util.stream.IntStream;

public class ImportantFunctions {

    public static void main(String[] args) {

    }

    static class topologicalSort {

        private int V;   // No. of vertices 
        private LinkedList<Integer> adj[]; // Adjacency List 

        //Constructor 
        topologicalSort(int v) {
            V = v;
            adj = new LinkedList[v];
            for (int i = 0; i < v; ++i) {
                adj[i] = new LinkedList();
            }
        }

        void topologicalSortUtil(int v, boolean visited[],
                Stack stack) {
            // Mark the current node as visited. 
            visited[v] = true;
            Integer i;

            // Recur for all the vertices adjacent to this 
            // vertex 
            Iterator<Integer> it = adj[v].iterator();
            while (it.hasNext()) {
                i = it.next();
                if (!visited[i]) {
                    topologicalSortUtil(i, visited, stack);
                }
            }

            // Push current vertex to stack which stores result 
            stack.push(new Integer(v));
        }

        void TopologicalSort() {
            Stack stack = new Stack();

            // Mark all the vertices as not visited 
            boolean visited[] = new boolean[V];
            for (int i = 0; i < V; i++) {
                visited[i] = false;
            }

            // Call the recursive helper function to store 
            // Topological Sort starting from all vertices 
            // one by one 
            for (int i = 0; i < V; i++) {
                if (visited[i] == false) {
                    topologicalSortUtil(i, visited, stack);
                }
            }

            // Print contents of stack 
            while (stack.empty() == false) {
                System.out.print(stack.pop() + " ");
            }
        }

    }

    static class DisjointUnionSet {

        int[] parent, rank, setSize;

        public DisjointUnionSet(int n) {
            parent = new int[n];
            rank = new int[n];
            setSize = new int[n];
            for (int i = 0; i < parent.length; i++) {
                parent[i] = i;
                setSize[i] = 1;
            }
            Arrays.fill(rank, 1);
        }

        public void Union(int x, int y) {
            int x1 = findset(x);
            int y1 = findset(y);
            if (x1 == y1) {
                return;
            }
            if (rank[x1] > rank[y1]) {
                parent[y1] = x1;
                setSize[x1] += setSize[y1];
                // setSize[y1] = setSize[x1];
            } else {
                parent[x1] = y1;
                setSize[y1] += setSize[x1];
                // setSize[x1] = setSize[y1];
                if (rank[x1] == rank[y1]) {
                    rank[y1]++;
                }
            }
        }

        public int findset(int x) {
            if (x == parent[x]) {
                return x;
            }
            return parent[x] = findset(parent[x]);
        }

        public boolean isSameSet(int x, int y) {
            return findset(x) == findset(y);
        }
    }
//    ----------------------------------------------------------------------------
// IMPORTANT || COMMON  METHODS

    static void Luckybfs() {
        ArrayList<Long> lucky = new ArrayList<>();
        Queue<String> q = new LinkedList<String>();
        q.add("4");
        q.add("7");

        while (!q.isEmpty()) {
            String cur = q.remove();
            lucky.add(Long.parseLong(cur));
            if (Long.parseLong(cur) > (long) (1e15)) {
                continue;
            }
            q.add(cur + "4");
            q.add(cur + "7");
        }
    }

    static boolean Parity(long x) {
        long y = x ^ (x >> 1);
        y = y ^ (y >> 2);
        y = y ^ (y >> 4);
        y = y ^ (y >> 8);
        y = y ^ (y >> 16);

        if ((y & 1) > 0) {
            return false;
        }
        return true;
    }

    //Prime numbers till n 
    static class NextPowerOf2 {

        public int nextPowerOf2(int num) {
            if (num == 0) {
                return 1;
            }
            if (num > 0 && (num & (num - 1)) == 0) {
                return num;
            }
            while ((num & (num - 1)) > 0) {
                num = num & (num - 1);
            }
            return num << 1;
        }

        public static void main(String args[]) {
            NextPowerOf2 np = new NextPowerOf2();
            System.out.println(np.nextPowerOf2(4));
        }
    }
// PAIR

    static class Pair {

        int x;
        int y;

        // Constructor 
        public Pair(int x, int y) {
            this.x = x;
            this.y = y;
        }

        // how to initialize
//            Pair level[] = new Pair[n];
//            for (int i = 0; i < n; i++) {
//                int x = in.nextInt();
//                int y = in.nextInt();
//
//                level[i] = new Pair(x, y);
//
//            }
        // How to sort
//              PairSort ps = new PairSort();
//              ps.sort(level, 1);
    }

    static class PairSort {

        Pair[] sort(Pair arr[], int i) {
            // Comparator to sort the pair according to second element 
            Arrays.sort(arr, new Comparator<Pair>() {
                @Override
                public int compare(Pair p1, Pair p2) {
                    // 1 sees the first elements ans sorts them in increasing order 
                    // 2 sees the second elements ans sorts them in increasing order 
                    // 3 sees the first elements ans sorts them in decreasing order 
                    // 4 sees the second elements ans sorts them in decreasing order 

                    //rememeber 1 - 2 means ascending order. small to big
                    switch (i) {

                        case 1:
                            return p1.x - p2.x;
                        case 2:
                            return p1.y - p2.y;
                        case 3:
                            return p2.x - p1.x;
                        default:
                            return p2.y - p1.y;

                    }

                }
            });

            return arr;
        }
    }

    static ArrayList<Integer> lucky() {
        ArrayList<Integer> l = new ArrayList<>();
        l.add(4);
        l.add(7);
//            0 2 6 14 30 
        int index = 0;
        for (int i = 2; i < 10; i++) {
            int a = (int) Math.pow(2, i);
            int a0 = (int) Math.pow(2, i - 1);
            for (int j = 0; j < a / 2; j++) {

                int x = 4 * (int) Math.pow(10, i - 1) + l.get(index + j);

                l.add(x);

            }

            for (int j = 0; j < a / 2; j++) {

                int x = 7 * (int) Math.pow(10, i - 1) + l.get(index + j);

                l.add(x);

            }

            index = index + a0;
        }
        return l;
    }

    static boolean[] sieve(int n) {

        boolean prime[] = new boolean[n + 1];
        for (int i = 0; i < n; i++) {
            prime[i] = true;
        }

        for (int p = 2; p * p <= n; p++) {

            if (prime[p] == true) {

                for (int i = p * p; i <= n; i += p) {
                    prime[i] = false;
                }
            }
        }

        return prime;

    }

    static class Graph {

        int V;
        static LinkedList<Integer> adj[];

        Graph(int V) {
            this.V = V;

            adj = new LinkedList[V];

            for (int i = 0; i < V; i++) {
                adj[i] = new LinkedList<>();
            }
        }

        int BFS(int in1, int in2) {

            int visited[] = new int[V];

            LinkedList<Integer> queue = new LinkedList<Integer>();
            visited[in1] = 1;
            queue.add(in1);

            while (queue.size() != 0) {
                int s = queue.poll();
                Iterator<Integer> i = adj[s].listIterator();
                while (i.hasNext()) {
                    int n = i.next();
                    if (visited[n] == 0) {
                        visited[n] = visited[s] + 1;
                        queue.add(n);
                    }

                    if (n == in2) {

                        return visited[n] - 1;
                    }
                }
            }
            return 0;
        }

        static void addEdge(int src, int dest) {

            adj[src].add(dest);

            adj[dest].add(src);
        }

    }

    static void dfs(Graph g, int s, boolean visited[]) {
        int n = g.V;
        visited = new boolean[n + 1];
        visited[s] = true;
        Stack<Integer> st = new Stack<>();
        st.add(s);

        while (!st.isEmpty()) {

            int u = st.pop();
            Iterator<Integer> i = g.adj[u].iterator();
            while (i.hasNext()) {
                int v = i.next();
                if (!visited[v]) {
                    visited[v] = true;
                    st.push(v);
                }
            }

        }
    }

    boolean compare(int num1, int num2) {

        String s1 = num1 + "";
        String s2 = num2 + "";
        int c = 0;
        if (s1.charAt(0) != s2.charAt(0)) {
            c++;
        }
        if (s1.charAt(1) != s2.charAt(1)) {
            c++;
        }
        if (s1.charAt(2) != s2.charAt(2)) {
            c++;
        }
        if (s1.charAt(3) != s2.charAt(3)) {
            c++;
        }

        return (c == 1);
    }

    int shortestPath(int num1, int num2) {

        ArrayList<Integer> pset = new ArrayList<>();

        boolean seive[] = sieve(9999);

        for (int j = 1000; j <= 9999; j++) {
            if (seive[j]) {
                pset.add(j);
            }
        }

        Graph g = new Graph(pset.size());

        for (int i = 0; i < pset.size(); i++) {
            for (int j = i + 1; j < pset.size(); j++) {
                if (compare(pset.get(i), pset.get(j))) {
                    g.addEdge(i, j);
                }
            }
        }

        int in1 = 0, in2 = 0;
        for (int j = 0; j < pset.size(); j++) {
            if (pset.get(j) == num1) {
                in1 = j;
            }
        }
        for (int j = 0; j < pset.size(); j++) {
            if (pset.get(j) == num2) {
                in2 = j;
            }
        }

        return g.BFS(in1, in2);
    }

    private static void shuffleArray(int[] array) {
        int index, temp;
        Random random = new Random();
        for (int i = array.length - 1; i > 0; i--) {
            index = random.nextInt(i + 1);
            temp = array[index];
            array[index] = array[i];
            array[i] = temp;

        }
    }

// GCD || HCF
    static long GCD(long x, long y) {
        return y == 0 ? x : GCD(y, x % y);
    }

    // LCM 
    static long LCM(long a, long b) {
        return (a * b) / GCD(a, b);
    }

    // Returns the divisors of a number
    static ArrayList Divisors(long n) {

//        SO. If you need number of divisors . better to use this
//        same idea as Seive of erasthothenes
//            long[] div = new long[(int) MAX];
//            for (int i = 1; i < MAX; i++) {
//                for (int j = i; j < MAX; j += i) {
//                    div[j]++;
//                }
//            }
        ArrayList<Long> div = new ArrayList<>();

        for (long i = 1; i <= Math.sqrt(n); i++) {
            if (n % i == 0) {
                div.add(i);

                if (n / i != i) {
                    div.add(n / i);
                }
            }
        }
        return div;
    }

    // PERMUTATIONS
    static void permute(int[] a, int k) {
        if (k == a.length) {

            for (int i = 0; i < k; i++) {

                System.out.print(a[i] + " ");

            }

            System.out.println("");
        } else {
            for (int i = k; i < a.length; i++) {
                int temp = a[k];
                a[k] = a[i];
                a[i] = temp;

                permute(a, k + 1);
                // swapping a[k] and a[i] .
                temp = a[k];
                a[k] = a[i];
                a[i] = temp;

            }
        }

        // How to run 
        //int seq[] = {1, 19, 6}; .. let seq be the array of numbers to be permuted
        //permute(seq, 0); // call this function . 
        //0 -> all permutations ( permutation start from 1st number) 
        //1-> permutations start from 2nd number 
        //2-> permutation start from 3rd number
        // n-1 that array will be returned
    }

    // Sorting the indices .
    public int[] SortIndices(int x[]) {

        int[] indices = IntStream.range(0, x.length)
                .boxed().sorted((i, j) -> x[i] - x[j])
                .mapToInt(ele -> ele).toArray();
        return indices;
    }

    public boolean checkAllAdjacentCellsOfMatrix(char m[][], int r, int c) {

        // NOT CORNER CELLS btw // only adj. means left right up down
        // I always get confused because of exceptions of -
        // i=0 .. here we cant check for left cell
        for (int i = 0; i < r; i++) {

            for (int j = 0; j < c; j++) {
                if (m[i][j] == 'S' && j != c - 1) {
                    if (m[i][j + 1] == 'W') {

                        return false;
                    }
                }
                if (m[i][j] == 'S' && j != 0) {
                    if (m[i][j - 1] == 'W') {

                        return false;
                    }
                }
                if (m[i][j] == 'S' && i != r - 1) {
                    if (m[i + 1][j] == 'W') {

                        return false;
                    }
                }
                if (m[i][j] == 'S' && i != 0) {
                    if (m[i - 1][j] == 'W') {

                        return false;
                    }
                }
            }

        }
        return true;

    }
}
