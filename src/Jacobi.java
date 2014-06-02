import java.io.*;
import java.util.Locale;
import java.util.Scanner;

import static java.lang.Math.*;

public class Jacobi {
    public static void main(String[] args) throws IOException {
        Locale.setDefault(Locale.US);

        Scanner sc = new Scanner(new FileReader("input.txt"));
        PrintWriter pw = new PrintWriter(new FileWriter("output.txt"));
        int n = sc.nextInt();

        double[][] aa = new double[n][n];
        for(int i = 0; i < n; i++)
            for(int j = 0; j < n; j++)
                aa[i][j] = sc.nextDouble();

        double[][] ba = new double[n][1];
        for(int i = 0; i < n; i++)
            ba[i][0] = sc.nextDouble();

        Matrix A = new Matrix(aa);
        Matrix b = new Matrix(ba);
        Matrix x = new Matrix(new double[n][1]);

        Matrix D = A.getDiagonalized();
        Matrix D_1 = D.inverseDiagonalized();

        Matrix B = D_1.multiply(D.substract(A));
        if(abs(B.determinant()) > 1)
            pw.println("Process might not converge!\n");

        Matrix g = D_1.multiply(b);

        for(int i = 0; i < Constants.JACOBI_ITERATIONS; i++)
            x = B.multiply(x).add(g);

        pw.println("The solution is:");
        x.trace(pw);

        pw.println("\nCheck the solution at:");
        pw.println(Matrix.wolframLink(A, x, b));

        pw.flush();
        pw.close();
    }
}