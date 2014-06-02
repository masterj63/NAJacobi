import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;

import static java.lang.Math.*;

class Matrix {
    private final int COLS, ROWS;
    private final double[][] a;

    private static StringBuilder wolframText(Matrix A, Matrix x, Matrix b) {
        String format = "%." + Constants.DIGS + "f";
        StringBuilder sb = new StringBuilder();

        sb.append("abs(");

        sb.append('{');
        for (int i = 0; i < A.ROWS; i++) {
            sb.append('{');
            for (int j = 0; j < A.COLS; j++) {
                sb.append(String.format(format, A.a[i][j]));
                trim(sb);
                if (j < A.COLS - 1)
                    sb.append(", ");
            }
            sb.append('}');
            if (i < A.ROWS - 1)
                sb.append(", ");
        }
        sb.append('}');

        sb.append(" * ");

        sb.append('{');
        for (int i = 0; i < x.ROWS; i++) {
            sb.append(String.format(format, x.a[i][0]));
            trim(sb);
            if (i < x.ROWS - 1)
                sb.append(", ");
        }
        sb.append('}');

        sb.append(" - ");

        sb.append('{');
        for (int i = 0; i < b.ROWS; i++) {
            sb.append(String.format(format, b.a[i][0]));
            trim(sb);
            if (i < b.ROWS - 1)
                sb.append(", ");
        }
        sb.append('}');

        sb.append(')');

        return sb;
    }

    private static void trim(StringBuilder sb){
        while(sb.charAt(sb.length() - 1) == '0')
            sb.delete(sb.length() - 1, sb.length());
        if(sb.charAt(sb.length() - 1) == '.')
            sb.delete(sb.length() - 1, sb.length());
    }

    private static Map<Character, String> charMap = new HashMap<Character, String>();

    static {
        charMap.put('(', "%28");
        charMap.put(')', "%29");
        charMap.put('{', "%7B");
        charMap.put('}', "%7D");
        charMap.put('+', "%2B");
        charMap.put(',', "%2C");
        charMap.put(' ', "");
    }

    static String wolframLink(Matrix A, Matrix x, Matrix b) {
        StringBuilder text = wolframText(A, x, b);
        StringBuilder link = new StringBuilder();

        link.append("http://www.wolframalpha.com/input/?i=");
        for (Character c : text.toString().toCharArray())
            if (charMap.containsKey(c))
                link.append(charMap.get(c));
            else
                link.append(c);

        return link.toString();
    }

    double determinant() {
        if (ROWS != COLS)
            throw new IllegalArgumentException("This matrix is not square");

        double det = 1.0d;

        double[][] a = new double[this.a.length][];
        for(int i = 0; i < a.length; i++)
            a[i] = this.a[i].clone();

        for (int i = 0; i < ROWS - 1; i++) {
            int maxInRow = i;
            for (int j = i; j < ROWS; j++)
                if (abs(a[i][j]) > abs(a[i][maxInRow]))
                    maxInRow = j;

            int maxInCol = i;
            for (int j = i; j < ROWS; j++)
                if (abs(a[j][i]) > abs(a[maxInCol][i]))
                    maxInCol = j;

            if (abs(a[maxInCol][i]) > abs(a[i][maxInRow])) {
                if((maxInCol - i) % 2 == 1)
                    det *= -1.0d;
                for (int j = 0; j < ROWS; j++) {
                    double t = a[i][j];
                    a[i][j] = a[maxInCol][j];
                    a[maxInCol][j] = t;
                }
            } else {
                if((maxInRow - i) % 2 == 1)
                    det *= -1.0d;
                for (int j = 0; j < ROWS; j++) {
                    double t = a[j][i];
                    a[j][i] = a[j][maxInRow];
                    a[j][maxInRow] = t;
                }
            }

            for(int j = i + 1; j < ROWS; j++){
                double d = a[j][i] / a[i][i];
                for(int k = i; k < ROWS; k++)
                    a[j][k] -= a[i][k] * d;
            }
        }

        for(int i = 0; i < ROWS; i++)
            det *= a[i][i];

        return det;
    }

    void trace(PrintWriter pw) {
        int ceilLen = 0;
        for (double[] a : this.a)
            for (double d : a)
                ceilLen = Math.max(ceilLen, String.format("%.2f", d).split("\\D")[0].length());
        String format = String.format("%%%d.%df ", ceilLen + 1 + Constants.TRACE_DIGITS, Constants.TRACE_DIGITS);
        for (double[] a : this.a) {
            for (double d : a)
                pw.printf(format, d);
            pw.println();
        }
    }

    Matrix(double[][] a) {
        COLS = a[0].length;

        for (double[] d : a)
            assert d.length == COLS;

        ROWS = a.length;
        this.a = a;
    }

    Matrix getDiagonalized() {
        if (ROWS != COLS)
            throw new IllegalArgumentException("This matrix is not square");

        double[][] a = new double[ROWS][ROWS];
        for (int i = 0; i < ROWS; i++)
            a[i][i] = this.a[i][i];

        return new Matrix(a);
    }

    Matrix inverseDiagonalized() {
        if (ROWS != COLS)
            throw new IllegalArgumentException("This matrix is not square");

        for (int i = 0; i < ROWS; i++)
            for (int j = 0; j < ROWS; j++)
                if (i != j && Math.abs(a[i][j]) > Constants.EPS)
                    throw new IllegalArgumentException("Not a diagonal matrix");

        for (int i = 0; i < ROWS; i++)
            if (Math.abs(a[i][i]) < Constants.EPS)
                throw new IllegalArgumentException("This diagonal matrix is not singular");

        double[][] a = new double[ROWS][ROWS];
        for (int i = 0; i < ROWS; i++)
            a[i][i] = 1.0d / this.a[i][i];
        return new Matrix(a);
    }

    Matrix multiply(Matrix m) {
        if (COLS != m.a.length)
            throw new IllegalArgumentException("Matrix size mismatch");

        int rows = ROWS;
        int cols = m.a[0].length;

        double[][] r = new double[rows][cols];

        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                for (int k = 0; k < COLS; k++)
                    r[i][j] += a[i][k] * m.a[k][j];

        return new Matrix(r);
    }

    Matrix add(Matrix m) {
        if (COLS != m.COLS || ROWS != m.ROWS)
            throw new IllegalArgumentException("Matrix size mismatch");

        double[][] a = new double[ROWS][COLS];
        for (int i = 0; i < ROWS; i++)
            for (int j = 0; j < COLS; j++)
                a[i][j] = this.a[i][j] + m.a[i][j];

        return new Matrix(a);
    }

    Matrix substract(Matrix m) {
        if (COLS != m.COLS || ROWS != m.ROWS)
            throw new IllegalArgumentException("Matrix size mismatch");

        double[][] a = new double[ROWS][COLS];
        for (int i = 0; i < ROWS; i++)
            for (int j = 0; j < COLS; j++)
                a[i][j] = this.a[i][j] - m.a[i][j];

        return new Matrix(a);
    }
}