/**
 * 
 * NewtonGregory
 * 
 * <n>: inteiro, numero de pontos tabelados (E)
 * <x>: vetor, vetor de x (E)
 * <y>: vetor, vetor de y=f(x) (E)
 * <point>: double, ponto a ser interpolado (E)
 * 
*/

double NewtonGregory(int n, double x[], double y[], double point) {

    double result = 0;

    double dif_div[n][n];
    memset(dif_div, 0, sizeof(dif_div));

    for (int i = 0; i < n; i++) {
        dif_div[i][0] = y[i];
    }

    for (int j = 1; j < n; j++) {
        for (int i = 0; i < n - j; i++) {
            dif_div[i][j] = dif_div[i + 1][j - 1] - dif_div[i][j - 1];
        }
    }

    double h = x[1] - x[0];
    double u = (point - x[0]) / h;
    double term = 1; 

    result += y[0];

    for (int i = 1; i < n; i++) {
        term *= (u - i + 1) / i;
        result += term * dif_div[0][i];
    }

    return result;

}
