/**
 * 
 * Newton
 * 
 * <n>: inteiro, numero de pontos tabelados (E)
 * <x>: vetor, vetor de x (E)
 * <y>: vetor, vetor de y=f(x) (E)
 * <point>: double, ponto a ser interpolado (E)
 *  
*/

double Newton (int n, double x[], double y[], double point) {

    double res = 0; // Valor do polinômio interpolado
    double aux;

    // Inicializa o array de diferenças divididas
    double dif_div[n][n];
    memset(dif_div, 0, sizeof(dif_div));

    // Calcula as diferenças divididas
    for (int i = 0; i < n; i++) {
        dif_div[i][0] = y[i];
    }

    for (int j = 1; j < n; j++) {
        for (int i = 0; i < n - j; i++) {
            dif_div[i][j] = (dif_div[i + 1][j - 1] - dif_div[i][j - 1]) / (x[i + j] - x[i]);
        }
    }

    // Calcula o valor do polinômio interpolado usando as diferenças divididas
    for (int i = 0; i < n; i++) {
        aux = dif_div[0][i];
        for (int j = 0; j < i; j++) {
            aux *= (point - x[j]);
        }
        res += aux;
    }

    return res;

}

// nao era pra chegar ate aqui...