#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MTXSIZE 105
#define DEBUG false
#define MAXGRAU 7

typedef double tabela[2][MTXSIZE];
typedef double vetorY[MTXSIZE];
typedef double vetorA[MAXGRAU];

// FUNCOES AUXILIARES

/**
 * cls
*/
void cls()
{
    // limpa a tela
    printf("\e[1;1H\e[2J");
}

// ROTINAS LEGADAS

void printVetor(int ordem, double matriz[MTXSIZE])
{
    for(int i = 0; i < ordem; i++){
        printf("%4.4lf ", matriz[i]);
    }
    printf("\n");
}

void extraiCol(int ordem, double matriz[MTXSIZE][MTXSIZE], int col, double vetor[MTXSIZE])
{
    for(int i = 0; i < ordem; i++){
        vetor[i] = matriz[i][col];
    }
}

void copiaMatriz(double matriz[MTXSIZE][MTXSIZE], double saida[MTXSIZE][MTXSIZE])
{
    memcpy(saida, matriz, MTXSIZE * MTXSIZE * sizeof(double));
    return;
}

void cortaMatriz(int ordem, double matriz[MTXSIZE][MTXSIZE], int lin, int col, double saida[MTXSIZE][MTXSIZE])
{
    if(ordem == 0){
        return;
    }

    double result[MTXSIZE][MTXSIZE];

    // indentificam se a linha ou coluna ja foram pulados para esse loop
    int pLin = 0, pCol;

    for(int i = 0; i < ordem; i++){
        // pula essa linha
        if(i == lin){
            pLin = 1;
            continue;
        }

        pCol = 0;

        for(int j = 0; j < ordem; j++){
            // pula essa coluna
            if(j == col){
                pCol = 1;
                continue;
            }

            result[i - pLin][j - pCol] = matriz[i][j];
        }
    }

    copiaMatriz(result, saida);
}

double Determinante (int ordem, double matriz[MTXSIZE][MTXSIZE]) {

    // saida
    if(ordem == 1){
        return matriz[0][0];
    }

    double mtAux[MTXSIZE][MTXSIZE], det = 0;
    int lin = 0;

    // determinante calculado baseado na primeira linha
    for(int col = 0; col < ordem; col++){

        // gera matriz auxiliar para calculo da menor principal
        cortaMatriz(ordem, matriz, lin, col, mtAux);
        
        // contabiliza esse cofator
        det += matriz[lin][col] * Determinante(ordem - 1, mtAux) * (lin + col % 2 == 0 ? 1 : -1);
        
    }
    
    return det;
    
}

int definida(int ordem, double matriz[MTXSIZE][MTXSIZE])
{
    double det;
    double mAux[MTXSIZE][MTXSIZE];
    copiaMatriz(matriz, mAux);

    while(ordem != 0){
        // checa determinante
        det = Determinante(ordem, mAux);
        if(det == 0){
            return 0;
        }

        // parte pra proxima submatriz
        cortaMatriz(ordem, mAux, ordem - 1, ordem - 1, mAux);
        ordem --;
    }

    return 1;
}

int SistemaTriangularInferior (int ordem, double matriz[MTXSIZE][MTXSIZE], double vetorI[MTXSIZE], double vetorS[MTXSIZE]) {

    // CONVERGENCIA
    if(!definida(ordem, matriz)){
        return 0;
    }

    vetorS[0] = vetorI[0]/matriz[0][0];
    double sum;
    
    if ( ordem > 1 ) {

        for ( int i = 1; i < ordem; i++) {
        	
        	sum = 0;

            for ( int j = 0; j < i; j++ ) {
                sum += matriz[i][j] * vetorS[j];
            }

            vetorS[i] = (vetorI[i] - sum) / matriz[i][i];

        }

    }

    return 1;

}

int SistemaTriangularSuperior (int ordem, double matriz[MTXSIZE][MTXSIZE], double vetorI[MTXSIZE], double vetorS[MTXSIZE]) {

    // CONVERGENCIA
    if(!definida(ordem, matriz)){
        return 0;
    }

    vetorS[ordem-1] = vetorI[ordem-1]/matriz[ordem-1][ordem-1]; double sum;
    
    if ( ordem > 1 ) {

        for ( int i = ordem-2; i >= 0; i-- ) {
        	
        	sum = 0;
        	
            for ( int j = i+1; j < ordem; j++ ) {
                sum += matriz[i][j] * vetorS[j];
            }

            vetorS[i] = (vetorI[i] - sum) / matriz[i][i];

        }

    }

    return 1;

}

int DecomposicaoLU(int ordem, double matriz[MTXSIZE][MTXSIZE], double vetorI[MTXSIZE], double vetorS[MTXSIZE]) {
    if (!definida(ordem, matriz)) {
        return 0;
    }

    double matrizL[MTXSIZE][MTXSIZE], matrizU[MTXSIZE][MTXSIZE];
    double vetorAux[MTXSIZE];
    double sum = 0;

    for (int i = 0; i < ordem; i++) {
        for (int j = i; j < ordem; j++) {
            sum = 0;
            for (int k = 0; k < i; k++) {
                sum += matrizL[i][k] * matrizU[k][j];
            }
            matrizU[i][j] = matriz[i][j] - sum;
        }
        for (int j = i; j < ordem; j++) {
            if (i == j) {
                matrizL[i][j] = 1;
                continue;
            }
            sum = 0.0;
            for (int k = 0; k < i; k++) {
                sum += matrizL[j][k] * matrizU[k][i];
            }
            matrizL[j][i] = (matriz[j][i] - sum) / matrizU[i][i];
        }
    }

    int buf;
    buf = SistemaTriangularInferior(ordem, matrizL, vetorI, vetorAux);
    buf = SistemaTriangularSuperior(ordem, matrizU, vetorAux, vetorS);

    return 1;
}

// ROTINAS DE INTERPOLACAO E AJUSTE DE CURVAS

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
double Newton(int n, tabela points, double point) {

    double res = 0;
    double aux;

    double dif_div[n][n];
    memset(dif_div, 0, sizeof(dif_div));

    for (int i = 0; i < n; i++) {
        dif_div[i][0] = points[1][i];
    }

    for (int j = 1; j < n; j++) {
        for (int i = 0; i < n - j; i++) 
            dif_div[i][j] = (dif_div[i + 1][j - 1] - dif_div[i][j - 1]) / (points[0][i + j] - points[0][i]);
        //
    }

    for (int i = 0; i < n; i++) {
        aux = dif_div[0][i];
        for (int j = 0; j < i; j++) 
            aux *= (point - points[0][j]);
        //
        res += aux;
    }

    return res;

}

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
double NewtonGregory(int n, tabela points, double point) {

    double result = 0;
    double dif_div[n][n];

    memset(dif_div, 0, sizeof(dif_div));

    for (int i = 0; i < n; i++)
        dif_div[i][0] = points[1][i];
    //

    for (int j = 1; j < n; j++) {
        for (int i = 0; i < n - j; i++) {
            dif_div[i][j] = dif_div[i + 1][j - 1] - dif_div[i][j - 1];
        }
    }

    double h = points[0][1] - points[0][0];
    double u = (point - points[0][0]) / h;
    double term = 1;
    result += points[1][0];

    for (int i = 1; i < n; i++) {
        term *= (u - i + 1) / i;
        result += term * dif_div[0][i];
    }

    return result;

}

/**
 * 
 * Retorna o coeficiente de determinação entre os pontos
 * tabelados e os pontos ajustados. 
 * <arg-1>: inteiro, número de pontos tabelados (E) 
 * <arg-2>: tabela, pontos tabelados (E) 
 * <arg-3>: vetorY, valores Y ajustado (E)
 * 
*/
double CoefDeterminacao(int npts, tabela pts, vetorY y) {
    double r;
    double sumY = 0, sumYQ = 0, sumErrQ = 0;

    for(int i = 0; i < npts; i++){
        sumY += pts[1][i];
        sumYQ += pow(pts[1][i], 2);
        sumErrQ += pow(pts[1][i] - y[i], 2);
    }
    
    r = 1 - ((npts * sumErrQ) / ((npts * sumYQ) - pow(sumY, 2)));
    return r;

}

/**
 * Ajusta os pontos tabelados a uma reta da forma y = a0 + a1x.
 * <arg-1>: inteiro, número de pontos tabelados (E)
 * <arg-2>: tabela, pontos tabelados (E)
 * <arg-3>: real, termo independente da equação de reta ajustada (a0) (S)
 * <arg-4>: real, coeficiente de grau 1 da equação de reta ajustada (a1) (S)
 * <arg-5>: vetorY, valores Y ajustados (S)
 * <arg-6>: real, coeficiente de determinação entre os pontos tabelados e os pontos ajustados (S)
*/
void AjusteReta(int npts, tabela pts, double *a0, double *a1, vetorY y, double *coef) {

    double sumX = 0, sumXQ = 0, sumY = 0, sumXY = 0;

    for(int i = 0; i < npts; i++){
        sumX += pts[0][i];
        sumXQ += pow(pts[0][i], 2);
        sumY += pts[1][i];
        sumXY += pts[0][i] * pts[1][i];        
    }  

    *a1 = ((npts * sumXY) - (sumX * sumY)) / ((npts * sumXQ) - pow(sumX, 2)); 
    *a0 = (sumY - (*a1 * sumX)) / npts; 


    for(int i = 0; i < npts; i++){
        y[i] = *a0 + (*a1 * pts[0][i]);
    }  

    *coef = CoefDeterminacao(npts, pts, y);
    return;

}

/**
 * AjustePolinomio
 * 
 * <arg-1>: inteiro, número de pontos tabelados (E)
 * <arg-2>: inteiro, grau do polinômio a ser ajustado (E)
 * <arg-3>: tabela, pontos tabelados (E)
 * <arg-4>: vetorA, coeficientes do polinômio ajustado (S)
 * <arg-5>: vetorY, valores Y ajustados (S)
 * <arg-6>: real, coeficiente de determinação entre os pontos tabelados e os pontos ajustados (S)
 * 
*/
void AjustePolinomio(int npts, int grau, tabela tab, vetorA a, vetorY y, double *coefDeterminacao) {
    
    double matriz[MTXSIZE][MTXSIZE];
    double vetorI[MTXSIZE];

    // Gera matriz e vetor independentedo sistema
    for (int i = 0; i <= grau; i++) {

        for (int j = 0; j <= grau; j++) {
            double sum = 0;
            for (int k = 0; k < npts; k++) 
                sum += pow(tab[0][k], i + j);
            //
            matriz[i][j] = sum;
        }

        double sum = 0;
        for (int k = 0; k < npts; k++)
            sum += tab[1][k] * pow(tab[0][k], i);
        //
        vetorI[i] = sum;

    }

    DecomposicaoLU(grau + 1, matriz, vetorI, a);

    // Calcula Y ajustados
    for (int k = 0; k < npts; k++) {
        y[k] = 0;
        for (int i = 0; i <= grau; i++)
            y[k] += a[i] * pow(tab[0][k], i);
        //
    }

    *coefDeterminacao = CoefDeterminacao(npts, tab, y);

}

/**
 * Ajusta os pontos tabelados a uma curva exponencial da forma y = ab^x
 * <arg-1>: inteiro, número de pontos tabelados (E)
 * <arg-2>: tabela, pontos tabelados (E)
 * <arg-3>: real, 1o coeficiente da equação exponencial ajustada (a) (S)
 * <arg-4>: real, 2o coeficiente da equação exponencial ajustada (b) (S)
 * <arg-5>: vetorY, valores Y ajustados (S)
 * <arg-6>: real, coeficiente de determinação entre os pontos tabelados e os pontos ajustados (S)
*/
void AjusteExponencial (int npts, tabela pts, double *a, double *b, vetorY y, double *coefDet) {

    double coef, a1, a0, yhat[MTXSIZE];
    double sumX = 0, sumXQ = 0, sumY = 0, sumXY = 0;

    for(int i = 0; i < npts; i++){
        pts[1][i] = log(pts[1][i]);

        sumX += pts[0][i];
        sumXQ += pow(pts[0][i], 2);
        sumY += pts[1][i];
        sumXY += pts[0][i] * pts[1][i];        
    }  

    a1 = ((npts * sumXY) - (sumX * sumY)) / ((npts * sumXQ) - pow(sumX, 2)); 
    a0 = (sumY - (a1 * sumX)) / npts; 
    *a = exp(a0);
    *b = exp(a1);

    for(int i = 0; i < npts; i++){
        y[i] = (*a) * pow(*b, pts[0][i]);
        yhat[i] = a0 + a1 *  pts[0][i];
    }  

    *coefDet = CoefDeterminacao(npts, pts, yhat);
    return;
    
}
// FLUXO DE PROGRAMA

int main() {

    int input = -1, num = 0, grau = 0;
    double point = 0;

    tabela tabelados; vetorY ajustados; vetorA a;
    double precisao = 0; double a0 = 0; double a1 = 0; double coef = 0;

    while(1){

        /* menu */ {
            printf("\n\tMETODOS DE SISTEMAS LINEARES\n\n");
            printf("\t01 - Interpolacao Metodo de Newton\n");
            printf("\t02 - Interpolacao Metodo de Newton-Gregory\n");
            printf("\t03 - Coeficiente de Determinação\n");
            printf("\t04 - Ajuste Linear\n");
            printf("\t05 - Ajuste Polinomial\n");
            printf("\t06 - Ajuste Exponencial \n");
            printf("\t00 - Sair\n\n");
        }

        do{ // ler input
            printf("\tDigite uma opcao:\t");
            scanf("%d", &input);
        }while(input < 0 || input > 6);

        switch(input){ // fluxo para cada metodo
		    
	        case 1:{ // Newton

	            cls();

                printf("\n\tINTERPOLACAO METODO DE NEWTON\n\n");
				printf("\tDigite o numero de pontos tabelados (Max: 100):\t"); scanf("%d", &num);

                printf("\n\tDigite os valores de x tabelados:\n\n\t{ ");
                for ( int i = 0; i < num; i++ ) {
                    scanf("%lf",&tabelados[0][i]);
				}

                printf("\n\tDigite os valores de y tabelados:\n\n\t{ ");
                for ( int i = 0; i < num; i++ ) {
                    scanf("%lf",&tabelados[1][i]);
				}

                printf("\n\tDigite o ponto a ser interpolado:\t"); scanf("%lf", &point);
                printf("\n\tResultado:\t%.4lf\n", Newton(num, tabelados, point));
                
            }
	        break;

            case 2:{ // Newton Gregory

	            cls();

                printf("\n\tINTERPOLACAO METODO DE NEWTON GREGORY\n\n");
				printf("\tDigite o numero de pontos tabelados (Max: 100):\t"); scanf("%d", &num);

                printf("\n\tDigite os valores de x tabelados:\n\n\t{ ");
                for ( int i = 0; i < num; i++ ) {
                    scanf("%lf",&tabelados[0][i]);
				}

                printf("\n\tDigite os valores de y tabelados:\n\n\t{ ");
                for ( int i = 0; i < num; i++ ) {
                    scanf("%lf",&tabelados[1][i]);
				}

                printf("\n\tDigite o ponto a ser interpolado:\t"); scanf("%lf", &point);
                printf("\n\tResultado:\t%.4lf\n", NewtonGregory(num, tabelados, point));
                
            }
	        break;

            case 3:{ // Coeficiente de Determinação

	            cls();

                printf("\n\tCOEFICIENTE DE DETERMINACAO\n\n");
				printf("\tDigite o numero de pontos tabelados (Max: 100):\t"); scanf("%d", &num);

                printf("\n\tDigite os valores de x tabelados:\n\n\t{ ");
                for ( int i = 0; i < num; i++ ) {
                    scanf("%lf",&tabelados[0][i]);
				}

                printf("\n\tDigite os valores de y tabelados:\n\n\t{ ");
                for ( int i = 0; i < num; i++ ) {
                    scanf("%lf",&tabelados[1][i]);
				}

                printf("\n\tDigite o vetor dos valores de y ajustados:\n\n\t{ ");
                for ( int i = 0; i < num; i++ ) {
                    scanf("%lf",&ajustados[i]);
                }            

                printf("\n\tResultado:\t%.4lf\n", CoefDeterminacao(num, tabelados, ajustados));
                
            }
	        break;

            case 4:{ // Ajuste Linear

	            cls();

                printf("\n\tAJUSTE LINEAR\n\n");
				printf("\tDigite o numero de pontos tabelados (Max: 100):\t"); scanf("%d", &num);

                printf("\n\tDigite os valores de x tabelados:\n\n\t{ ");
                for ( int i = 0; i < num; i++ ) 
                    scanf("%lf",&tabelados[0][i]);
				//

                printf("\n\tDigite os valores de y tabelados:\n\n\t{ ");
                for ( int i = 0; i < num; i++ ) 
                    scanf("%lf",&tabelados[1][i]);
				//

                AjusteReta(num, tabelados, &a0, &a1, ajustados, &coef);

                printf("\n\tResultados:\t termo independente da eq. da reta = %.4lf", a0);
                printf("\n\tCoeficiente de primeiro grauda eq. da reta = %.4lf",a1); 
                printf("\n\tCoeficiente de determinacao entre os pontos tabelados e ajustados = %.4lf\n",coef);

                printf("\n\tValores Y ajustados:\n\n\t{ ");
                for ( int i = 0; i < num; i++ )
                    printf("%.4lf ", ajustados[i]);
                //
                
            }
	        break;

            case 5:{ // Ajuste Polinomial

	            cls();

                printf("\n\tAJUSTE POLINOMIAL\n\n");
				printf("\tDigite o numero de pontos tabelados (Max:100):\t"); scanf("%d", &num);
                printf("\tDigite o grau do polinomio desejado (Max:5):\t"); scanf("%d", &grau);

                printf("\n\tDigite os valores de x tabelados:\n\n\t{ ");
                for ( int i = 0; i < num; i++ ) 
                    scanf("%lf",&tabelados[0][i]);
				//

                printf("\n\tDigite os valores de y tabelados:\n\n\t{ ");
                for ( int i = 0; i < num; i++ ) 
                    scanf("%lf",&tabelados[1][i]);
			    //

                AjustePolinomio(num, grau, tabelados, a,ajustados, &coef);
                printf("\n\tResultados:\t coeficientes do polinomio ajustado = { ");
                for ( int i = 0; i <= grau; i++ )
                    printf("%.4lf ", a[i]);
                //
                printf("\n\tValores Y ajustados:\n\n\t{ ");
                for ( int i = 0; i < num; i++ )
                    printf("%.4lf ", ajustados[i]);
                //
                printf("}\n\tCoeficiente de determinacao entre os pontos tabelados e ajustados = %.4lf\n",coef);

            }
	        break;

            case 6:{ // Ajuste Exponencial

	            cls();

                printf("\n\tAJUSTE EXPONENCIAL\n\n");
				printf("\tDigite o numero de pontos tabelados (Max: 100):\t"); scanf("%d", &num);

                printf("\n\tDigite os valores de x tabelados:\n\n\t{ ");
                for ( int i = 0; i < num; i++ ) 
                    scanf("%lf",&tabelados[0][i]);
				//

                printf("\n\tDigite os valores de y tabelados:\n\n\t{ ");
                for ( int i = 0; i < num; i++ ) 
                    scanf("%lf",&tabelados[1][i]);
				//

                AjusteExponencial(num, tabelados, &a0, &a1, ajustados, &coef);
                printf("}\n\tCoeficiente de determinacao entre os pontos tabelados e ajustados = %.4lf\n",coef);
                printf("\n\tResultados:\t coeficientes da eq. exponencial ajustada = %.4lf %.4lf", a0, a1);
                printf("\n\tValores Y ajustados:\n\n\t{ ");
                for ( int i = 0; i < num; i++ )
                    printf("%.4lf ", ajustados[i]);
                //
                
            }

	        break;
	        
		}

        if(input == 0){ // controle de saida / pos operacao
            printf("\n");
            break;
        }
        
        else{
            num = 0;
            printf("\n\tPressione qualquer tecla para continuar...\t");
            scanf("%d", &input);
            cls();
        }
        
    }

    // nao era pra chegar ate aqui...
    return 0;
    
}