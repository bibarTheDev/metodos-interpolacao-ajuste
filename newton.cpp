#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MTXSIZE 105
#define DEBUG true
double Determinante (int ordem, double matriz[MTXSIZE][MTXSIZE]);

// FUNCOES AUXILIARES

/**
 * cls
*/
void cls()
{
    // limpa a tela
    printf("\e[1;1H\e[2J");
}

// ROTINA NEWTON

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
float Newton (int n, double x[], double y[], double point) {

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

// FLUXO DE PROGRAMA

int main()
{

    int input = -1, num = 0;
    double x[MTXSIZE]; double y[MTXSIZE]; double point = 0;
    
    double precisao = 0;

    while(1){

        /* menu */ {
            printf("\n\tMETODOS DE SISTEMAS LINEARES\n\n");
            printf("\t01 - Interpolacao Metodo de Newton\n");
            printf("\t02 - Interpolacao Metodo de Newton-Gregory\n");
            printf("\t03 - \n");
            printf("\t04 - \n");
            printf("\t05 - \n");
            printf("\t06 - \n");
            printf("\t00 - \n\n");
        }

        do{ // ler input
            printf("\tDigite uma opcao:\t");
            scanf("%d", &input);
        }while(input < 0 || input > 10);

        switch(input){ // fluxo para cada metodo
		    
	        case 1:{ // Newton

	            cls();

                printf("\n\tINTERPOLACAO METODO DE NEWTON\n\n");
				printf("\tDigite o numero de pontos tabelados (Max: 100):\t"); scanf("%d", &num);

                printf("\n\tDigite o vetor dos valores de x:\n\n\t{ ");
                for ( int i = 0; i < num; i++ ) {
                    scanf("%lf",&x[i]);
				}

                printf("\n\tDigite o vetor dos valores de y:\n\n\t{ ");
                for ( int i = 0; i < num; i++ ) {
                    scanf("%lf",&y[i]);
				}

                printf("\n\tDigite o ponto a ser interpolado:\t"); scanf("%lf", &point);

                printf("\n\tResultado:\t%.4lf\n", Newton(num, x, y, point));
                
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