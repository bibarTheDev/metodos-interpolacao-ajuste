void printVetor(int ordem, double matriz[MTXSIZE])
{
    for(int i = 0; i < ordem; i++){
        printf("%4.4lf ", matriz[i]);
    }
    printf("\n");
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
void AjusteReta(int npts, double pts[2][MTXSIZE], double *a0, double *a1, double y[MTXSIZE], double coefDet[MTXSIZE])
{
    double coef;
    double sumX = 0, sumXQ = 0, sumY = 0, sumXY = 0;

    for(int i = 0; i < npts; i++){
        sumX += pts[0][i];
        sumXQ += pow(pts[0][i], 2);
        sumY += pts[1][i];
        sumXY += pts[0][i] * pts[1][i];        
    }  

    if(DEBUG) printf("%lf | %lf | %lf | %lf\n", sumX, sumXQ, sumY, sumXY);

    *a1 = ((npts * sumXY) - (sumX * sumY)) / ((npts * sumXQ) - pow(sumX, 2)); 
    *a0 = (sumY - (*a1 * sumX)) / npts; 

    if(DEBUG) printf("%lf | %lf\n", *a0, *a1);

    for(int i = 0; i < npts; i++){
        y[i] = *a0 + (*a1 * pts[0][i]);
    }  

    coef = CoefDeterminacao(npts, pts, y);

    if(DEBUG) printVetor(npts, y);
    if(DEBUG) printf("%lf\n", coef);

    return;
}