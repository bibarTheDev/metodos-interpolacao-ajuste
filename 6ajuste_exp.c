void printVetor(int ordem, double matriz[MTXSIZE])
{
    for(int i = 0; i < ordem; i++){
        printf("%4.4lf ", matriz[i]);
    }
    printf("\n");
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
void AjusteExponencial(int npts, double pts[2][MTXSIZE], double *a, double *b, double y[MTXSIZE], double coefDet[MTXSIZE])
{
    double coef, a1, a0, yhat[MTXSIZE];
    double sumX = 0, sumXQ = 0, sumY = 0, sumXY = 0;

    for(int i = 0; i < npts; i++){
        pts[1][i] = log(pts[1][i]);

        sumX += pts[0][i];
        sumXQ += pow(pts[0][i], 2);
        sumY += pts[1][i];
        sumXY += pts[0][i] * pts[1][i];        
    }  

    if(DEBUG) printf("%lf | %lf | %lf | %lf\n", sumX, sumXQ, sumY, sumXY);

    a1 = ((npts * sumXY) - (sumX * sumY)) / ((npts * sumXQ) - pow(sumX, 2)); 
    a0 = (sumY - (a1 * sumX)) / npts; 

    if(DEBUG) printf("%lf | %lf\n", a0, a1);

    *a = exp(a0);
    *b = exp(a1);

    if(DEBUG) printf("%lf | %lf\n", *a, *b);

    for(int i = 0; i < npts; i++){
        y[i] = (*a) * pow(*b, pts[0][i]);
        yhat[i] = a0 + a1 *  pts[0][i];
    }  

    *coefDet = CoefDeterminacao(npts, pts, yhat);

    if(DEBUG) printVetor(npts, y);
    if(DEBUG) printf("%lf\n", *coefDet);

    return;
}