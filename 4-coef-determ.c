/**
 * Retorna o coeficiente de determinação entre os pontos
 * tabelados e os pontos ajustados. 
 * <arg-1>: inteiro, número de pontos tabelados (E) 
 * <arg-2>: tabela, pontos tabelados (E) 
 * <arg-3>: vetorY, valores Y ajustado (E)
*/
double CoefDeterminacao(int npts, double pts[MTXSIZE], double y[MTXSIZE])
{
    double r;
    double sumY = 0, sumYQ = 0, sumErrQ = 0;

    for(int i = 0; i < npts; i++){
        sumY += pts[i];
        sumYQ += pow(pts[i], 2);
        sumErrQ += pow(pts[i] - y[i], 2);
    }
    
    printf("%lf | %lf | %lf\n", sumY, sumYQ, sumErrQ);

    r = 1 - ((npts * sumErrQ) / ((npts * sumYQ) - pow(sumY, 2)));

    return r;
}