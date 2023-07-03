
void extraiCol(int ordem, double matriz[MTXSIZE][MTXSIZE], int col, double vetor[MTXSIZE])
{
    for(int i = 0; i < ordem; i++){
        vetor[i] = matriz[i][col];
    }
}

//////////////////////

void copiaMatriz(double matriz[MTXSIZE][MTXSIZE], double saida[MTXSIZE][MTXSIZE])
{
    memcpy(saida, matriz, MTXSIZE * MTXSIZE * sizeof(double));
    return;
}

//////////////////////

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

//////////////////////

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

//////////////////////

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

//////////////////////

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

//////////////////////

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

//////////////////////

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

//////////////////////

void AjustePolinomio(int npts, int grau, double tabela[2][MTXSIZE], double vetorA[MTXSIZE], double vetorY[MTXSIZE], double& coefDeterminacao) {
    
    double matriz[MTXSIZE][MTXSIZE];
    double vetorI[MTXSIZE];

    // Gera matriz e vetor independentedo sistema
    for (int i = 0; i <= grau; i++) {

        for (int j = 0; j <= grau; j++) {
            double sum = 0;
            for (int k = 0; k < npts; k++) 
                sum += pow(tabela[0][k], i + j);
            //
            matriz[i][j] = sum;
        }

        double sum = 0;
        for (int k = 0; k < npts; k++)
            sum += tabela[1][k] * pow(tabela[0][k], i);
        //
        vetorI[i] = sum;

    }

    DecomposicaoLU(grau + 1, matriz, vetorI, vetorA);

    // Calcula Y ajustados
    for (int k = 0; k < npts; k++) {
        vetorY[k] = 0;
        for (int i = 0; i <= grau; i++)
            vetorY[k] += vetorA[i] * pow(tabela[0][k], i);
        //
    }

    coefDeterminacao = CoefDeterminacao(npts, tabela, vetorY);

}

////////////////////