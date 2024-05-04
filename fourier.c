#include<stdlib.h>
#include<stdio.h>
#include<string.h> 
#include<math.h>
#include<time.h>

#define neg2PI -3.1415926535*2

char openDirPath[256];

/// x + iy
typedef struct _complex COMPLEX;

void compCpy(COMPLEX* a, COMPLEX b){
    // null check
    if(a == NULL){
        return;
    }

    a->x = b.x;
    a->y = b.y;
}
COMPLEX compAdd(COMPLEX a, COMPLEX b){
    COMPLEX out;
    
    out.x = a.x + b.x;
    out.y = a.y + b.y;

    return out;
}

COMPLEX compSub(COMPLEX a, COMPLEX b){
    COMPLEX out;

    out.x = a.x - b.x;
    out.y = a.y - b.y;

    return out;
}

COMPLEX compMalF(double a, COMPLEX b){
    COMPLEX out;

    out.x = a * b.x;
    out.y = a * b.y;

    return out;
}

 COMPLEX compMal(COMPLEX a, COMPLEX b){
    COMPLEX out;

    out.x = a.x * b.x - a.y * b.y;
    out.y = a.x * b.y + a.y * b.x;

    return out;
}

COMPLEX compExp(double x){
    COMPLEX out;

    out.x = cos(x);
    out.y = sin(x);

    return out;
}

double imabs(COMPLEX x){
    double  a = x.x,
            b = x.y,
            out;

    out = sqrt(a*a + b*b);

    return out;
}

int checkPow2(int n){
    return (n>1)&&!(n&(n-1));
}

int invBit(int i, int n){
    int o = 0;

    for(; n > 0; n--){
        o = (o<<1) | (i&1);
        i >>= 1;
    }

    return o;
}

void replaceIdx(unsigned int **idx, unsigned int n, short *bitN){
    int t = n;

    *idx = (unsigned int *)malloc(sizeof(unsigned int)*n);

    *bitN = -1;
    while(t > 0){
        t >>= 1;
        (*bitN)++;
    }

    for(int i = 0; i < n; i++){
        (*idx)[i] = invBit(i, (int)*bitN);
    }
}

COMPLEX butterfly(COMPLEX a, COMPLEX b, int n, int k){
    return compAdd(
        a,
        compMal(
            b,
            compExp(neg2PI/n * k)
        )
    );
}

int main(){
//* debug
clock_t startTime = clock();
//*/
    FILE* fpInput;
    FILE* fpOutput;
    double SampleInterval;
    float buf;
    short bitN;
    unsigned int *idx = NULL;
    unsigned int fN = 1;
    unsigned int n = 1;

///// read params /////
    fopen_s(&fpInput, "params.txt", "r");
    
    if(fpInput == NULL){
        printf_s("params.txt can't open.");
        return -1;
    }
    
    fscanf_s(fpInput, "%lf", &SampleInterval);
    
    fclose(fpInput);

///// read data from the input file /////
    // number of data
    fopen_s(&fpInput, "input.txt", "r");

    if(fpInput == NULL){
        printf_s("input.txt can't open.");
        return -1;
    }

    fscanf_s(fpInput, "%f", &buf);
    while(fgetc(fpInput) != '\n'){
        fscanf_s(fpInput, "%f", &buf);
        fN++;
    }
    {
        char str[256];
        while(fgets(str, 256, fpInput) != NULL){
            n++;
        }
    }
    printf("%d", n);
    fclose(fpInput);

    // check power of two datas
    if(!checkPow2(n)){
        printf_s("number of input data is not power of two.");
    }

    SampleInterval = 1 / SampleInterval;
    SampleInterval = SampleInterval / n;


    // read data
    fopen_s(&fpInput, "input.txt", "r");

    if(fpInput == NULL){
        printf_s("input.txt can't open.");
        return -1;
    }

    COMPLEX* f = (COMPLEX*)malloc(sizeof(COMPLEX)*fN*n);
    for(int i = 0; i < n; i++){
        for(int j = 0; j < fN; j++){
            fscanf(fpInput, "%f", &buf);

            f[i*fN+j].x = buf;
            f[i*fN+j].y = 0;
        }
    }
    fclose(fpInput);

    // bit inversion
    replaceIdx(&idx, n, &bitN);

    int s = 1,
        e = n / 2;
    e = e > n ? n : e;

///// FFT /////
    for(int i = 0; i < fN; i++){
        for(int j = 0; j < bitN-1; j++){
            int pearN = n >> (j+1), // number of pear
                deg = 1 << j;   // degree 2^j
            for(int k = 0; k < pearN; k++){
                int s = k*deg<<1;
                for(int l = 0; l < deg; l++){ // (deg+1)th degree DFT
                    COMPLEX a, b;

                    compCpy(&a, f[idx[s+l]*fN + i]);
                    compCpy(&b, f[idx[s+l+deg]*fN + i]);

                    compCpy(
                        &f[idx[s+l]*fN + i],
                        butterfly(
                            a,
                            b,
                            deg<<1,
                            l
                        )
                    );

                    compCpy(
                        &f[idx[s+l+deg]*fN + i],
                        butterfly(
                            a,
                            b,
                            deg<<1,
                            l + deg
                        )
                    );
                }
            }
        }

        int pearN = n >> (bitN); // number of pear
        int deg = 1 << bitN-1;   // degree 2^i
        for(int k = 0; k < pearN; k++){
            int s = k*deg<<1;
            for(int l = 0; l < deg; l++){ // (deg+1)th degree DFT
                COMPLEX a, b;

                compCpy(&a, f[idx[s+l]*fN + i]);
                compCpy(&b, f[idx[s+l+deg]*fN + i]);

                compCpy(
                    &f[idx[s+l]*fN + i],
                    butterfly(
                        a,
                        b,
                        bitN-1,
                        l
                    )
                );
            }
        }
    }

///// output to the output file /////
    // output complex
    fopen_s(&fpOutput, "output1.txt", "w");
    for(int i = s; i < e; i++){
        char *str;
        int bytes = 309+(621)*fN + 2; // max double e+308 + (e+308 * 2 + ('i') + ('+') + (' ')) * fN + ('\n') + ('\0')

        str = (char*)malloc(sizeof(char)*bytes); // max double e+308 + (e+308 * 2 + ('i') + ('+') + (' ')) * fN + ('\n') + ('\0')

        sprintf_s(str, bytes, "%lf", SampleInterval * i);
        for(int j = 0; j < fN; j++){
            sprintf_s(str, bytes, "%s %lf+%lfi", str, f[idx[i]*fN + j].x, f[idx[i]*fN + j].y);
        }
        sprintf_s(str, bytes, "%s\n", str);

        fprintf_s(fpOutput, "%s", str);
        free(str);
    }
    fclose(fpOutput);

    // output imabs
    fopen_s(&fpOutput, "output2.txt", "w");
    for(int i = s; i < e; i++){
        char *str;
        int bytes = 309+310*fN + 2; // max double e+308 + (e+308 + (' ')) * fN + ('\n') + ('\0')

        str = (char*)malloc(sizeof(char)*bytes);

        sprintf_s(str, bytes, "%lf", SampleInterval * i);
        for(int j = 0; j < fN; j++){
            sprintf_s(str, bytes, "%s %lf", str, imabs(f[idx[i]*fN + j]));
        }
        sprintf_s(str, bytes, "%s\n", str);

        fprintf_s(fpOutput, "%s", str);
        free(str);
    }
    fclose(fpOutput);
    
    free(f);
    free(idx);

//* debug
clock_t endTime = clock();
printf("\n%ld", endTime - startTime);
//*/

    return 0;
}