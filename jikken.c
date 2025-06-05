#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<math.h>

#define BUFSIZE 1024 //ファイルから読み込む一行の最大文字数
#define MAX_SEQ_NUM 30 //一つの転写因子に対して与えられる結合部位配列の最大数
#define MAX_GENE_NUM 8 /*与えられるプロモータ領域の最大遺伝子数*/

char g_motif[MAX_SEQ_NUM][BUFSIZE]; //転写因子の結合部位配列を保存する配列

struct promoter{
  char name[BUFSIZE];
  char seq[BUFSIZE];
}g_pro[MAX_GENE_NUM]; //遺伝子のプロモータ領域を保存する構造体

//グローバル変数はローカル変数と区別するため、全部大文字にするかg_を先頭につけるのが一般的

int read_multi_seq(char* filename){
  int seq_num = 0;
  char buffer[BUFSIZE];
  FILE *fp = fopen(filename,"r");

  if(fp == NULL){
    printf("motif_region_file open error.\n");
    exit(1); //ファイルが開けなかった場合プログラムを終了
  }

  while(fscanf(fp, "%s", buffer) != EOF){ //プログラムから一行ずつ読み込む
    if(buffer[strlen(buffer)-1]=='\n'){
      buffer[strlen(buffer)-1]='\0'; //改行を切り落とす
    }
    strcpy(g_motif[seq_num],buffer); //結合部位配列を保存
    seq_num++;
  }
  return seq_num;
}

int read_promoter(char *filename){
  int gene_num = 0;  
  char buffer[BUFSIZE];
  FILE *fp = fopen(filename,"r");

  if(fp == NULL){
    printf("scorefile open error.\n");
    exit(1);
  }

  while(fscanf(fp, "%s", buffer) != EOF){
    if(buffer[strlen(buffer)-1]=='\n'){
      buffer[strlen(buffer)-1]='\0';
    }
    
    if(buffer[0]=='>'){
      strcpy(g_pro[gene_num].name,buffer+1); 
    }else{
      strcpy(g_pro[gene_num].seq,buffer);
      gene_num++;
    }    
  }
  return gene_num;
}

int main(int argc, char* argv[]){
  int seq_num = read_multi_seq(argv[1]); //１番目の引数で指定した転写因子の複数の結合部位配列を読み込む

  printf("motif region:\n");
  for(int i = 0; i < seq_num; i++){
    printf("%s\n",g_motif[i]); //読み込んだ転写因子の結合部位配列を表示
  }
  printf("\n");

  int gene_num = read_promoter(argv[2]);  //２番目の引数で指定した遺伝子のプロモータ領域を読み込む
  
  printf("promoter_sequence:\n");
  for(int i = 0; i < gene_num; i++){
    printf(">%s\n", g_pro[i].name); //読み込んだプロモータ領域を表示
    printf("%s\n", g_pro[i].seq);
  }

  int FLQ[4][10] = {0};
  int i, j, s, t = 0;
  
      for(s=0; s<seq_num; s++){
        for(t=0; t<10; t++){
          
      if (g_motif[s][t]=='A') {
        FLQ[0][t] += 1;
      }
      if (g_motif[s][t]=='C') {
        FLQ[1][t] += 1;
      }
      if (g_motif[s][t]=='G') {
        FLQ[2][t] += 1;
      }
      if (g_motif[s][t]=='T') {
        FLQ[3][t] += 1;
      }
        }
      }

for (i=0; i<4; i++){
    for (j=0; j<10; j++){
    printf(" %d", FLQ[i][j]);
  }
  printf("\n");
}

i, j = 0;
  
      for (i=0; i<4; i++){
        for(j=0; j<10; j++){
          FLQ[i][j] += 1;
        }
      }

int row_sum[10] = {0};

      for (j=0; j<10; j++){
        for(i=0; i<4; i++){
           row_sum[j] += FLQ[i][j];
        }
      }

float p[4][10] = {0};

      for (i=0; i<4; i++){
        for(j=0; j<10; j++){
          p[i][j] = FLQ[i][j]/(seq_num+4.0);
        }
      }

float flq_sum = 2.0*(7519429+4637676);
float q[4] = {0};
q[0] = 7519429/flq_sum;
q[1] = 4637676/flq_sum;  
q[2] = 4637676/flq_sum;
q[3] = 7519429/flq_sum; 

i, j = 0;
  
double OddFLQ[4][10] = {0};
for(i=0; i<4; i++){
  for(j=0; j<10; j++){
     OddFLQ[i][j] = log10(p[i][j]/q[i]);
  }
}

printf("p=\n");
    
    for (i=0; i<4; i++){
    for (j=0; j<10; j++){
    printf(" %f", p[i][j]);
  }
  printf("\n");
}
for (i=0; i<4; i++){
    for (j=0; j<10; j++){
    printf(" %f", OddFLQ[i][j]);
  }
  printf("\n");
}

  printf("Motif: MATa1\n\n");  // モチーフ名

  const int motif_len = 10;
  double threshold = 2.0;

  for (int g = 0; g < gene_num; g++) {
    char* prom = g_pro[g].seq;
    int prom_len = strlen(prom);

    for (int i = 0; i <= prom_len - motif_len; i++) {
      double score = 0.0;
      int valid = 1;

      for (int j = 0; j < motif_len; j++) {
        char base = prom[i + j];
        int row;
        if (base == 'A') row = 0;
        else if (base == 'C') row = 1;
        else if (base == 'G') row = 2;
        else if (base == 'T') row = 3;
        else {
          valid = 0;
          break;
        }
        score += OddFLQ[row][j];
      }

      if (valid && score >= threshold) {
        printf("pro:%s\n", g_pro[g].name);
        printf("pos:%d\n", i);
        printf("hit(");
        for (int j = 0; j < motif_len; j++) {
          printf("%c", prom[i + j]);
        }
        printf(")= %.2f\n\n", score);
      }
    }
  }
  return 0;
}