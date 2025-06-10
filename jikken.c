#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h> 

#define BUFSIZE 1024 //ファイルから読み込む一行の最大文字数
#define MAX_SEQ_NUM 30 //一つの転写因子に対して与えられる結合部位配列の最大数
#define MAX_GENE_NUM 8 /*与えられるプロモータ領域の最大遺伝子数*/
// 全プロモーター配列の長さ合計をMAX_GENE_NUM * BUFSIZE と仮定し、
// それをモチーフ長で割った最大スコア数。安全のため少し大きめに設定。
#define MAX_TOTAL_SCORES 10000


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
    fclose(fp); // ファイルを閉じる
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
    fclose(fp); // ファイルを閉じる
    return gene_num;
}

int main(int argc, char* argv[]){
    // コマンドライン引数の数をチェック
    if (argc != 3) {
        printf("使用法: %s <モチーフファイル> <プロモーターファイル>\n", argv[0]);
        return 1;
    }

    int seq_num = read_multi_seq(argv[1]); //１番目の引数で指定した転写因子の複数の結合部位配列を読み込む

    printf("motif region:\n");
    // モチーフの長さをここで取得する
    int motif_len = 0; 
    if (seq_num > 0) {
        motif_len = strlen(g_motif[0]); // 読み込んだ最初のモチーフの長さを取得
        // 全てのモチーフが同じ長さであることを前提とします。
        // 異なる長さのモチーフが混在する場合は、別途処理が必要ですが、
        // この課題では統一されていると仮定します。
        printf("検出されたモチーフの長さ: %d\n\n", motif_len); // 追加: 検出されたモチーフ長を表示
        for(int i = 0; i < seq_num; i++){
            printf("%s\n",g_motif[i]); //読み込んだ転写因子の結合部位配列を表示
        }
    } else {
        printf("エラー: モチーフ配列が読み込まれませんでした。\n");
        return 1;
    }
    printf("\n");

    int gene_num = read_promoter(argv[2]);  //２番目の引数で指定した遺伝子のプロモータ領域を読み込む
    
    printf("promoter_sequence:\n");
    for(int i = 0; i < gene_num; i++){
        printf(">%s\n", g_pro[i].name); //読み込んだプロモータ領域を表示
        printf("%s\n", g_pro[i].seq);
    }
    printf("\n");

    // FLQのサイズをmotif_lenに合わせる（動的配列ではないので、BUFSIZEなど十分大きいサイズで定義）
    // モチーフの最大長はBUFSIZEなので、FLQもBUFSIZEで定義します。
    int FLQ[4][BUFSIZE] = {0}; // BUFSIZEは1024なので十分なサイズ

    int i, j, s, t = 0; // ループ変数i,jは、後で再利用されるため、初期化をここでしない。
    
    // 頻度表の計算
    for(s=0; s<seq_num; s++){
        for(t=0; t<motif_len; t++){ // ここを motif_len に変更
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

    // 頻度表の出力
    printf("--- 頻度表 ---\n"); // 追加：見出し
    printf("  "); // 列番号のヘッダー
    for (j=0; j<motif_len; j++){ // ここを motif_len に変更
        printf("%d ", j + 1);
    }
    printf("\n");

    for (i=0; i<4; i++){
        for (j=0; j<motif_len; j++){ // ここを motif_len に変更
            printf(" %d", FLQ[i][j]);
        }
        printf("\n\n");
    }

    // 擬似頻度1の加算
    for (i=0; i<4; i++){
        for(j=0; j<motif_len; j++){ // ここを motif_len に変更
            FLQ[i][j] += 1;
        }
    }

    int row_sum[BUFSIZE] = {0}; // サイズをBUFSIZEに
    for (j=0; j<motif_len; j++){ // ここを motif_len に変更
        for(i=0; i<4; i++){
            row_sum[j] += FLQ[i][j];
        }
    }

    float p[4][BUFSIZE] = {0}; // サイズをBUFSIZEに

    // p_i(x|M) の計算 (擬似頻度を加味した確率)
    for (i=0; i<4; i++){
        for(j=0; j<motif_len; j++){ // ここを motif_len に変更
            p[i][j] = (float)FLQ[i][j]/(seq_num+4.0); // FLQ[i][j] を float にキャストして正確な浮動小数点除算
        }
    }

    float total_background_sum = 7519429.0f + 4637676.0f + 4637676.0f + 7519429.0f; // floatリテラルを使用
    float q[4] = {0};
    q[0] = 7519429.0f/total_background_sum; // Aのバックグラウンド確率
    q[1] = 4637676.0f/total_background_sum; // Cのバックグラウンド確率
    q[2] = 4637676.0f/total_background_sum; // Gのバックグラウンド確率
    q[3] = 7519429.0f/total_background_sum; // Tのバックグラウンド確率

    double OddFLQ[4][BUFSIZE] = {0}; // サイズをBUFSIZEに
    for(i=0; i<4; i++){
        for(j=0; j<motif_len; j++){ // ここを motif_len に変更
            OddFLQ[i][j] = log10(p[i][j]/q[i]);
        }
    }
    printf("p=\n"); // p行列の出力
    for (i=0; i<4; i++){
        for (j=0; j<motif_len; j++){ // ここを motif_len に変更
            printf(" %f", p[i][j]);
        }
        printf("\n");
    }
    printf("\n"); // 追加：改行

    printf("--- 対数オッズスコア行列 ---\n"); // 追加：見出し
    printf("  "); // 列番号のヘッダー
    for (j=0; j<motif_len; j++){ // ここを motif_len に変更
        printf("%d ", j + 1);
    }
    printf("\n");
    for (i=0; i<4; i++){
        for (j=0; j<motif_len; j++){ // ここを motif_len に変更
            printf(" %f", OddFLQ[i][j]);
        }
        printf("\n");
    }
    printf("\n"); // 追加：改行

    printf("Motif: MATa1\n\n");  // モチーフ名

    double threshold = 2.0; // 閾値

    // スコア分布計算のための変数
    double all_scores[MAX_TOTAL_SCORES];
    int score_count = 0;
    double sum_scores = 0.0;
    double sum_sq_diff_scores = 0.0; // 差の二乗の合計
    
    // 結合部位の探索
    for (int g = 0; g < gene_num; g++) {
        char* prom = g_pro[g].seq;
        int prom_len = strlen(prom);

        for (int i_prom_pos = 0; i_prom_pos <= prom_len - motif_len; i_prom_pos++) { // ループ変数名を変更して重複を避ける
            double score = 0.0;
            int valid = 1;

            for (int j_motif_pos = 0; j_motif_pos < motif_len; j_motif_pos++) { // ループ変数名を変更して重複を避ける
                char base = prom[i_prom_pos + j_motif_pos];
                int row;
                if (base == 'A') row = 0;
                else if (base == 'C') row = 1;
                else if (base == 'G') row = 2;
                else if (base == 'T') row = 3;
                else {
                    valid = 0;
                    break;
                }
                score += OddFLQ[row][j_motif_pos];
            }

            // 有効なスコアであれば全て配列に格納 (閾値にかかわらず)
            if (valid) {
                if (score_count < MAX_TOTAL_SCORES) { // 配列の範囲チェック
                    all_scores[score_count] = score;
                    sum_scores += score; // 平均計算のためにスコアを合計
                    score_count++;
                } else {
                    fprintf(stderr, "警告: スコア配列が上限に達しました。一部のスコアは記録されません。\n");
                    // このプロモーターの残りの部分をスキップするか、プログラムを終了するかは要検討
                    // ここでは、警告を出して処理を続けます
                    break; 
                }
            }

            // 閾値を超えた結合部位は出力
            if (valid && score >= threshold) {
                printf("pro:%s\n", g_pro[g].name);
                printf("pos:%d\n", i_prom_pos);
                printf("hit(");
                for (int j_print = 0; j_print < motif_len; j_print++) {
                    printf("%c", prom[i_prom_pos + j_print]);
                }
                printf(")= %.2f\n\n", score);
            }
        }
    }

    // スコア分布の平均と標準偏差の計算と表示
    printf("--- スコア分布の統計情報 ---\n");
    if (score_count > 0) {
        double mean = sum_scores / score_count; // 平均の計算
        printf("平均 (Mean): %.4f\n", mean);

        // 標準偏差の計算
        for (int k = 0; k < score_count; k++) {
            sum_sq_diff_scores += (all_scores[k] - mean) * (all_scores[k] - mean);
        }
        // 分散の計算 (不偏分散ではなく、標本分散を使う)
        double variance = sum_sq_diff_scores / score_count; 
        double std_dev = sqrt(variance); // 標準偏差
        printf("標準偏差 (Standard Deviation): %.4f\n", std_dev);
    } else {
        printf("有効なスコアがありませんでした。\n");
    }
    printf("\n");

    return 0;
}