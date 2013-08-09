
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <sched.h>

/*
 * 
 */
#define MAXFILENAME 500
#define MAXCOMMAND 1000
#define MAXNAMES 500
#define COMMENT_SIZE 1000
#define GAP_LIST "-.#*~"

#define AMINOACID_NUMBER 20
char *amino_acid_order = "ARNDCQEGHILKMFPSTWYV";
double amino_acid_values[] = {0.307, 0.106, 0.240, 0.479, 0.165,  0.248,  0.449,  0.265, 0.202, 0.198,  0.246,  0.111,  0.265,  0.246, 0.327, 0.281, 0.242, 0.172, 0.185, 0.200};

struct Parameters {
    
    int debug;
    int noplot;
    char *gnuplot;
    char *input_file;
    int window_size;
    float threshold;
    char *graph_file;
    char *result_file;
    char *data_file;
    char *prompt;
};
typedef struct Parameters Parameters;

struct Sequence {
    char **seq;
    char **seq_names;
    int *seq_len;
    int *wseq_len;
    int *seq_id;
    int nseq;
    int min_len;
    int max_len;
    
};
typedef struct Sequence Sequence;


Parameters* init();
void free_parameters(Parameters *par);
Parameters* load_parameters(int argc, char **argv);
void help();

Sequence *declare_sequence(int nseq, int max_len, int min_len);
void free_sequence(Sequence *LS);
Sequence *read_fasta_sequences(char *fname, int window_size);
int count_nseqs(char *fname, char x);
int count_n_char_in_file(char *fname);
int fscanf_seq_name ( FILE *fp, char *sname);
int is_gap(char x);
int is_in_set(char r, char *list);

void sliding(Sequence *S, Parameters *par, int cont);
void gnuplot_sld(Parameters *par);

float pnormal(double x);
float return_value(char residue, float threshold);

int main(int argc, char** argv) {

    Parameters *par;
    Sequence *S;
    FILE *fp;
    
    int cont=0;
    
    par = load_parameters(argc, argv);
    S = read_fasta_sequences(par->input_file, par->window_size);
    
    for(cont=0; cont<S->nseq; cont++){
        printf("%s - %s\n", S->seq_names[0], S->seq[0]);
        if(S->nseq > 1){
            sprintf(par->data_file, "data.%d.txt", cont);
            if(par->noplot==0) {
                sprintf(par->graph_file, "plot.%d.png", cont);
            } 
        }
        sliding(S, par, cont);
        if(par->noplot==0) {
            gnuplot_sld(par);
        }
        if(S->nseq > 1){
            remove(par->data_file);
        }
    }

    free_sequence(S);
    free_parameters(par);
    
    return (EXIT_SUCCESS);
}

Parameters* init(){
    
    Parameters *par;
    
    par = calloc(1, sizeof(Parameters));
    par->debug = 0;
    par->noplot = 0;
    par->gnuplot = (char *) calloc(FILENAME_MAX, sizeof(char));
    par->input_file = (char *) calloc(FILENAME_MAX, sizeof(char));
    sprintf(par->input_file, "protein.txt");
    par->threshold = 0.225;
    par->window_size = 7;
    
    par->graph_file = (char *) calloc(FILENAME_MAX, sizeof(char));
    sprintf(par->graph_file, " ");

    par->result_file = (char *) calloc(FILENAME_MAX, sizeof(char));
    sprintf(par->result_file, "results.%d", getpid());
    
    par->data_file = (char *) calloc(FILENAME_MAX, sizeof(char));
    sprintf(par->data_file, "sliding.%d", getpid());
    
    par->prompt = (char *) calloc(FILENAME_MAX, sizeof(char));
    
    return par;
     
}

void free_parameters(Parameters *par){

    if (!par) return;
    
    free(par->input_file);
    free(par->gnuplot);
    free(par->result_file);
    free(par->graph_file);
    free(par->data_file);
    free(par->prompt);
    
    free(par);
}

Parameters* load_parameters(int argc, char **argv){
    
    int i=0, found=0;
    Parameters *par;
    
    if(argc < 2 && argc > 8){
        fprintf(stderr, "\nERROR - Wrong parameters number\n");
        exit(EXIT_FAILURE);  
    }
    
    if(strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "-H") == 0){
        help();
    }
    
    par = init();
        
    for(i=1; i<argc; i++){
        if(strcmp(argv[i], "-in") == 0){
            if(argv[i+1][0] != '-'){
                sprintf(par->input_file, argv[i+1]);
                i++;
                found=1;
            }
            else {
                fprintf(stderr, "\nERROR - Wrong input file name\n");
                exit(EXIT_FAILURE);  
            }
        }
        else if(strcmp(argv[i], "-t") == 0){
            if(argv[i+1][0] != '-'){
                par->threshold =atof(argv[i+1]);
                i++;          
            }
            else {
                fprintf(stderr, "\nERROR - Wrong threshold value\n");
                exit(EXIT_FAILURE);  
            }
        }
        else if(strcmp(argv[i], "-w") == 0){
            if(argv[i+1][0] != '-'){
                par->window_size =atoi(argv[i+1]);
                i++;          
            }
            else {
                fprintf(stderr, "\nERROR - Wrong window size value\n");
                exit(EXIT_FAILURE);  
            }
        }
        else if(strcmp(argv[i], "-rf") == 0){
            if(argv[i+1][0] != '-'){
                sprintf(par->result_file, argv[i+1]);
                i++;
            }
            else {
                fprintf(stderr, "\nERROR - Wrong result file name\n");
                exit(EXIT_FAILURE);  
            }
        }
        else if(strcmp(argv[i], "-df") == 0){
            if(argv[i+1][0] != '-'){
                sprintf(par->data_file, argv[i+1]);
                i++;
            }
            else {
                fprintf(stderr, "\nERROR - Wrong data file name\n");
                exit(EXIT_FAILURE);  
            }
        }
        else if(strcmp(argv[i], "-gf") == 0){
            if(argv[i+1][0] != '-'){
                sprintf(par->graph_file, argv[i+1]);
                i++;
            }
            else {
                fprintf(stderr, "\nERROR - Wrong graph file name\n");
                exit(EXIT_FAILURE);  
            }
        }
        else if(strcmp(argv[i], "-debug") == 0){
            par->debug = 1;
        }
        else if(strcmp(argv[i], "-noplot") == 0){
            par->noplot = 1;
        }
        else {
            fprintf(stderr, "\nERROR - Wrong parameter\n");
            exit(EXIT_FAILURE);  
        }
    }
    
    if(found == 0){
        fprintf(stderr, "\nERROR - No input file\n");
        exit(EXIT_FAILURE);  
    }
    
    sprintf(par->gnuplot, "gnuplot");
 /*  if(par->noplot == 0){
        par->gnuplot = getenv ("GNUPLOT_BIN");
        if (par->gnuplot == NULL){
            printf("%s\n", par->gnuplot);
            
            printf("%s\n", par->gnuplot);

            fprintf(stderr, "\nERROR - Gnuplot not found. Option -noplot = 1\n");
            fprintf(stderr, "Write the path of the gnuplot bin in your environment file using the variable GNUPLOT_BIN=path_to_gnuplot\n");
            par->noplot = 1;

        
        
    }
       }*/ 
    
    if(par->debug != 0){
        fprintf(stdout, "\nDEBUG:\n");
        fprintf(stdout, "Window size      : %f\n", par->window_size);
        fprintf(stdout, "Threshold value  : %f\n", par->threshold);
        fprintf(stdout, "Input file name  : %s\n", par->input_file);
        fprintf(stdout, "Graph file name  : %s\n", par->graph_file);;
        fprintf(stdout, "Result file name : %s\n", par->result_file);
        fprintf(stdout, "Data file name   : %s\n", par->data_file);
        fprintf(stdout, "Gnuplot bin      : %s\n", par->gnuplot);
    }

    return par;    
}

void help(){
    
    fprintf(stdout, "\nUsage: ampaC -in <input fasta file> [other options]\n\n");
    fprintf(stdout, "Available options:\n");
    fprintf(stdout, "-in      Fasta sequence input file\n");
    fprintf(stdout, "-t       Threshold value (default: 7)\n");
    fprintf(stdout, "-w       Window size value (default: 0.225)\n");
    fprintf(stdout, "-rf      File where store the program result\n");
    fprintf(stdout, "-df      File where store the produced plot data (Only Nseq == 1)\n");
    fprintf(stdout, "-gf      File where store the generated plot (Only Nseq == 1)\n");
    fprintf(stdout, "-noplot  Skip plot creation step\n");
    fprintf(stdout, "-help    This help information\n");
 
    
    exit(EXIT_SUCCESS);      
}

Sequence *declare_sequence(int nseq, int max_len, int min_len){
    
    Sequence *LS;
    int i=0;
   
    LS=calloc(1, sizeof(Sequence));
    
    LS->seq = (char **) calloc(nseq, sizeof (char *));
    LS->seq_names = (char **) calloc(nseq, sizeof (char *));
    for(i=0; i<nseq; i++){
        LS->seq[i] = (char *) calloc(max_len+1, sizeof (char));
        LS->seq_names[i] = (char *) calloc(MAXNAMES+1, sizeof (char));
    }
    
    LS->seq_len= (int *) calloc(nseq, sizeof (int));
    LS->wseq_len= (int *) calloc(nseq, sizeof (int));
    LS->seq_id= (int *) calloc(nseq, sizeof (int));
    LS->max_len=max_len;
    LS->min_len=min_len;
    LS->nseq=nseq;

    return LS;
}

void free_sequence(Sequence *LS){

    if (!LS) return;
    int i=0;
    
    for(i=0; i<LS->nseq; i++){
        free(LS->seq[i]);
        free(LS->seq_names[i]);
    }
    free(LS->seq);
    free(LS->seq_names);
    
    free(LS->seq_len);
    free(LS->wseq_len);
    free(LS->seq_id);
    free(LS);
}

Sequence *read_fasta_sequences(char *fname, int window_size){

    FILE *fp;
    Sequence *S;
    int nseq=0, wlen=0, max_len_seq=0, min_len_seq=0, max=0, coor=0, current=0;
    int l=0, a=0, p=0, i=0, clen=0;
    int c, nul;
    int window;
    char *sub, *name;

    nseq=count_nseqs(fname, '>');

    if (nseq==0){
        return NULL;
    }
    void test(char *fname);
    min_len_seq=max=count_n_char_in_file(fname);
    sub = (char *) calloc(max+1, sizeof (char));
    name = (char *) calloc(10000, sizeof(char));
    if((fp = fopen(fname, "r")) == NULL){
        fprintf(stderr, "\nERROR - Cannot open the file: %s\n", fname);
         exit(EXIT_FAILURE);
    }
    c=fgetc(fp);
    while (c!=EOF){
        if (c=='>'){
            nul=fscanf_seq_name (fp, name);
            while((c=fgetc(fp))!='\n' && c!=EOF);
            while((c=fgetc(fp))!='>' && c!=EOF){
                if(isalnum(c) || is_gap(c)){
                    sub[clen++]=c;
                }
            }
            max_len_seq=(clen> max_len_seq)?clen: max_len_seq;
            min_len_seq=(clen< min_len_seq)?clen: min_len_seq;
            clen=0;
        }
        else {
            c=fgetc (fp);
        }
    }
    fclose (fp);
    free(name);
    free(sub);
    
    window = ((window_size + 1)/2)-1;
    wlen = (window * 2);
    max_len_seq += wlen;
    
    S = declare_sequence(nseq, max_len_seq, min_len_seq);   
 
    if((fp = fopen(fname, "r")) == NULL){
        fprintf(stderr, "\nERROR - Cannot open the file: %s\n", fname);
         exit(EXIT_FAILURE);
    }
    c=fgetc(fp);
    coor++;
    while (c!=EOF){
        if (c=='>'){
            coor+=fscanf_seq_name(fp, S->seq_names[current]);
            l=strlen (S->seq_names[current]);

            if (S->seq_names[current][l-1]==','|| S->seq_names[current][l-1]==';'){
                S->seq_names[current][l-1]='\0';
            }
            //seq_names[current]=translate_name(seq_names[current]);
            a=0;
            while ((c=fgetc(fp))!='\n' && c!=EOF && a<(COMMENT_SIZE-1)){
                coor++;
            }
            coor++;
            p=0;
            for(i=0; i<window; i++){
                S->seq[current][p] ='X';
                p++;
            }
            
            while ((c=fgetc(fp))!='>' && c!=EOF){
                coor++;
                if (!isspace(c)){
                    if(is_in_set(c, amino_acid_order)){
                        S->seq[current][p]=c;
                        p++;
                    }
                }
            }
            
            for(i=0; i<window; i++){
                S->seq[current][p] ='X';
                p++;
            }
            coor++;
            
            S->seq[current][p]='\0';
            S->seq_len[current]=p - wlen;
            S->wseq_len[current] = p;
            S->seq_id[current]=current;
            
            current++;
        }
        else {
            c=fgetc(fp);
            coor++;
        }
    }
    
    return S;
}

int count_nseqs(char *fname, char x){
    FILE *fp;
    int n, c;

    n=0;
    if((fp = fopen(fname, "r")) == NULL){
        fprintf(stderr, "\nERROR -  Cannot open the file: %s\n", fname);
        exit(EXIT_FAILURE);
    }
    while((c=fgetc(fp))!=EOF){
        n+=(c==x);
    }
    fclose (fp);
    return n;
}

int count_n_char_in_file(char *fname){
    int  c, n;
    FILE *fp;

    n=0;
    if((fp = fopen(fname, "r")) == NULL){
        fprintf(stderr, "\nERROR - Cannot open the file: %s\n", fname);
         exit(EXIT_FAILURE);
    }
    while((c=fgetc(fp))!=EOF){
        n++;
    }
    fclose (fp);
    return n;
}

int fscanf_seq_name ( FILE *fp, char *sname){

    static char *name;
    int r, nul;

    if (!sname){
        return 0;
    }
    if (!name){
        name = (char *) calloc(10000, sizeof (char));
    }
    nul=fscanf (fp, "%s", name);
    r=strlen (name);
    if(strlen(name)>MAXNAMES){
        fprintf(stderr, "\nWARNING: Seq Name Too long: [%s]. Truncated to %d", name, MAXNAMES);
    }
    name[MAXNAMES]='\0';
    sprintf(sname, "%s", name);

    return r;
}

int is_gap(char x){
    return (is_in_set( x, GAP_LIST));
}

int is_in_set(char r, char *list){

    char s[2];
    s[0]=r;

    if (strstr(list, s)!=NULL){
        return 1;
    }
    else {
        return 0;
    }
}

void sliding(Sequence *S, Parameters *par, int cont){
    
    FILE *kd, *rs,*tmp;
    
    int center, length, half;
    int i=0, j=0;
    int acc1=0;
    int acc2=0;
    int bacindex=0;
    float amvalue=0.0;
    float bacvalue=0.0;
    float PV = 0.0;
    float sum = 0.0;
    float prob=0.0;
    int rprob=0;
    int position=0;
    float APV=0.0;
    float manvalue=0.0;
    
    char *window, *residue;
    int init=0, last=0;
    
    if((kd = fopen(par->data_file, "w")) == NULL){
        fprintf(stderr, "\nERROR -  Cannot open the file: %s\n", par->data_file);
        exit(EXIT_FAILURE);
    }
    if((rs = fopen(par->result_file, "a")) == NULL){
        fprintf(stderr, "\nERROR -  Cannot open the file: %s\n", par->result_file);
        exit(EXIT_FAILURE);
    }
    
    fprintf(rs, "Protein: %s\n", S->seq_names[cont]);
    
    half = ((par->window_size+1)/2);
    
    fprintf(kd, "# The window size is currently set at %d.\n", par->window_size);
    fprintf(kd, "# Here are the AMSI values for this protein:\n");
    fprintf(kd, "# %s\n", S->seq_names[cont]);
    fprintf(kd, "# Pos\\ APV\n");
    fprintf(kd, "# ---\t-------------------\n");
    
    window = (char*) calloc(par->window_size +1, sizeof(char));
    residue = (char*) calloc(1, sizeof(char));
    
    for(i=0; i<=S->wseq_len[cont]-par->window_size; i++) {
        strncpy(window, S->seq[cont]+i, par->window_size);
        //window = substr($protein, $i, $window_size);
        sum=0;
        for(j=0; j<par->window_size; j++){
            PV = 0;
            strncpy(residue, window+j, 1);
            //my $residue = substr($window, $j, 1);
            
            if ((PV = return_value(residue[0], par->threshold)) == -1){
                fprintf(stderr, "ERROR - Wrong residue\n");
                PV=0.0;
            }
            sum+=PV;
        }
    
    
        center = i + half;
        position = center - 3;
        APV=(float) sum/par->window_size;
        fprintf(kd, "%d\t%.3f\n", position, APV);
        amvalue=amvalue+APV;

        if(acc1==3) {
            if(acc2>=10){
                init=position-acc2-2;
                last=position-1;
                bacvalue=((bacvalue)/(last-init+1));
               // prob = &Math::CDF::pnorm((bacvalue-0.2584)/0.02479);
                prob = pnormal((bacvalue-0.2584)/0.02479);
                prob*=100;
                //rbacvalue = (float) round(bacvalue * 1000) / 1000;                            
                rprob = round(prob);
                fprintf(rs, "Antimicrobial stretch found in %d to %d. Propensity value %.3f (%d%) \n", init, last, bacvalue, rprob);
                bacindex++;                                                        
                acc2=0;
                acc1=0;
                bacvalue=0.0;
                prob=0.0;
            }
            if(acc2<10){
                acc1=0;
                acc2=0;
                bacvalue=0.0;
            }
        }
        if(acc1!=3){
            if(acc2==0) {
                acc1=0;
                bacvalue=0.0;
            }                                                
            if(APV<=par->threshold) {
                acc2++;
                bacvalue+=APV;
            }
            if(APV>par->threshold){
                acc1++;
                bacvalue+=APV;
            }
        }
        

        if(position==(S->seq_len[cont]-((par->window_size)-1)) && acc2>=10) {
            init=position-acc2-2;
            last=position-1;
            bacvalue=((bacvalue)/(last-init+1));
            bacindex++;
            //prob = &Math::CDF::pnorm((bacvalue-0.2584)/0.02479);
            prob = pnormal((bacvalue-0.2584)/0.02479);
            prob*=100;
            //rbacvalue = nearest(.001, bacvalue);                            
            //rprob = nearest(1, prob);  
            rprob = round(prob);
            fprintf(rs, "Antimicrobial stretch found in %d to %d. Propensity value %.3f (%d%) \n", init, last, bacvalue, rprob);    
        }
    }

    fprintf(rs, "# This protein has %d bactericidal stretches \n", bacindex);
    
    manvalue = amvalue/(S->wseq_len[cont]-6);
    fprintf(rs, "# This protein has a mean antimicrobial value of %.3f\n\n", manvalue);  
    
    free(window);
    free(residue);

    fclose(kd); 
    fclose(rs);
    
}

float pnormal(double x){
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;

    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);

    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);

    return 0.5*(1.0 + sign*y);
}

float return_value(char residue, float threshold){
    int i=0;
    
    for(i=0; i<strlen(amino_acid_order); i++){
        if(amino_acid_order[i]==residue){
            return (float) amino_acid_values[i];
        }
    }
    
    if(residue == 'X'){
        return threshold;
    }
    
    return (float) -1;
}

void gnuplot_sld(Parameters *par){
    
    char *cmdline, *filename;
    FILE *gp;
    int terminal=0;
    
    cmdline = (char *) calloc(MAXCOMMAND+1, sizeof(char));
    filename = (char *) calloc(FILENAME_MAX, sizeof(char));
    sprintf(filename, "tmp_plot.sh");
    
    if(strcmp(par->graph_file, " ")==0) { 
        sprintf(par->graph_file, "plot.%d.png", getpid());
        terminal=1;
    }

    if((gp = fopen(filename, "w")) == NULL){
        fprintf(stderr, "Error - No gnuplot\n");
    }
    fprintf(gp, "#!/usr/bin/gnuplot\n\n");
    fprintf(gp, "set term png\n");
    fprintf(gp,"set title \"Antimicrobial Profile\"\n");
    fprintf(gp,"set xlabel \"Position\"\n");
    fprintf(gp,"set ylabel \"Score\"\n");
    fprintf(gp,"set grid\n");
    fprintf(gp,"set data style lines\n");
    fprintf(gp,"set output '%s'\n", par->graph_file);
    fprintf(gp, "plot \"%s\" using 1:2 title \"Antimicrobial profile \" with lines\n", par->data_file);
    
    fclose(gp);
    
    sprintf(cmdline, "%s %s", par->gnuplot, filename);
    system(cmdline);
    remove(filename);
    if(terminal){
        sprintf(filename, "/usr/bin/eog");
        if((gp = fopen(filename, "r")) == NULL){
            fprintf(stderr, "Error - No eog\n");  
            sprintf(filename, "/usr/bin/gthumb");
            if((gp = fopen(filename, "r")) == NULL){
                fprintf(stderr, "Error - No gthumb\n"); 
            }
            else {
                sprintf(cmdline, "%s %s &", filename, par->graph_file);
                system(cmdline);
            }
        }
        else {
            sprintf(cmdline, "%s %s &", filename, par->graph_file);
            system(cmdline);
        }
        
    }
    
    free(filename);
    free(cmdline);   
}
