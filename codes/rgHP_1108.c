/*
 * This file calculates the population average values repeatly.
 * Average data are collected from last 10% duration of the simulations. 
 * Please check key parameters for importatnt explanations for the code. 
 ************************************************************************************************
 # Execution #
rgHP_0711.c
./a.out
 # Plot with gnuplot #
gnuplot
plot 'summary.txt' using 3:1 title 'rising-tide' with lines lc rgb 'orange',\
'summary.txt' using 3:2 title 'bet-hedging' with lines lc rgb 'skyblue'
 ************************************************************************************************
 * Key parameters:
 * s1-s5: switches for different submodels, s4==0 for models in the main text, s4==1 for variable cost model
 * num_gen, N, N_par: duration of simulation and population sizes
 * trait 1 variables are termed as "size", whereas trait 2 variables are termed as "patt"
 * trait dissimilarity (d)= (ini_par_size- ini_size) and (ini_par_patt- ini_patt), min_size_diff and min_patt_diff in the output files
 */
#include <stdio.h>
#include <stdlib.h>
#include "../dSFMT-src-2.2.3/dSFMT.c"
#include <time.h>
#include <math.h>

// Functions
    double Mean_array(double p[], int length);
    double rand_normal_dist (double mean, double sd, double u1, double u2);
    double CDF_norm(double x, double mu, double sd);                            // SD has to be larger than 0
    double typ1_error(int n, double x[], double mu[], double sd[]);
    double typ2_error(int n, double x[], double mu_alt[], double sd_alt[]);
    double nest_survive(int n);

// Main function
int main (void){
    // Switches
        int s1=1;           // Using single or two traits for recognition
        int s2=1;           // Which trait is used (only for single trait, s1==1)
        int s3=1;           // Enable coevolution of acceptance thresholds and undesirable senders (creating an actual undesirable sender population)
        int s4=0;           // Specify the way of reproduction (0: brood parasitism, 1: kin recognition)
        int s5=1;           // Another threshold located at the 1% CDF of the distribution of desirable recipients (accepted only when traits are above this threshld)

    // Basic variables
        int i,j,k,l;                        // For loop counters
        double u1, u2;                      // Storing random numbers
        int num_gen=        3000;           // Number of generation simulated
        int N=              1000;           // Number of adults and nests in each generation
        int N_par=          200;            // Number of undesirable senders in each generation
        int num_geno=       4;              // Number of genotype information (2 traits+ 2 thresholds)
        int num_trait=      6;              // Number of genotype and phenotype information (genotype+ 2 phenotypes)
        int num_offspring=  3;              // Number of offspring produced in each nest (desirable sender)
        double mut_rate=    0.001;          // Mutation rate
        int off_log, par_log, nest_log, idx, off_count, abn, ini_check, reini_count, end_gen, off_pass, par_pass;
        double mean_size, mean_patt, mean_par_size, mean_par_patt, mean_p_size, mean_p_patt, mean_size_thr, mean_patt_thr, upr_size, lwr_size, upr_patt, lwr_patt, par_effort, tmp1, tmp2, typ1_err, typ2_err, nest_survival, rej_rate, acp_rate;

    // Properties of the population
        double sd_ratio=    1.2;            // Square-rooted ratio of difference in phenotypic variations between desirable senders and undesirable senders
        double sd_size=     0.5/sd_ratio;   // Phenotypic variation of desirable sender (trait 1)
        double sd_par_size= 0.5*sd_ratio;   // Phenotypic variation of undesirable sender (trait 1)
        double sd_patt=     0.5/sd_ratio;   // Phenotypic variation of desirable sender (trait 2)
        double sd_par_patt= 0.5*sd_ratio;   // Phenotypic variation of undesirable sender (trait 2)
        double mut_size=    1.0;            // Scaling factor of the mutations in trait 1, 1.0
        double mut_patt=    1.0;            // Scaling factor of the mutations in trait 2
        double mut_par_size=0.25;           // Scaling factor of the mutations in trait 1 (undesirable sender)
        double mut_par_patt=0.25;           // Scaling factor of the mutations in trait 2 (undesirable sender)
        double mut_size_thr=1.0;            // Scaling factor of the mutations in thresholds of trait 1, 0.3
        double mut_patt_thr=1.0;            // Scaling factor of the mutations in thresholds of trait 2

    // Initial conditions
        // Trait values
            double ini_size=    10.0;                       // Initial value of trait 1 (average of the population)
            double ini_sd_size= 0.0;                        // Genetic variation (standard deviation) of trait 1 in the population (assuming no variation)
            double ini_patt=    10.0;                       // Initial average value of trait 2
            double ini_sd_patt= 0.0;                        // Genetic variation of trait 2 
            double ini_par_size=15.0;                       // Initial average value of trait 1 of the undesirable senders (20)
            double ini_sd_psize= sqrt(ini_par_size);        // Genetic variation of trait 1 (undesirable sender)
            double ini_par_patt=15.0;                       // Initial average value of trait 2 of the undesirable senders
            double ini_sd_ppatt= sqrt(ini_par_patt);        // Genetic variation of trait 2 (undesirable sender)
        // Thresholds
            double ini_size_thr=2*sd_size;                  // Initial value of threshold of trait 1 (threshold 1), represented as deviation from average trait value
            double ini_sd_sthr= sqrt(ini_size_thr);         // Genetic variation of threshold 1
            double ini_patt_thr=2*sd_patt;                  // Initial value of threshold 2
            double ini_sd_pthr= sqrt(ini_patt_thr);         // Genetic variation of threshold 2
            if(ini_size_thr< 1) ini_sd_sthr= pow(ini_size_thr,2);
            if(ini_patt_thr< 1) ini_sd_pthr= pow(ini_patt_thr,2);
            ini_size_thr+= ini_size;                        // Converting the deviation to actual trait values
            ini_patt_thr+= ini_patt;                        // Converting the deviation to actual trait values
            double buttom_thr= 9.030688;                    // Calculated from R (qnorm(0.01, 10, 0.5/1.2)), for s5==1 only

    // Creating temporal spaces
        double Pop[N][num_geno];                            // The population matrix with N individuals and num_geno loci
        double Pop_par[N_par][num_geno];                    // The population matrix with N_par individuals and num_geno loci
        double Offspring[num_offspring*N][num_trait];       // The matrix storing the surviving offspring in the current breeding event
        double Offs_log[num_offspring*N][num_trait];        // The matrix storing all the produced offspring for undesirable sender to change its trait values
        double Offs_par[N][num_trait];                      // The matrix storing the surviving undesirable sender offspring in the current breeding event
        double temporary[num_offspring*N];                  // A temporary array storing offspring phenotype values
        double temp_off[num_offspring][num_trait];          // The matrix storing all offspring of the focal nest
        double temp_par[num_trait];                         // The array storing the undesirable sender of the focal nest
        double list_thr[s1];                                // The array storing the thresholds for calculating errors
        double list_mean[s1];                               // The array storing average desirable sender trait values for errors
        double list_sd[s1];                                 // The array storing standard deviations of desirable sender traits for errors
        double list_par_mean[s1];                           // The array storing average undesirable sender trait values for errors
        double list_par_sd[s1];                             // The array storing standaed deviations of undesirable sender traits for errors
        // Index
            int i_size=         0;  // genotype
            int i_size_thr=     1;  // behavioral criteria
            int i_patt=         2;  // genotype
            int i_patt_thr=     3;  // behavioral criteria
            int i_p_size=       4;  // phenotype
            int i_p_patt=       5;  // phenotype
        // undesirable sender params
            double par_size, par_patt;
    // Repeating 
        int ii, jj, kk, log_index;
        int rep=                100;
        double log_ratio=       0.1;
        int length_log=         round(num_gen*log_ratio);
        int ini_log=            num_gen- length_log;
        int min_length=         num_gen/3;
        double avg_size, avg_patt, avg_size_thr, avg_patt_thr, avg_par_size, avg_par_patt;
        // Log matrices
            double typ1_log[length_log];
            double typ2_log[length_log];
            double size_log[length_log];
            double patt_log[length_log];
            double size_thr_log[length_log];
            double patt_thr_log[length_log];
            double par_size_log[length_log];
            double par_patt_log[length_log];
            double rej_log[length_log];
            double acp_log[length_log];
   	// Initialization of random number genertor
		int seed;
		dsfmt_t dsfmt;
		seed= time(NULL);
		if(seed==0)seed= 1;
		dsfmt_init_gen_rand(&dsfmt,seed);

	// Output
		FILE *out;
		out= fopen("summary.txt","w");
        if(s1==1 && s2==1 && s3==0 && s4==0)    fprintf(out,"sd_size\tsd_par_size\tmin_size_dif\ttyp1_err\ttyp2_err\t\
            avg_size\tavg_size_thr\tavg_par_size\tend_gen\n");
        if(s1==1 && s2==2 && s3==0 && s4==0)    fprintf(out,"sd_patt\tsd_par_patt\tmin_patt_dif\ttyp1_err\ttyp2_err\t\
            avg_patt\tavg_patt_thr\tavg_par_patt\tend_gen\n");
        if(s1==2 && s3==0 && s4==0)             fprintf(out,"sd_size\tsd_par_size\tmin_size_dif\tsd_patt\tsd_par_patt\tmin_patt_dif\ttyp1_err\ttyp2_err\t\
            avg_size\tavg_size_thr\tavg_patt\tavg_patt_thr\tavg_par_size\tavg_par_patt\tend_gen\n");
        if(s1==1 && s2==1 && s3==1 && s4==0)    fprintf(out,"mut_par_size\ttyp1_err\ttyp2_err\tavg_size\tavg_size_thr\tavg_par_size\tend_gen\n");
        if(s1==1 && s2==2 && s3==1 && s4==0)    fprintf(out,"mut_par_patt\ttyp1_err\ttyp2_err\tavg_patt\tavg_patt_thr\tavg_par_patt\tend_gen\n");
        if(s1==2 && s3==1 && s4==0)             fprintf(out,"mut_par_size\tmut_par_patt\ttyp1_err\ttyp2_err\t\
            avg_size\tavg_size_thr\tavg_par_size\tavg_patt\tavg_patt_thr\tavg_par_patt\tend_gen\n");
        if(s1==1 && s2==1 && s4==1)             fprintf(out,"mut_par_size\ttyp1_err\ttyp2_err\trej_rate\tacp_rate\tavg_size\tavg_size_thr\tavg_par_size\tend_gen\n");
        if(s1==1 && s2==2 && s4==1)             fprintf(out,"mut_par_patt\ttyp1_err\ttyp2_err\trej_rate\tacp_rate\tavg_patt\tavg_patt_thr\tavg_par_patt\tend_gen\n");
        if(s1==2 && s4==1)                      fprintf(out,"mut_par_size\tmut_par_patt\ttyp1_err\ttyp2_err\t\
            rej_rate\tacp_rate\tavg_size\tavg_size_thr\tavg_par_size\tavg_patt\tavg_patt_thr\tavg_par_patt\tend_gen\n");

    // Main loop
    for (jj=0; jj<5; jj++){
        // 1x reprtitions for s3==0
        if(s1==1 && s2==1 && s3==0) ini_par_size= ini_size+ 0.5*jj;
        if(s1==1 && s2==2 && s3==0) ini_par_patt= ini_patt+ 0.5*jj;
        if(s1==2 && s3==0){
            ini_par_patt= 0.5*jj;
            ini_par_patt= 0.5*(4-jj);
        }
        // 5x repetitions for s3==1
        if(s1==1 && s2==1 && s3==1) mut_par_size= 0.25;
        if(s1==1 && s2==2 && s3==1) mut_par_patt= 0.25;
        if(s1==2 && s3==1){
            mut_par_size= 0.25;
            mut_par_patt= 0.25;
        }
    for (ii=0; ii< rep; ii++){
        abn= log_index= 0;
        for (i=0; i< length_log; i++){
            typ1_log[i]= NAN;
            typ2_log[i]= NAN;
            size_log[i]= NAN;
            patt_log[i]= NAN;
            size_thr_log[i]= NAN;
            patt_thr_log[i]= NAN;
            if(s3==1){
                par_size_log[i]= NAN;
                par_patt_log[i]= NAN;
            }
        }
        // Initializing
            for (i= 0; i< N; i++){
                Pop[i][i_size]=     rand_normal_dist(ini_size,ini_sd_size,dsfmt_genrand_open_open(&dsfmt),dsfmt_genrand_open_open(&dsfmt));
                Pop[i][i_size_thr]= rand_normal_dist(ini_size_thr,ini_sd_sthr,dsfmt_genrand_open_open(&dsfmt),dsfmt_genrand_open_open(&dsfmt));
                Pop[i][i_patt]=     rand_normal_dist(ini_patt,ini_sd_patt,dsfmt_genrand_open_open(&dsfmt),dsfmt_genrand_open_open(&dsfmt));
                Pop[i][i_patt_thr]= rand_normal_dist(ini_patt_thr,ini_sd_pthr,dsfmt_genrand_open_open(&dsfmt),dsfmt_genrand_open_open(&dsfmt));
            }
            if(s3==0){
                par_size= ini_par_size;
                par_patt= ini_par_patt;
            }
            if(s3==1){
                for(i=0; i< N_par; i++){
                    Pop_par[i][i_size]= rand_normal_dist(ini_par_size,ini_sd_psize,dsfmt_genrand_open_open(&dsfmt),dsfmt_genrand_open_open(&dsfmt));
                    Pop_par[i][i_patt]= rand_normal_dist(ini_par_patt,ini_sd_ppatt,dsfmt_genrand_open_open(&dsfmt),dsfmt_genrand_open_open(&dsfmt));
                }
            }
        // Start simulation
        for (i= 0; i< num_gen; i++){
            off_count= off_log= par_log= off_pass= par_pass= 0;
            // Erasing the offspring matrix
                for(j=0; j< num_offspring*N; j++){
                    for(k=0; k< num_trait; k++) {
                        Offspring[j][k]=    NAN;
                        Offs_log[j][k]=     NAN;
                    }
                }
                if(s3==1){
                    for(j=0; j< N; j++){
                        for(k=0; k< num_trait; k++) Offs_par[j][k]= NAN;
                    }
                }
            // For each nest
                for(j=0; j< N; j++){
                // Produce offspring
                    for(k=0; k< num_offspring; k++){
                            u1= dsfmt_genrand_open_open(&dsfmt);
                        if(u1< mut_rate){
                            temp_off[k][i_size_thr]= rand_normal_dist(Pop[j][i_size_thr],mut_size_thr,dsfmt_genrand_open_open(&dsfmt),dsfmt_genrand_open_open(&dsfmt));
                        }
                        else temp_off[k][i_size_thr]= Pop[j][i_size_thr];
                            u1= dsfmt_genrand_open_open(&dsfmt);
                        if(u1< mut_rate){
                            temp_off[k][i_patt_thr]= rand_normal_dist(Pop[j][i_patt_thr],mut_patt_thr,dsfmt_genrand_open_open(&dsfmt),dsfmt_genrand_open_open(&dsfmt));
                        }
                        else temp_off[k][i_patt_thr]= Pop[j][i_patt_thr];
                        // Determine the phenotype
                        temp_off[k][i_p_size]= rand_normal_dist(Pop[j][i_size],sd_size,dsfmt_genrand_open_open(&dsfmt),dsfmt_genrand_open_open(&dsfmt));
                        temp_off[k][i_p_patt]= rand_normal_dist(Pop[j][i_patt],sd_patt,dsfmt_genrand_open_open(&dsfmt),dsfmt_genrand_open_open(&dsfmt));
                        for(l=0; l< num_trait; l++) Offs_log[off_count][l]= temp_off[k][l];
                        off_count+= 1;
                    }
                    // Placing the undesirable sender egg
                        if(s3==0){
                            temp_par[i_p_size]= rand_normal_dist(par_size,sd_par_size,dsfmt_genrand_open_open(&dsfmt),dsfmt_genrand_open_open(&dsfmt));
                            temp_par[i_p_patt]= rand_normal_dist(par_patt,sd_par_patt,dsfmt_genrand_open_open(&dsfmt),dsfmt_genrand_open_open(&dsfmt));
                        }
                        if(s3==1){
                            idx= round(dsfmt_genrand_open_open(&dsfmt)*(N_par-1));
                            temp_par[i_size_thr]= NAN;
                            temp_par[i_patt_thr]= NAN;
                            // Determine the phenotype
                            temp_par[i_p_size]= rand_normal_dist(Pop_par[idx][i_size],sd_par_size,dsfmt_genrand_open_open(&dsfmt),dsfmt_genrand_open_open(&dsfmt));
                            temp_par[i_p_patt]= rand_normal_dist(Pop_par[idx][i_patt],sd_par_patt,dsfmt_genrand_open_open(&dsfmt),dsfmt_genrand_open_open(&dsfmt));
                        }       
                    // Rejection
                        // Accepting undesirable sender kills the entire nest
                        if(s4==0){
                            if(s1== 1 && s2==1){
                                upr_size= Pop[j][i_size_thr];
                                // If the undesirable sender fails
                                if( (temp_par[i_p_size]> upr_size) || (s5==1 && temp_par[i_p_size]<= buttom_thr) ){
                                    for(k=0; k< num_offspring; k++){
                                        // If an offspring survives
                                        if(upr_size>= temp_off[k][i_p_size]){
                                            if(s5==1){
                                                if(temp_off[k][i_p_size]> buttom_thr){
                                                    for(l=0; l< num_trait; l++) Offspring[off_log][l]= temp_off[k][l];
                                                    off_log+= 1;
                                            }}
                                            else{
                                                for(l=0; l< num_trait; l++) Offspring[off_log][l]= temp_off[k][l];
                                                off_log+= 1;
                                        }}
                                }}
                                // The undesirable sender survives
                                else{
                                    if(s3==1) for(k=0; k< num_trait; k++) Offs_par[par_log][k]= temp_par[k];
                                    par_log+= 1;
                                }
                            }
                            if(s1== 1 && s2==2){
                                upr_patt= Pop[j][i_patt_thr];
                                // If the undesirable sender fails
                                if( (temp_par[i_p_patt]> upr_patt) || (s5==1 && temp_par[i_p_patt]<= buttom_thr) ){
                                    for(k=0; k< num_offspring; k++){
                                        // If an offspring survives
                                        if(upr_patt>= temp_off[k][i_p_patt]){
                                            if(s5==1){
                                                if(temp_off[k][i_p_patt]> buttom_thr){
                                                    for(l=0; l< num_trait; l++) Offspring[off_log][l]= temp_off[k][l];
                                                    off_log+= 1;
                                            }}
                                            else{
                                                for(l=0; l< num_trait; l++) Offspring[off_log][l]= temp_off[k][l];
                                                off_log+= 1;
                                        }}
                                }}
                                // The undesirable sender survives
                                else{
                                    if(s3==1) for(k=0; k< num_trait; k++) Offs_par[par_log][k]= temp_par[k];
                                    par_log+= 1;
                                }
                            }
                            if(s1==2){
                                upr_size= Pop[j][i_size_thr];
                                upr_patt= Pop[j][i_patt_thr];
                                // If the undesirable sender fails
                                if( (temp_par[i_p_size]> upr_size) || (temp_par[i_p_patt]> upr_patt) || (s5==1 && temp_par[i_p_size]<= buttom_thr) || (s5==1 && temp_par[i_p_patt]<= buttom_thr) ){
                                    for(k=0; k< num_offspring; k++){
                                        if(upr_size>= temp_off[k][i_p_size] && upr_patt>= temp_off[k][i_p_patt]){
                                            if(s5==1){
                                                if(temp_off[k][i_p_size]> buttom_thr && temp_off[k][i_p_patt]> buttom_thr){
                                                    for(l=0; l< num_trait; l++) Offspring[off_log][l]= temp_off[k][l];
                                                    off_log+= 1; 
                                            }}
                                            else{
                                                for(l=0; l< num_trait; l++) Offspring[off_log][l]= temp_off[k][l];
                                                off_log+= 1; 
                                        }}
                                }}
                                // The undesirable sender survives
                                else{
                                    if(s3==1) for(k=0; k< num_trait; k++) Offs_par[par_log][k]= temp_par[k];
                                    par_log+= 1;
                                }
                            }
                        }
                        if(s4==1){
                            nest_log= 0;
                            if(s1==1 && s2==1){
                                upr_size= Pop[j][i_size_thr];
                                // Calculate the offspring that is not rejected
                                if(temp_par[i_p_size]<= upr_size){
                                    if(s5==1){
                                        if(temp_par[i_p_size]> buttom_thr) nest_log+= 1;
                                    }
                                    else nest_log+= 1;
                                }
                                for(k=0; k< num_offspring; k++){
                                    if(upr_size>= temp_off[k][i_p_size]){
                                        if(s5==1){
                                            if(temp_off[k][i_p_size]> buttom_thr) nest_log+= 1;
                                        }
                                        else nest_log+= 1;
                                    }
                                }
                                // Calculate the survival rate
                                nest_survival= nest_survive(nest_log);
                                // Acceptance/ rejection and survival processes
                                // The undesirable sender
                                u1= dsfmt_genrand_open_open(&dsfmt);
                                if(temp_par[i_p_size]<= upr_size){
                                    if(s5==1){
                                        if(temp_par[i_p_size]> buttom_thr){
                                            par_pass+= 1;
                                            if(u1< nest_survival) {
                                                if(s3==1) for(k=0; k< num_trait; k++) Offs_par[par_log][k]= temp_par[k];
                                                par_log+= 1;
                                    }}}
                                    else{
                                        par_pass+= 1;
                                        if(u1< nest_survival) {
                                            if(s3==1) for(k=0; k< num_trait; k++) Offs_par[par_log][k]= temp_par[k];
                                            par_log+= 1;
                                    }}
                                }
                                // The offspring
                                for(k=0; k< num_offspring; k++){
                                    u1= dsfmt_genrand_open_open(&dsfmt);
                                    if(upr_size>= temp_off[k][i_p_size]){
                                        if(s5==1){
                                            if(temp_off[k][i_p_size]> buttom_thr){
                                                off_pass+= 1;
                                                if(u1< nest_survival) {
                                                    for(l=0; l< num_trait; l++) Offspring[off_log][l]= temp_off[k][l];
                                                    off_log+= 1;
                                        }}}
                                        else{
                                            off_pass+= 1;
                                            if(u1< nest_survival) {
                                                for(l=0; l< num_trait; l++) Offspring[off_log][l]= temp_off[k][l];
                                                off_log+= 1;
                                        }}
                                }}
                            }
                            if(s1==1 && s2==2){
                                upr_patt= Pop[j][i_patt_thr];
                                // Calculate the offspring that is not rejected
                                if(temp_par[i_p_patt]<= upr_patt){
                                    if(s5==1){
                                        if(temp_par[i_p_patt]> buttom_thr) nest_log+= 1;
                                    }
                                    else nest_log+= 1;
                                }
                                for(k=0; k< num_offspring; k++){
                                    if(upr_patt>= temp_off[k][i_p_patt]){
                                        if(s5==1){
                                            if(temp_off[k][i_p_patt]> buttom_thr) nest_log+= 1;
                                        }
                                        else nest_log+= 1;
                                    }
                                }
                                nest_survival= nest_survive(nest_log);
                                // Acceptance/ rejection and survival processes
                                u1= dsfmt_genrand_open_open(&dsfmt);
                                if(temp_par[i_p_patt]<= upr_patt){
                                    if(s5==1){
                                        if(temp_par[i_p_patt]> buttom_thr){
                                            par_pass+= 1;
                                            if(u1< nest_survival){
                                                if(s3==1) for(k=0; k< num_trait; k++) Offs_par[par_log][k]= temp_par[k];
                                                par_log+= 1;
                                    }}}
                                    else{
                                        par_pass+= 1;
                                        if(u1< nest_survival){
                                            if(s3==1) for(k=0; k< num_trait; k++) Offs_par[par_log][k]= temp_par[k];
                                            par_log+= 1;
                                    }}
                                }
                                for(k=0; k< num_offspring; k++){
                                    u1= dsfmt_genrand_open_open(&dsfmt);
                                    if(upr_patt>= temp_off[k][i_p_patt]){
                                        if(s5==1){
                                            if(temp_off[k][i_p_patt]> buttom_thr){
                                                off_pass+= 1;
                                                if(u1< nest_survival){
                                                    for(l=0; l< num_trait; l++) Offspring[off_log][l]= temp_off[k][l];
                                                    off_log+= 1;
                                                }
                                            }
                                        }
                                        else{
                                            off_pass+= 1;
                                            if(u1< nest_survival){
                                                for(l=0; l< num_trait; l++) Offspring[off_log][l]= temp_off[k][l];
                                                off_log+= 1;
                                        }}
                                }}
                            }
                            if(s1==2){
                                upr_size= Pop[j][i_size_thr];
                                upr_patt= Pop[j][i_patt_thr];
                                // Calculate the offspring that is not rejected
                                if(temp_par[i_p_size]<= upr_size && temp_par[i_p_patt]<= upr_patt){
                                    if(s5==1){
                                        if(temp_par[i_p_size]> buttom_thr && temp_par[i_p_patt]> buttom_thr) nest_log+= 1;
                                    }
                                    else nest_log+= 1;
                                }
                                for(k=0; k< num_offspring; k++){
                                    if(upr_size>= temp_off[k][i_p_size] && upr_patt>= temp_off[k][i_p_patt]){
                                        if(s5==1){
                                            if(temp_off[k][i_p_size]> buttom_thr && temp_off[k][i_p_patt]> buttom_thr) nest_log+= 1;
                                        }
                                        else nest_log+= 1;
                                    }
                                }
                                nest_survival= nest_survive(nest_log);
                                // Acceptance/ rejection and survival processes
                                u1= dsfmt_genrand_open_open(&dsfmt);
                                if(temp_par[i_p_size]<= upr_size && temp_par[i_p_patt]<= upr_patt){
                                    if(s5==1){
                                        if(temp_par[i_p_size]> buttom_thr && temp_par[i_p_patt]> buttom_thr){
                                            par_pass+= 1;
                                            if(u1< nest_survival){
                                                if(s3==1) for(k=0; k< num_trait; k++) Offs_par[par_log][k]= temp_par[k];
                                                par_log+= 1;
                                    }}}
                                    else{
                                        par_pass+= 1;
                                        if(u1< nest_survival){
                                            if(s3==1) for(k=0; k< num_trait; k++) Offs_par[par_log][k]= temp_par[k];
                                            par_log+= 1;
                                    }}
                                }
                                for(k=0; k< num_offspring; k++){
                                    u1= dsfmt_genrand_open_open(&dsfmt);
                                    if(upr_size>= temp_off[k][i_p_size] && upr_patt>= temp_off[k][i_p_patt]){
                                        if(s5==1){
                                            if(temp_off[k][i_p_size]> buttom_thr && temp_off[k][i_p_patt]> buttom_thr){
                                                off_pass+= 1;
                                                if(u1< nest_survival){
                                                    for(l=0; l< num_trait; l++) Offspring[off_log][l]= temp_off[k][l];
                                                    off_log+= 1;
                                        }}}
                                        else{
                                            off_pass+= 1;
                                            if(u1< nest_survival){
                                                for(l=0; l< num_trait; l++) Offspring[off_log][l]= temp_off[k][l];
                                                off_log+= 1;
                                        }}
                            }}}
                        }
                }
            // Creating the new population
                if(off_log>0){
                    for(j=0; j< N; j++){
                        idx= round(dsfmt_genrand_open_open(&dsfmt)*(off_log-1)); // array starts from 0, while off_log starts from 1
                        for(k=0; k< num_geno; k++) Pop[j][k]= Offspring[idx][k];
                    }
                }
                else abn= 1;
            // Changing the trait value of undesirable sender
                if(s3==0){
                    if(s1== 1&& s2==1 && off_log> 0){
                        for(j=0; j< off_count; j++) temporary[j]= Offs_log[j][i_p_size];
                        mean_p_size= Mean_array(temporary, off_count);
                    }
                    if(s1== 1&& s2==2 && off_log> 0){
                        for(j=0; j< off_count; j++) temporary[j]= Offs_log[j][i_p_patt];
                        mean_p_patt= Mean_array(temporary, off_count);
                    }
                    if(s1== 2 && off_log> 0){
                        for(j=0; j< off_count; j++) temporary[j]= Offs_log[j][i_p_size];
                        mean_p_size= Mean_array(temporary, off_count);
                        for(j=0; j< off_count; j++) temporary[j]= Offs_log[j][i_p_patt];
                        mean_p_patt= Mean_array(temporary, off_count);
                    }
                }
                if(s3==1){
                    if(par_log> 0){
                        for(j=0; j< N_par; j++){
                            idx= round(dsfmt_genrand_open_open(&dsfmt)*(par_log-1));
                            for(k=0; k< num_geno; k++) Pop_par[j][k]= Offs_par[idx][k];
                    }}
                    else abn= 1;
                }

            // Calculating average values
            if(off_log>0){
                // Average trait values and threshoulds
                    for(j=0; j< N; j++) temporary[j]= Pop[j][i_size];
                    mean_size= Mean_array(temporary, N);
                    for(j=0; j< N; j++) temporary[j]= Pop[j][i_patt];
                    mean_patt= Mean_array(temporary, N);
                    for(j=0; j< N; j++) temporary[j]= Pop[j][i_size_thr];
                    mean_size_thr= Mean_array(temporary, N);
                    for(j=0; j< N; j++) temporary[j]= Pop[j][i_patt_thr];
                    mean_patt_thr= Mean_array(temporary, N);
                    if(s3==1){
                        for(j=0; j<N_par; j++) temporary[j]= Pop_par[j][i_size];
                        mean_par_size= Mean_array(temporary, N_par);
                        for(j=0; j<N_par; j++) temporary[j]= Pop_par[j][i_patt];
                        mean_par_patt= Mean_array(temporary, N_par);
                    }
                // Errors
                if(s4==0){
                    typ2_err= (double)par_log/N;
                    typ1_err= (double)(N-off_log/num_offspring-par_log)/N;
                }
                if(s4==1){
                    typ2_err= (double)par_pass/N;
                    typ1_err= (double)(N-off_pass/num_offspring)/N;
                }
                log_index= i%length_log;
                typ1_log[log_index]= typ1_err;
                typ2_log[log_index]= typ2_err;
                size_log[log_index]= mean_size;
                patt_log[log_index]= mean_patt;
                size_thr_log[log_index]= mean_size_thr;
                patt_thr_log[log_index]= mean_patt_thr;
                if(s3==1){
                    par_size_log[log_index]= mean_par_size;
                    par_patt_log[log_index]= mean_par_patt;
                }
                else{
                    par_size_log[log_index]= par_size;
                    par_patt_log[log_index]= par_patt;
                }
                // Saving additional acceptance/ rejection data (after survival process)
                if(s4==1){
                    rej_log[log_index]= (double)(N-off_log/num_offspring)/N;
                    acp_log[log_index]= (double)par_log/N;
                }
            }
            // Record current generation number
            end_gen= i;
            // Terminating the simulation due to no desirable senders
            if(abn==1) i= num_gen;
        }
        // Printing
        if(end_gen>= min_length){
            typ1_err=       Mean_array(typ1_log, length_log);
            typ2_err=       Mean_array(typ2_log, length_log);
            avg_size=       Mean_array(size_log, length_log);
            avg_patt=       Mean_array(patt_log, length_log);
            avg_size_thr=   Mean_array(size_thr_log, length_log);
            avg_patt_thr=   Mean_array(patt_thr_log, length_log);
            avg_par_size=   Mean_array(par_size_log, length_log);
            avg_par_patt=   Mean_array(par_patt_log, length_log);
            if(s4==1){
                rej_rate=   Mean_array(rej_log, length_log);
                acp_rate=   Mean_array(acp_log, length_log);
            }
        // If NAN is printed, the simulation lasted less than length_log generations
        if(s1==1 && s2==1 && s3==0 && s4==0)    fprintf(out,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\n",
            sd_size, sd_par_size, (ini_par_size- ini_size), typ1_err, typ2_err, avg_size, avg_size_thr, avg_par_size, end_gen);
        if(s1==1 && s2==2 && s3==0 && s4==0)    fprintf(out,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\n",
            sd_patt, sd_par_patt, (ini_par_patt- ini_patt), typ1_err, typ2_err, avg_patt, avg_patt_thr, avg_par_patt, end_gen);
        if(s1==2 && s3==0 && s4==0)             fprintf(out,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\n",
            sd_size, sd_par_size, (ini_par_size- ini_size), sd_patt, sd_par_patt, (ini_par_patt- ini_patt), typ1_err, typ2_err,
            avg_size, avg_size_thr, avg_patt, avg_patt_thr, avg_par_size, avg_par_patt, end_gen);
        if(s1==1 && s2==1 && s3==1 && s4==0)    fprintf(out,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\n",
            mut_par_size, typ1_err, typ2_err, avg_size, avg_size_thr, avg_par_size, end_gen);
        if(s1==1 && s2==2 && s3==1 && s4==0)    fprintf(out,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\n",
            mut_par_patt, typ1_err, typ2_err, avg_patt, avg_patt_thr, avg_par_patt, end_gen);
        if(s1==2 && s3==1 && s4==0)             fprintf(out,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\n",
            mut_par_size, mut_par_patt, typ1_err, typ2_err, avg_size, avg_size_thr, avg_par_size, avg_patt, avg_patt_thr, avg_par_patt, end_gen);
        if(s1==1 && s2==1 && s4==1)             fprintf(out,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\n",
            mut_par_size, typ1_err, typ2_err, rej_rate, acp_rate, avg_size, avg_size_thr, avg_par_size, end_gen);
        if(s1==1 && s2==2 && s4==1)             fprintf(out,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\n",
            mut_par_patt, typ1_err, typ2_err, rej_rate, acp_rate, avg_patt, avg_patt_thr, avg_par_patt, end_gen);
        if(s1==2 && s4==1)                      fprintf(out,"%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\n",
            mut_par_size, mut_par_patt, typ1_err, typ2_err, rej_rate, acp_rate,
            avg_size, avg_size_thr, avg_par_size, avg_patt, avg_patt_thr, avg_par_patt, end_gen);
        }
        else ii-=1;
    }}
    // Free memory
    fclose(out);
    return 0;
}
////////////////////// Functions are defined in below ////////////////////////
double Mean_array(double p[], int length){
    int i;
    double temp, sum;
    temp= sum= 0.0;
    if(length>0) {
        for (i=0; i<length; ++i) sum+= p[i];
        temp= sum/(double)length;
        return temp;
    }
    else return NAN;
}
double rand_normal_dist (double mean, double sd, double u1, double u2)
{
    // Using Box-Muller method to generate pseudo-normal distributed numbers in [0,1]
    // Weird bug happens if random number is generated within the function
    double z1;
	z1= sqrt(-2* log(u1))* cos(2* M_PI* u2);
	return z1*sd+ mean;
}
double CDF_norm(double x, double mu, double sd){
    double z1;
    z1= (x-mu)/sqrt(2)/sd;
    return (1+erf(z1))/2;
}
double typ1_error(int n, double x[], double mu[], double sd[]){
    int i;
    double tmp= 1.0;
    double p[n];
    for(i=0; i<n; i++) p[i]= CDF_norm(x[i], mu[i], sd[i]);
    for(i=0; i<n; i++) tmp= tmp*p[i];
    return 1-tmp;
}
double typ2_error(int n, double x[], double mu_alt[], double sd_alt[]){
    int i;
    double tmp= 1.0;
    double p[n];
    for(i=0; i<n; i++) p[i]= CDF_norm(x[i], mu_alt[i], sd_alt[i]);
    for(i=0; i<n; i++) tmp= tmp*p[i];
    return tmp;
}
double nest_survive(int n){
    double tmp= 0.0;
    if (n==1) tmp= 0.95;
    if (n==2) tmp= 0.90;
    if (n==3) tmp= 0.85;
    if (n==4) tmp= 0.25;
    return tmp;
}