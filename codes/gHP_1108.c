/*
 * This file generates the time series of trait values from single simulation. 
 * Please check key parameters for importatnt explanations for the code. 
 ************************************************************************************************
 # Please paste these lines to the terminal to execute the code #
gHP_1108.c
./a.out
 # Plot with gnuplot #
gnuplot
set terminal qt 0
set yrange [5:15]
set xrange [0:50]
plot 'summary.txt' using 1:2 title 'desirable recipient' with lines lc rgb 'orange',\
'summary.txt' using 1:6 title 'undesirable recipient' with lines lc rgb 'skyblue',\
'summary.txt' using 1:3 title 'threshold' with lines lc rgb 'black'

set terminal qt 1
set yrange [5:15]
set xrange [0:50]
plot 'summary.txt' using 1:4 title 'desirable recipient' with lines lc rgb 'orange',\
'summary.txt' using 1:8 title 'undesirable recipient' with lines lc rgb 'skyblue',\
'summary.txt' using 1:5 title 'threshold' with lines lc rgb 'black'
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
    double SD_array(double p[], int length);
    double rand_normal_dist (double mean, double sd, double u1, double u2);
    double CDF_norm(double x, double mu, double sd);                            // SD has to be larger than 0
    double typ1_error(int n, double x[], double mu[], double sd[]);
    double typ2_error(int n, double x[], double mu_alt[], double sd_alt[]);
    double nest_survive(int n);

// Main function
int main (void){
    // Switches
        int s1=2;           // Using single or two traits for recognition
        int s2=1;           // Which trait is used (only for single trait, s1==1)
        int s3=1;           // Enable coevolution of acceptance thresholds and undesirable senders (creating an actual undesirable sender population)
        int s4=0;           // Specify the way of reproduction (0: brood parasitism, 1: non-kin parasitism)
        int s5=1;           // Another threshold located at the 1% CDF of the distribution of desirable recipients (accepted only when traits are above this threshld)

    // Basic variables
        int i,j,k,l;                        // For loop counters
        double u1, u2;                      // Storing random numbers
        int num_gen=        500;            // Number of generation simulated
        int N=              1000;           // Number of adults and nests in each generation
        int N_par=          200;            // Number of undesirable senders in each generation
        int num_geno=       4;              // Number of genotype information (2 traits+ 2 thresholds)
        int num_trait=      6;              // Number of genotype and phenotype information (genotype+ 2 phenotypes)
        int num_offspring=  3;              // Number of offspring produced in each nest (desirable sender)
        double mut_rate=    0.001;          // Mutation rate
        int off_log, par_log, nest_log, idx, off_count, ini_check, reini_count, par_pass, off_pass;
        double mean_size, mean_patt, mean_par_size, mean_par_patt, mean_p_size, mean_p_patt, mean_size_thr, mean_patt_thr, upr_size, lwr_size, upr_patt, lwr_patt, par_effort, tmp1, tmp2, typ1_err, typ2_err, nest_survival, print_psize_sd, print_ppatt_sd;

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
        // Trait values (with ini_sd allowing genotypic variations within the population)
        double ini_size=    10.0;                       // Initial value of trait 1 (average of the population)
        double ini_sd_size= 0.0;                        // Genetic variation (standard deviation) of trait 1 in the population (assuming no variation)
        double ini_patt=    10.0;                       // Initial average value of trait 2
        double ini_sd_patt= 0.0;                        // Genetic variation of trait 2 
        double ini_par_size=15.0;                       // Initial average value of trait 1 of the undesirable senders (20)
        double ini_sd_psize= sqrt(ini_par_size);        // Genetic variation of trait 1 (undesirable sender)
        double ini_par_patt=15.0;                       // Initial average value of trait 2 of the undesirable senders
        double ini_sd_ppatt= sqrt(ini_par_size);        // Genetic variation of trait 2 (undesirable sender)
        // Thresholds (with ini_sd allowing variation of thresholds among individuals)
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
        double Pop[N][num_geno];                        // The population matrix with N individuals and num_geno loci
        double Pop_par[N_par][num_geno];                // The population matrix with N_par individuals and num_geno loci
        double Offspring[num_offspring*N][num_trait];   // The matrix storing the surviving offspring in the current breeding event
        double Offs_log[num_offspring*N][num_trait];    // The matrix storing all the produced offspring for undesirable sender to change its trait values
        double Offs_par[N][num_trait];                  // The matrix storing the surviving undesirable sender offspring in the current breeding event
        double temporary[num_offspring*N];              // A temporary array storing offspring phenotype values
        double temp_off[num_offspring][num_trait];      // The matrix storing all offspring of the focal nest
        double temp_par[num_trait];                     // The array storing the undesirable sender of the focal nest
        double list_thr[s1];                            // The array storing the thresholds for calculating errors
        double list_mean[s1];                           // The array storing average desirable sender trait values for errors
        double list_sd[s1];                             // The array storing standard deviations of desirable sender traits for errors
        double list_par_mean[s1];                       // The array storing average undesirable sender trait values for errors
        double list_par_sd[s1];                         // The array storing standaed deviations of undesirable sender traits for errors
        // Index
            int i_size=         0;
            int i_size_thr=     1;
            int i_patt=         2;
            int i_patt_thr=     3;
            int i_p_size=       4;  // phenotype
            int i_p_patt=       5;  // phenotype
        // undesirable sender params
            double par_size;
            double par_patt;
   	// Initialization of random number genertor
		int seed;
		dsfmt_t dsfmt;
		seed= time(NULL);
		if(seed==0)seed= 1;
		dsfmt_init_gen_rand(&dsfmt,seed);
	// Output
		FILE *out, *fpop, *fpop_par;
		out= fopen("summary.txt","w");
        fprintf(out,"gen\tsize\tsize_thr\tpatt\tpatt_thr\tpar_size\tpar_size_sd\tpar_patt\tpar_patt_sd\toff_size\tpar_size\ttyp1_err\ttyp2_err\n");
        fpop= fopen("final_pop.txt", "w");
        fprintf(fpop,"idvl\tg_size\tsize_thr\tg_patt\tpatt_thr\n");
        fpop_par= fopen("final_pop_par.txt","w");
        fprintf(fpop_par, "idvl\tg_size\tg_patt\n");
    // Initializing
        ini_check=reini_count=0;
        while (ini_check==0){
            tmp1= tmp2= 0;
            // Fill in the values
            for (i= 0; i< N; i++){
                Pop[i][i_size]=     rand_normal_dist(ini_size,ini_sd_size,dsfmt_genrand_open_open(&dsfmt),dsfmt_genrand_open_open(&dsfmt));
                Pop[i][i_size_thr]= rand_normal_dist(ini_size_thr,ini_sd_sthr,dsfmt_genrand_open_open(&dsfmt),dsfmt_genrand_open_open(&dsfmt));
                Pop[i][i_patt]=     rand_normal_dist(ini_patt,ini_sd_patt,dsfmt_genrand_open_open(&dsfmt),dsfmt_genrand_open_open(&dsfmt));
                Pop[i][i_patt_thr]= rand_normal_dist(ini_patt_thr,ini_sd_pthr,dsfmt_genrand_open_open(&dsfmt),dsfmt_genrand_open_open(&dsfmt));
            }
            for (i=0; i< N; i++){
                tmp1+= Pop[i][i_size_thr]/N;
                tmp2+= Pop[i][i_patt_thr]/N;
            }
            if ((tmp1/ini_size_thr)> 0.0 && (tmp2/ini_patt_thr)> 0.0) ini_check=1;
            reini_count+= 1;
            printf("%lf\t%lf\t%d\n", tmp1, tmp2, reini_count);
        }
        if(s3==0){
            par_size= ini_par_size;
            par_patt= ini_par_patt;
        }
        if(s3==1){
            for(i=0; i< N_par; i++){
                Pop_par[i][i_size]= rand_normal_dist(ini_par_size,ini_sd_psize,dsfmt_genrand_open_open(&dsfmt),dsfmt_genrand_open_open(&dsfmt));
                Pop_par[i][i_patt]= rand_normal_dist(ini_par_patt,ini_sd_ppatt,dsfmt_genrand_open_open(&dsfmt),dsfmt_genrand_open_open(&dsfmt));
        }}
    // Main loop
    for (i= 0; i< num_gen; i++){
        off_count= off_log= par_log= par_pass= off_pass= 0;
        // For each nest
            for(j=0; j< N; j++){
            // Produce offspring
                for(k=0; k< num_offspring; k++){
                    if(dsfmt_genrand_open_open(&dsfmt)< mut_rate){
                        temp_off[k][i_size_thr]= rand_normal_dist(Pop[j][i_size_thr],mut_size_thr, dsfmt_genrand_open_open(&dsfmt),dsfmt_genrand_open_open(&dsfmt));
                    }
                    else temp_off[k][i_size_thr]= Pop[j][i_size_thr];
                    if(dsfmt_genrand_open_open(&dsfmt)< mut_rate){
                        temp_off[k][i_patt_thr]= rand_normal_dist(Pop[j][i_patt_thr],mut_patt_thr, dsfmt_genrand_open_open(&dsfmt),dsfmt_genrand_open_open(&dsfmt));
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
                        if( (temp_par[i_p_size]> upr_size) || (temp_par[i_p_patt]> upr_patt) 
                            || (s5==1 && temp_par[i_p_size]<= buttom_thr) || (s5==1 && temp_par[i_p_patt]<= buttom_thr) ){
                            for(k=0; k< num_offspring; k++){
                                // If an offspring survives
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
                // Survival response based on the number of accepted offspring (regardless of its identity)
                if(s4==1){
                    upr_size= Pop[j][i_size_thr];
                    upr_patt= Pop[j][i_patt_thr];
                    nest_log= 0;
                    if(s1==1 && s2==1){
                        // Count the number of accepted offspring
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
                        nest_survival= nest_survive(nest_log);
                        // Survival processes of the undesirable sender
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
                        // Survival processes of the native offspring
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
                        // Count the number of accepted offspring
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
                        // Survival processes of the undesirable sender
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
                        // Survival processes of the native offspring
                        for(k=0; k< num_offspring; k++){
                            u1= dsfmt_genrand_open_open(&dsfmt);
                            if(upr_patt>= temp_off[k][i_p_patt]){
                                if(s5==1){
                                    if(temp_off[k][i_p_patt]> buttom_thr){
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
                        }}
                    }
                    if(s1==2){
                        // Count the number of accepted offspring
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
                        // Survival processes of the undesirable sender
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
                        // Survival processes of the native offspring
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
                    }}
                }}
            }
        // Creating the new population
            if(off_log>0){
                for(j=0; j< N; j++){
                    idx= round(dsfmt_genrand_open_open(&dsfmt)*(off_log-1));
                    for(k=0; k< num_geno; k++) Pop[j][k]= Offspring[idx][k];
                }
            }
            else i= num_gen;
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
                else i=num_gen;
            }
        // Erasing the offspring matrix
            for(j=0; j< num_offspring*N; j++){
                for(k=0; k< num_trait; k++) {
                    Offspring[j][k]= NAN;
                    Offs_log[j][k]= NAN;
                }
            }
            if(s3==1){
                for(j=0; j< N; j++){
                    for(k=0; k< num_trait; k++) Offs_par[j][k]= NAN;
                }
            }
        // Printing
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
                print_ppatt_sd= SD_array(temporary, N_par);
                for(j=0; j<N_par; j++) temporary[j]= Pop_par[j][i_patt];
                mean_par_patt= Mean_array(temporary, N_par);
                print_ppatt_sd= SD_array(temporary, N_par);
            }
            // Calculate the errors
            if(off_log>0){
                if(s4==0){
                    typ2_err= (double)par_log/N;
                    typ1_err= (double)(N-off_log/num_offspring-par_log)/N;
                }
                if(s4==1){
                    typ2_err= (double)par_pass/N;
                    typ1_err= (double)(N-off_pass/num_offspring)/N;
                }
                if(s3==0) fprintf(out,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%d\t%lf\t%lf\n", i, mean_size, mean_size_thr, mean_patt, mean_patt_thr, mean_par_size, mean_par_patt, off_log, par_log, typ1_err, typ2_err);
                if(s3==1 && par_log> 0) fprintf(out,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%d\t%d\t%lf\t%lf\n", i, mean_size, mean_size_thr, mean_patt, mean_patt_thr, mean_par_size, print_psize_sd, mean_par_patt, print_ppatt_sd, off_log, par_log, typ1_err, typ2_err);
            }
        // Printing the final population
            if(i>= (num_gen-1)){
                for(j=0; j<N; j++){
                    fprintf(fpop,"%d\t%lf\t%lf\t%lf\t%lf\n", j, Pop[j][i_size], Pop[j][i_size_thr], Pop[j][i_patt], Pop[j][i_patt_thr]);
                }
                for(j=0; j<N_par; j++){
                    fprintf(fpop_par,"%d\t%lf\t%lf\n", j, Pop_par[j][i_size], Pop_par[j][i_patt]);
                }
            }
            if (i== num_gen) printf("off_log is %d, par_log is %d.\n", off_log, par_log);
    }
    // Free memory
    fclose(out);
    fclose(fpop);
    fclose(fpop_par);
    return 0;
}
////////////////////// Functions are defined in below ////////////////////////
double Mean_array(double p[], int length){
    int i;
    double temp;
    temp= 0.0;
    if(length> 0){
        for (i=0; i<length; ++i) temp+= p[i]/(double)length;
        return temp;
    }
    else return NAN;
}

double SD_array(double p[], int length){
    if(length> 1){
        int i;
        double avg, temp;
        avg= temp= 0.0;
        for (i=0; i<length; ++i) avg+= p[i]/(double)length;
        for (i=0; i<length; ++i) temp+= pow(p[i]-avg, 2);
        return sqrt(temp/(length-1));
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
    if(sd>0){
        double z1;
        z1= (x-mu)/sqrt(2)/sd;
        return (1+erf(z1))/2;
    }
    else return NAN;
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