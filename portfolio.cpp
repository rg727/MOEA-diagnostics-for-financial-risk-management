
/* main-port-DPS.cpp

  Adapted by Andrew Hamilton from Lake Problem DPS, Jul 2018
  University of North Carolina at Chapel Hill
  andrew.hamilton@unc.edu

  Lake Problem DPS {
    Riddhi Singh, May, 2014
    The Pennsylvania State University
    rus197@psu.edu

    Adapted by Tori Ward, July 2014
    Cornell University
    vlw27@cornell.edu

    Adapted by Jonathan Herman and David Hadka, Sept-Dec 2014
    Cornell University and The Pennsylvania State University

    Adapted by Julianne Quinn, July 2015 as DPS problem
    Cornell University
    jdq8@cornell.edu
  }


*/


#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <sstream>
#include <ctime>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/moment.hpp>
#include <boost/accumulators/statistics/tail_quantile.hpp>
#include "boostutil.h"
//#include "borg.h"
#include "moeaframework.h"
//#include "mpi.h"

#define DPS_RUN_TYPE 2          // 0: 2dv version; 1: full DPS with RBFs, maxDebt formulation; 2: full DPS with RBFs, minRev formulation
#define BORG_RUN_TYPE 5      // 0: single run no borg; 1: borg run, serial; 2: borg parallel for cluster;
#define NUM_YEARS 20                   //20yr sims
#define NUM_SAMPLES 1000
#define NUM_DECISIONS_TOTAL 2           // each year, have to choose value snow contract + adjusted revenue
#define NUM_LINES_STOCHASTIC_INPUT 999999    //Input file samp.txt has 1M rows, 6 cols.
#define NUM_VARIABLES_STOCHASTIC_INPUT 6            //6 cols in input: swe,powIndex,revRetail,revWholesale,sswp,pswp
#define INDEX_STOCHASTIC_REVENUE 2   // 2 = revRetail, 3 = revWholesale
#define INDEX_STOCHASTIC_SNOW_PAYOUT 4    // 4 = sswp
#define INDEX_STOCHASTIC_POWER_INDEX 1  // 1 = power index
#if INDEX_STOCHASTIC_REVENUE == 2
#define MEAN_REVENUE 128.48255822159567     // mean revenue in absense of any financial risk mgmt. Make sure this is consistent with current input HHSamp revenue column.
#elif INDEX_STOCHASTIC_REVENUE == 3
#define MEAN_REVENUE  70.08967184742373     // mean revenue in absense of any financial risk mgmt. Make sure this is consistent with current input HHSamp revenue column.
#endif
#define NORMALIZE_SNOW_CONTRACT_SIZE 4.0
#define NORMALIZE_REVENUE 250.0
#define NORMALIZE_CONTINGENCY_FUND 150.0
#define NORMALIZE_POWER_PRICE 350.0
#define BUFFER_MAX_SIZE 5000
#define EPS 0.0000000000001
#define NUM_OBJECTIVES 4
#define EPS_ANNREV 0.075
#define EPS_MAXDEBT 0.225
#define EPS_MINREV 0.225
#define EPS_MAXCOMPLEXITY 0.05
#define EPS_MAXFUND 0.225
#if DPS_RUN_TYPE<2
#define NUM_CONSTRAINTS 1
#else
#define NUM_CONSTRAINTS 0
#endif
#define EPS_CONS1 0.05
#define NUM_RBF 4       // number of radial basis functions
#define NUM_INPUTS_RBF 3    // inputs: fund balance, power index, rev+hedge cash flow
#if DPS_RUN_TYPE>0
#define NUM_DV (2 * NUM_RBF * NUM_INPUTS_RBF) + (NUM_DECISIONS_TOTAL * (NUM_RBF + 2))
#else
#define NUM_DV 2
#endif
#define MIN_SNOW_CONTRACT 0.05          // DPS_RUN_TYPE==0 only: if contract slope dv < $0.05M/inch, act as if 0.
#define MIN_MAX_FUND 0.05               // DPS_RUN_TYPE==0 only: if max fund dv < $0.05M, act as if 0.
#define NUM_PARAM 6         // cost_fraction, discount_rate, delta_interest_fund, delta_interest_debt, lambda, lambda_prem_shift
#define NUM_PARAM_SAMPLES 1  // number of LHC samples in param file. Last line is values for SFPUC, Oct 2016.


namespace ublas = boost::numeric::ublas;
namespace tools = boost::math::tools;
namespace accumulator = boost::accumulators;
using namespace std;

void normalizeWeights(ublas::vector<double> & f_dv_w);
double policyAdjustedRevenue(const double f_fund_balance, const double f_debt, const double f_power_price_index,
                             const double f_cash_in, const ublas::vector<double> & f_dv_d,
                             const ublas::vector<double> & f_dv_c, const ublas::vector<double> & f_dv_b,
                             const ublas::vector<double> & f_dv_w, const ublas::vector<double> & f_dv_g);
double policySnowContractValue(const double f_fund_balance, const double f_debt, const double f_power_price_index,
                               const ublas::vector<double> & f_dv_d,
                               const ublas::vector<double> & f_dv_c, const ublas::vector<double> & f_dv_b,
                               const ublas::vector<double> & f_dv_w, const ublas::vector<double> & f_dv_g) ;
double policyAdjustedRevenue_2dv(const double f_fund_balance, const double f_cash_in, const double f_max_fund_size);
double policySnowContractValue_2dv(const double f_value);
double policyMaxFund_2dv(const double f_value);

double stochastic_input[NUM_LINES_STOCHASTIC_INPUT][NUM_VARIABLES_STOCHASTIC_INPUT];                   // Stochastic variables
double param_LHC_sample[NUM_PARAM][NUM_PARAM_SAMPLES];
#if (BORG_RUN_TYPE == 0)|(BORG_RUN_TYPE == 5)
double problem_dv[NUM_DV];
double problem_objs[NUM_OBJECTIVES];
double problem_constraints[NUM_CONSTRAINTS];
double pareto[NUM_DV + 2*(NUM_OBJECTIVES + NUM_CONSTRAINTS)][10000];
double N_pareto = 0;
#endif

ublas::vector<double> annualized_adjusted_revenue(NUM_SAMPLES);     // objectives
ublas::vector<double> max_hedge_complexity(NUM_SAMPLES);
ublas::vector<double> max_fund_balance(NUM_SAMPLES);
#if DPS_RUN_TYPE<2
ublas::vector<double> debt_steal(NUM_SAMPLES);
#else
ublas::vector<double> min_adjusted_revenue(NUM_SAMPLES);
#endif

ublas::vector<double> revenue(NUM_YEARS);                          // state variables
ublas::vector<double> payout_snow_contract(NUM_YEARS);      // snow contract (sswp or sput)
ublas::vector<double> power_price_index(NUM_YEARS+1);             // power price index
ublas::vector<double> discount_factor(NUM_YEARS);

ublas::vector<double> adjusted_revenue(NUM_YEARS);                     // decisions
ublas::vector<double> fund_balance(NUM_YEARS + 1);
ublas::vector<double> fund_withdrawal(NUM_YEARS);
ublas::vector<double> debt(NUM_YEARS + 1);
ublas::vector<double> value_snow_contract(NUM_YEARS);           // value of contract 1
ublas::vector<double> dv_d(NUM_DECISIONS_TOTAL);                 // decision variables, dps params
ublas::vector<double> dv_c(NUM_RBF * NUM_INPUTS_RBF);
ublas::vector<double> dv_b(NUM_RBF * NUM_INPUTS_RBF);
ublas::vector<double> dv_w(NUM_RBF * NUM_DECISIONS_TOTAL);
ublas::vector<double> dv_g(NUM_DECISIONS_TOTAL);

double cost_fraction;   // params to be looped over with LHC sample, doing borg each time. Declare globally.
double avg_surplus_revenue;
double discount_rate;
double interest_fund;
double interest_debt;
double lambda_prem_shift;
int LHC_set;
//int NFE;
//int NFE_counter = 0;
//unsigned int seed_borg;
unsigned int seed_sample;      // use same seed for each function evaluation, so always comparing same simulations. should be less noisy.
int lines_to_use[NUM_SAMPLES];

#if DPS_RUN_TYPE<2
typedef accumulator::accumulator_set<double, accumulator::stats<accumulator::tag::tail_quantile<accumulator::right> > > accumulator_t;
#endif

// problem for borg search
void portfolioProblem(double *problem_dv, double *problem_objs, double *problem_constraints) {
    //NFE_counter += 1;
//    if ((NFE_counter % 100) == 0){printf("%d\n", NFE_counter);}

    // initialize variables
    zero(annualized_adjusted_revenue);
    zero(max_hedge_complexity);
    zero(max_fund_balance);
    zero(revenue);
    zero(power_price_index);
    zero(payout_snow_contract);
    zero(discount_factor);
    zero(dv_d);
    zero(dv_c);
    zero(dv_b);
    zero(dv_w);
    zero(dv_g);
#if DPS_RUN_TYPE<2
    zero(debt_steal);
#else
    zero(min_adjusted_revenue);
#endif

#if DPS_RUN_TYPE>0
    // get dvs
    for (int i = 0; i < dv_d.size(); i++){
        dv_d(i) = problem_dv[i];                      // cutoffs: dv_d = [CFMAX, XMIN1, XMIN2]
    }
    for (int i = 0; i < dv_c.size(); i++){
        dv_c(i) = problem_dv[i + dv_d.size()];        // centers for RBFs
    }
    for (int i = 0; i < dv_b.size(); i++){
        dv_b(i) = max(EPS, problem_dv[i + dv_d.size() + dv_c.size()]);          // radii for RBFs
    }
    for (int i = 0; i < dv_w.size(); i++){
        dv_w(i) = problem_dv[i + dv_d.size() + dv_c.size() + dv_b.size()];        // weights for RBFs
    }
    for (int i = 0; i < dv_g.size(); i++){
        dv_g(i) = problem_dv[i + dv_d.size() + dv_c.size() + dv_b.size() + dv_w.size()];                      // const addition
    }
    // normalize weights
    normalizeWeights(dv_w);
#else
    double fixed_max_fund = policyMaxFund_2dv(problem_dv[0]);
    double fixed_value_snow_contract = policySnowContractValue_2dv(problem_dv[1]);
#endif

    // create discounting factor
    double discount_normalization = 0.0;      // discounting normalization, 1/sum_(discount_factor)
    for (int i = 0; i < NUM_YEARS; i++){
        discount_factor(i) = pow(discount_rate, i+1);
        discount_normalization += discount_factor(i);
    }
    discount_normalization = 1.0 / discount_normalization;

#if DPS_RUN_TYPE<2
    accumulator_t debt_q95(accumulator::tag::tail<accumulator::right>::cache_size = NUM_SAMPLES);    // accumulator object for calculating upper 5th quantile of debt
#endif

    double net_payout_snow_contract = 0.;
    double cash_in = 0.;

    // run revenue model simulation
    for (int s = 0; s < NUM_SAMPLES; s++) {
        // randomly generated revenues
#if NUM_SAMPLES > 1
        int index = lines_to_use[s];
#else
        int index = 1;
#endif
//        printf("\n\nSample %d  %d\n", s, index);

        // get the random revenue from the States of the world file
        //each line of SOW file covers 20 years of revenue
        power_price_index(0) = stochastic_input[index - 1][INDEX_STOCHASTIC_POWER_INDEX];
        for (int i = 0; i < NUM_YEARS; i++) {
            revenue(i) = (stochastic_input[index + i][INDEX_STOCHASTIC_REVENUE] - MEAN_REVENUE * cost_fraction);
            payout_snow_contract(i) = stochastic_input[index + i][INDEX_STOCHASTIC_SNOW_PAYOUT];
            power_price_index(i+1) = stochastic_input[index + i][INDEX_STOCHASTIC_POWER_INDEX];
//            printf("%f  %f  %f\n", revenue(i), payout_snow_contract(i), payout_power_contract(i));
        }

        // initial simulation variables
        zero(value_snow_contract);
        zero(fund_balance);
        zero(fund_withdrawal);
        zero(adjusted_revenue);
        zero(debt);

        max_hedge_complexity(s) = 0;
#if DPS_RUN_TYPE==0
        if (fixed_value_snow_contract > EPS) {
            max_hedge_complexity(s) = 1;
        }
#endif

        //calculate new revenues, contingency fund balance, objectives
        for (int i = 0; i < NUM_YEARS; i++) {

            // find next policy-derived index insurance and CF withdrawal
#if  DPS_RUN_TYPE>0
            value_snow_contract(i) = policySnowContractValue(fund_balance(i), debt(i),
                                                             power_price_index(i), dv_d, dv_c, dv_b, dv_w, dv_g);
            net_payout_snow_contract = value_snow_contract(i) * (payout_snow_contract(i) - lambda_prem_shift);
            cash_in = revenue(i) + net_payout_snow_contract - (debt(i) * interest_debt);
            adjusted_revenue(i) = policyAdjustedRevenue(fund_balance(i) * interest_fund, debt(i) * interest_debt,
                                                        power_price_index(i + 1), cash_in,
                                                        dv_d, dv_c, dv_b, dv_w, dv_g);
            if (abs(value_snow_contract(i)) > EPS) {
                max_hedge_complexity(s) = 1;
            }

#else
            net_payout_snow_contract = fixed_value_snow_contract * (payout_snow_contract(i) - lambda_prem_shift);
            cash_in = revenue(i) + net_payout_snow_contract - (debt(i) * interest_debt);
            adjusted_revenue(i) = policyAdjustedRevenue_2dv(fund_balance(i) * interest_fund, cash_in, fixed_max_fund);
#endif
            fund_withdrawal(i) = adjusted_revenue(i) - cash_in;
#if DPS_RUN_TYPE<2
            if (adjusted_revenue(i) < 0.){
                debt(i + 1) = -adjusted_revenue(i);
                adjusted_revenue(i) = 0;
            }
#endif
            fund_balance(i + 1) = fund_balance(i) * interest_fund - fund_withdrawal(i);

            annualized_adjusted_revenue(s) += adjusted_revenue(i) * discount_factor(i);

        }
        annualized_adjusted_revenue(s) = discount_normalization *
                                         (annualized_adjusted_revenue(s) +
                                          ((fund_balance(NUM_YEARS) * interest_fund * discount_factor(0)) -
                                           (debt(NUM_YEARS) * interest_debt * discount_factor(0))) *
                                          discount_factor(NUM_YEARS - 1));
#if DPS_RUN_TYPE<2
        debt_q95(vmax(debt));                                                               //q95(max(debt)) constraint
        debt_steal(s) = debt(NUM_YEARS) - debt(NUM_YEARS - 1);
#else
        min_adjusted_revenue(s) = vmin(adjusted_revenue);
#endif

        max_fund_balance(s) = vmax(fund_balance);

    }

    // aggregate objectives
    problem_objs[0] = -1 * vsum(annualized_adjusted_revenue) / NUM_SAMPLES; // max: average annualized adjusted_revenue, across samp
#if DPS_RUN_TYPE<2
    problem_objs[1] = accumulator::quantile(debt_q95, accumulator::quantile_probability = 0.95);    //minimize 95th percentile of max debt
#else
    problem_objs[1] = -1 * vsum(min_adjusted_revenue) / NUM_SAMPLES;
#endif
    problem_objs[2] = 1 * vsum(max_hedge_complexity) / NUM_SAMPLES; // min: avg_avg_hedging complexity
    problem_objs[3] = 1 * vsum(max_fund_balance) / NUM_SAMPLES; // min: max_fund_balance

#if DPS_RUN_TYPE<2
    // check constraint
    problem_constraints[0] = max(0.0, (vsum(debt_steal) / NUM_SAMPLES) - EPS_CONS1);
#endif

//    printf("\n\n\n%f   %f   %f   %f   %f\n\n\n", problem_objs[0], problem_objs[1], problem_objs[2], problem_objs[3], problem_constraints[0]);

    annualized_adjusted_revenue.clear();
    max_hedge_complexity.clear();
    max_fund_balance.clear();
    fund_balance.clear();
    revenue.clear();
    power_price_index.clear();
    payout_snow_contract.clear();
    discount_factor.clear();
    value_snow_contract.clear();
    fund_withdrawal.clear();
    adjusted_revenue.clear();
    dv_d.clear();
    dv_c.clear();
    dv_b.clear();
    dv_w.clear();
    dv_g.clear();
#if DPS_RUN_TYPE<2
    debt_steal.clear();
    debt.clear();
#else
    min_adjusted_revenue.clear();
#endif

}

// normalize RBF weights. for each decision, weights should sum to 1.
void normalizeWeights(ublas::vector<double> & f_dv_w){
    for (int j = 0; j < NUM_DECISIONS_TOTAL; j++){
        double total = 0.0;
        for(int i = j * NUM_RBF; i < (j + 1) * NUM_RBF; i++) {
            total += f_dv_w(i);
        }
        if (total != 0){
            for (int i = j * NUM_RBF; i < (j + 1) * NUM_RBF; i++){
                f_dv_w(i) = f_dv_w(i) / total;
//                printf("%f   ",dv_w(i));
            }
        }
    }

}

// calculate adjusted revenue at end of year, using fund balance, power price, and cash flow as inputs.
double policyAdjustedRevenue(const double f_fund_balance, const double f_debt, const double f_power_price_index,
                             const double f_cash_in, const ublas::vector<double> & f_dv_d,
                             const ublas::vector<double> & f_dv_c, const ublas::vector<double> & f_dv_b,
                             const ublas::vector<double> & f_dv_w, const ublas::vector<double> & f_dv_g) {
    // decision between 0 and 1 based on RBF and state variable inputs
    int decision_order = 0;
    double value = 0;
    for (int i = 0; i < NUM_RBF; i++){
        value =+ f_dv_w(i) * exp(
                - pow((f_fund_balance - f_debt + NORMALIZE_CONTINGENCY_FUND)/(2 * NORMALIZE_CONTINGENCY_FUND)
                      - f_dv_c(NUM_INPUTS_RBF * i), 2) / f_dv_b(NUM_INPUTS_RBF * i)
                - pow(f_power_price_index/NORMALIZE_POWER_PRICE - f_dv_c(NUM_INPUTS_RBF * i + 1), 2)
                  / f_dv_b(NUM_INPUTS_RBF * i + 1)
                - pow((f_cash_in + NORMALIZE_REVENUE)/(2 * NORMALIZE_REVENUE) - f_dv_c(NUM_INPUTS_RBF * i + 2), 2)
                  / f_dv_b(NUM_INPUTS_RBF * i + 2) ) ;
//        printf("%f  %f  %f  %f  %f", f_dv_w(i), f_dv_c(2*i), f_dv_b(2*i), f_dv_c(2*i+1), f_dv_b(2*i+1));
//        printf("%f  ", value);
    }
    value =+ f_dv_g(decision_order);
    // now scale back to [-NORMALIZE_REVENUE,NORMALIZE_REVENUE]
    value = (value * 2 * NORMALIZE_REVENUE) - NORMALIZE_REVENUE;
    // ensure that cant withdraw more than fund balance
    value = min(value , f_cash_in + f_fund_balance);
    // ensure that cant deposit more than space left in fund balance
    value = max(value, f_cash_in + f_fund_balance - f_dv_d(decision_order) * NORMALIZE_CONTINGENCY_FUND);
    // ensure that cant deposit more than available cash flow
    value = max(value, f_cash_in - max(f_cash_in, 0.));
//    printf("%f    %f    %f\n", f_fund_balance, f_cash_in, value);
    return value;
}

// calculate value multiplier for snow contract this year. input = fund balance.
double policySnowContractValue(const double f_fund_balance, const double f_debt,  const double f_power_price_index,
                               const ublas::vector<double> & f_dv_d,
                               const ublas::vector<double> & f_dv_c, const ublas::vector<double> & f_dv_b,
                               const ublas::vector<double> & f_dv_w, const ublas::vector<double> & f_dv_g)  {
    int decision_order = 1;
    double value = 0;
    double dum = f_dv_d[0];
    // decision between 0 and 1 based on RBF and state var inputs
    for (int i = 0; i < NUM_RBF; i++){
        value =+ f_dv_w(decision_order * NUM_RBF + i) * exp(
                - pow((f_fund_balance - f_debt + NORMALIZE_CONTINGENCY_FUND)/(2 * NORMALIZE_CONTINGENCY_FUND)
                      - f_dv_c(NUM_INPUTS_RBF * i), 2) / f_dv_b(NUM_INPUTS_RBF * i)
                - pow(f_power_price_index/NORMALIZE_POWER_PRICE - f_dv_c(NUM_INPUTS_RBF * i + 1), 2)
                  / f_dv_b(NUM_INPUTS_RBF * i + 1));
//        printf("%f\n", value);
    }
    // scale back to [0, NORMALIZE_SNOW_CONTRACT_SIZE]
    value = min((value + f_dv_g(decision_order)) * NORMALIZE_SNOW_CONTRACT_SIZE, NORMALIZE_SNOW_CONTRACT_SIZE);
    // enforce minimum contract size
    if (value < f_dv_d(decision_order) * NORMALIZE_SNOW_CONTRACT_SIZE){
        value = 0.0;
    }
//    printf("\n");

    return value;
}

// calculate 2dv version (DPS_RUN_TYPE==0) adjusted revenue at end of year, using fund balance, cash flow, and max fund size as inputs.
double policyAdjustedRevenue_2dv(const double f_fund_balance, const double f_cash_in, const double f_max_fund_size) {
    double adjusted_revenue = 0;
    if (f_cash_in < 0.0){
        if (f_fund_balance < EPS){
            adjusted_revenue = f_cash_in;
        }else{
            adjusted_revenue = min(f_cash_in + f_fund_balance, 0.0);
        }
    }else{
        if (f_fund_balance > (f_max_fund_size - EPS)){
            adjusted_revenue = f_cash_in + (f_fund_balance - f_max_fund_size);
        }else{
            adjusted_revenue = max(f_cash_in - (f_max_fund_size - f_fund_balance), 0.0);
        }
    }
    return adjusted_revenue;
}

// calculate 2dv version (DPS_RUN_TYPE==0) snow contract value, ensuring >= MIN_SNOW_CONTRACT
double policySnowContractValue_2dv(const double f_value) {
    double value_edit;
    if (f_value < MIN_SNOW_CONTRACT){
        value_edit = 0.;
    }else{
        value_edit = f_value;
    }
    return value_edit;
}

// calculate 2dv version (DPS_RUN_TYPE==0) snow contract value, ensuring >= MIN_SNOW_CONTRACT
double policyMaxFund_2dv(const double f_value) {
    double value_edit;
    if (f_value < MIN_MAX_FUND){
        value_edit = 0.;
    }else{
        value_edit = f_value;
    }
    return value_edit;
}

int main(int argc, char* argv[]) {

    clock_t begin = clock();

    // get stochastic inputs
    for (int i = 0; i < NUM_LINES_STOCHASTIC_INPUT; i++) {
        for (int j = 0; j < NUM_VARIABLES_STOCHASTIC_INPUT; j++) {
            stochastic_input[i][j] = 0.0;
        }
    }

    FILE *myfile;
    myfile = fopen("./HHsamp06262019.txt", "r");

    int linenum = 0;
    char testbuffer[BUFFER_MAX_SIZE];

    if (myfile == NULL) {
        perror("Error opening file");
    } else {
        char buffer[BUFFER_MAX_SIZE];
        fgets(buffer, BUFFER_MAX_SIZE, myfile);       // eat header line
        while (fgets(buffer, BUFFER_MAX_SIZE, myfile) != NULL) {
            linenum++;
            if (buffer[0] != '#') {
                char *pStart = testbuffer;
                char *pEnd;
                for (int i = 0; i < BUFFER_MAX_SIZE; i++) {
                    testbuffer[i] = buffer[i];
                }
                for (int cols = 0; cols < NUM_VARIABLES_STOCHASTIC_INPUT; cols++) {
                    stochastic_input[linenum - 1][cols] = strtod(pStart, &pEnd);
                    pStart = pEnd;
//                    printf("%f ",stochastic_input[linenum-1][cols]);
                }
//                printf("\n");
            }
        }

    }
    fclose(myfile);

    // read in LHC dv
    FILE *myfile2;
    myfile2 = fopen("./param_SFPUC_withLamPremShift.txt", "r");
    char testbuffer2[BUFFER_MAX_SIZE];
    linenum = 0;

    if (myfile2 == NULL) {
        perror("Error opening file");
    } else {
        char buffer2[BUFFER_MAX_SIZE];
        fgets(buffer2, BUFFER_MAX_SIZE, myfile2);       // eat header line
        while (fgets(buffer2, BUFFER_MAX_SIZE, myfile2) != NULL) {
            linenum++;
            if (buffer2[0] != '#') {
                char *pStart2 = testbuffer2;
                char *pEnd2;
                for (int i = 0; i < BUFFER_MAX_SIZE; i++) {
                    testbuffer2[i] = buffer2[i];
                }
                for (int i = 0; i < NUM_PARAM; i++) {
                    param_LHC_sample[i][linenum - 1] = strtod(pStart2, &pEnd2);
                    pStart2 = pEnd2;
//                    if (linenum == 1){
//                        printf("%f  %d\n",param_LHC_sample[i][linenum-1], i);
//                    }
                }
            }
        }
    }
//    printf("\n");
    fclose(myfile2);

#if BORG_RUN_TYPE > 0
    // setting random seeds
    //seed_borg = atoi(argv[1]);
    seed_sample = 1;      // use same seed for sample each time, so always comparing same simulations
    //NFE = atoi(argv[3]);
#else
    seed_sample = atoi(argv[1]);
#endif

    srand(seed_sample);
    for (int s = 0; s < NUM_SAMPLES; s++) {
        //choose NUM_SAMPLES number of samples (starting year out of NUM_YEARS). Cant be 0, since need power_price_index[t-1], and cant be less than NUM_YEARS from end
        lines_to_use[s] = rand() % (NUM_LINES_STOCHASTIC_INPUT - NUM_YEARS - 1) + 1;
//        printf("%d\n", lines_to_use[s]);
    }

#if BORG_RUN_TYPE == 2
    // interface with Borg-MS
    BORG_Algorithm_ms_startup(&argc, &argv);
    BORG_Algorithm_ms_max_evaluations(NFE);
    BORG_Algorithm_output_frequency(NFE / 200);
#endif

#if BORG_RUN_TYPE == 0      // get dv and params from argv, not borg, and just run once
    LHC_set = atoi(argv[2]);

    // read snow contract put strike premiums
    FILE *myfile4;
    myfile4 = fopen(argv[3], "r");
    char testbuffer4[BUFFER_MAX_SIZE];
    linenum = 0;

    if (myfile4 == NULL) {
        perror("Error opening file");
    } else {
        char buffer4[BUFFER_MAX_SIZE];
//        fgets(buffer4, BUFFER_MAX_SIZE, myfile4);       // eat header line
        while (fgets(buffer4, BUFFER_MAX_SIZE, myfile4) != NULL) {
            linenum++;
            if (buffer4[0] != '#') {
                char *pStart4 = testbuffer4;
                char *pEnd4;
                for (int i = 0; i < BUFFER_MAX_SIZE; i++) {
                    testbuffer4[i] = buffer4[i];
                }
                for (int i = 0; i < (NUM_DV + NUM_OBJECTIVES + NUM_CONSTRAINTS); i++) {
                    pareto[i][linenum - 1] = strtod(pStart4, &pEnd4);
                    pStart4 = pEnd4;
                }
            }
        }
    }
//    printf("\n");
    fclose(myfile4);
    N_pareto = linenum;

    // read snow contract put strike premiums
    ofstream retest_write;
    retest_write.open(argv[4], ios::out | ios::trunc);

    for (int i = 0; i < N_pareto; i++){
        // decision variables
        for (int j = 0; j < NUM_DV; j++) {
            problem_dv[j] = pareto[j][i];
        }

        // params from LHC sample
        cost_fraction = param_LHC_sample[0][LHC_set];             // fraction of MEAN_REVENUE that is must-meet costs
        double delta = param_LHC_sample[1][LHC_set];              // discount rate, as %/yr
        double Delta_interest_fund = param_LHC_sample[2][LHC_set];    // interest rate on reserve funds, as %/yr, markdown below delta (all negative)
        double Delta_interest_debt = param_LHC_sample[3][LHC_set];    // interest rate charged on debt, as %/yr, markup above delta (all positive)
        lambda_prem_shift = param_LHC_sample[5][LHC_set];         // shift in snow contract premium to apply, based on lambda parameter (relative to based premiums from lambda=0.25)

        // calculated params for LHC sensitivity analysis, used in portfolioProblem
        avg_surplus_revenue = MEAN_REVENUE * (1. - cost_fraction);
        discount_rate = 1. / (delta / 100. + 1.);
        interest_fund = (Delta_interest_fund + delta) / 100. + 1.;
        interest_debt = (Delta_interest_debt + delta) / 100. + 1.;

        // run dps with given dv
        portfolioProblem(problem_dv, problem_objs, problem_constraints);

        retest_write << std::fixed << std::setprecision(10);
        for (int j = 0; j < (NUM_DV + NUM_OBJECTIVES + NUM_CONSTRAINTS); ++j) {
            retest_write << pareto[j][i] << "\t";
        }
        for (int j = 0; j < NUM_OBJECTIVES; ++j) {
            retest_write << problem_objs[j] << "\t";
        }
        for (int j = 0; j < NUM_CONSTRAINTS; ++j) {
            retest_write << problem_constraints[j] << "\t";
        }
        retest_write << '\n';

    }
    retest_write.close();


#elif BORG_RUN_TYPE > 0
    // loop over uncertain parameters in LHC sample, borg each time

    for (int p = 0; p < NUM_PARAM_SAMPLES; ++p){
        LHC_set = p;
        // params from LHC sample
        cost_fraction = param_LHC_sample[0][p];             // fraction of MEAN_REVENUE that is must-meet costs
        double delta = param_LHC_sample[1][p];              // discount rate, as %/yr
        double Delta_interest_fund = param_LHC_sample[2][p];    // interest rate on reserve funds, as %/yr, markdown below delta (all negative)
        double Delta_interest_debt = param_LHC_sample[3][p];    // interest rate charged on debt, as %/yr, markup above delta (all positive)
        lambda_prem_shift = param_LHC_sample[5][p];         // shift in snow contract premium to apply, based on lambda parameter (relative to based premiums from lambda=0.25)


        // calculated params for LHC sensitivity analysis, used in portfolioProblem
        avg_surplus_revenue = MEAN_REVENUE * (1. - cost_fraction);
        discount_rate = 1. / (delta / 100. + 1.);
        interest_fund = (Delta_interest_fund + delta) / 100. + 1.;
        interest_debt = (Delta_interest_debt + delta) / 100. + 1.;
        }

        // Define the problem with decisions, objectives, constraints and the evaluation function
        //BORG_Problem problem = BORG_Problem_create(NUM_DV, NUM_OBJECTIVES, NUM_CONSTRAINTS, portfolioProblem);

//#if DPS_RUN_TYPE>0
//        // Set all the parameter bounds and epsilons
//        for (int i = 0; i < dv_d.size(); i++){
//            BORG_Problem_set_bounds(problem, i, 0, 1);            //threshold params  (dv_d)
//        }
//        for (int i = 0; i < dv_c.size(); i++){
//            BORG_Problem_set_bounds(problem, dv_d.size() + i, -1, 1);             // dv_c
//        }
//        for (int i = 0; i < dv_b.size(); i++){
//            BORG_Problem_set_bounds(problem, dv_d.size() + dv_c.size() + i, EPS, 1);             // dv_b (needs to be > 0 or numerical issues)
//        }
//        for (int i = 0; i < dv_w.size(); i++){
//            BORG_Problem_set_bounds(problem, dv_d.size() + dv_c.size() + dv_b.size() + i, 0, 1);             // dv_w
//        }
//        for (int i = 0; i < dv_g.size(); i++){
//            BORG_Problem_set_bounds(problem, dv_d.size() + dv_c.size() + dv_b.size() + dv_w.size() + i, 0, 1);            //const addition (dv_g)
//        }
//#else
//        BORG_Problem_set_bounds(problem, 0, 0, NORMALIZE_CONTINGENCY_FUND);            //bounds for contingency max size
//        BORG_Problem_set_bounds(problem, 1, 0, NORMALIZE_SNOW_CONTRACT_SIZE);            //bounds for swap contract slope
//#endif
//
//        BORG_Problem_set_epsilon(problem, 0, EPS_ANNREV); // avg_annualized_adjusted_revenue (units $M, so 0.01=$10,000)
//#if DPS_RUN_TYPE<2
//        BORG_Problem_set_epsilon(problem, 1, EPS_MAXDEBT); // q95_max_debt (units fraction of avg_surplus_revenue, so 0.02=2%)
//#else
//        BORG_Problem_set_epsilon(problem, 1, EPS_MINREV); // q95_max_debt (units fraction of avg_surplus_revenue, so 0.02=2%)
//#endif
//        BORG_Problem_set_epsilon(problem, 2, EPS_MAXCOMPLEXITY); // max_hedge_complexity (unitless)
//        BORG_Problem_set_epsilon(problem, 3, EPS_MAXFUND); // max_fund_balance (units $M)
//
//#if BORG_RUN_TYPE == 1
//
//        //This is set up to run only one seed at a time
//        char output_filename[256];
//        FILE *output_file = NULL;
//#if DPS_RUN_TYPE==1
//        sprintf(output_filename, "./sets/PortDPS_DPS_maxDebt.set");
//#elif DPS_RUN_TYPE==2
//        sprintf(output_filename, "./sets/PortDPS_DPS_minRev.set");
//#else
//        sprintf(output_filename, "./sets/PortDPS_2dv.set");
//#endif
//        BORG_Random_seed(seed_borg);
//        BORG_Archive result = BORG_Algorithm_run(problem, NFE); // this actually runs the optimization
//
//        //If this is the master node, print out the final archive
//        if (result != NULL) {
//            output_file = fopen(output_filename, "w");
//            if (!output_file) {
//                BORG_Debug("Unable to open final output file\n");
//            }
//            BORG_Archive_print(result, output_file);
//            BORG_Archive_destroy(result);
//            fclose(output_file);
//        }
//
//        BORG_Problem_destroy(problem);
//
//#elif BORG_RUN_TYPE == 2
//        //This is set up to run only one seed at a time
//        char output_filename[256];
//        char runtime[256];
//        FILE *output_file = NULL;
//#if DPS_RUN_TYPE==1
//        sprintf(output_filename, "./sets/PortDPS_DPS_maxDebt_samp%d_seedS%d_seedB%d.set", NUM_SAMPLES, seed_sample, seed_borg);
//        sprintf(runtime, "./runtime/PortDPS_DPS_maxDebt_samp%d_seedS%d_seedB%d.runtime", NUM_SAMPLES, seed_sample, seed_borg);
//#elif DPS_RUN_TYPE==2
//        sprintf(output_filename, "./sets/PortDPS_DPS_minRev_samp%d_seedS%d_seedB%d.set", NUM_SAMPLES, seed_sample, seed_borg);
//        sprintf(runtime, "./runtime/PortDPS_DPS_minRev_samp%d_seedS%d_seedB%d.runtime", NUM_SAMPLES, seed_sample, seed_borg);
//#else
//        sprintf(output_filename, "./sets/PortDPS_2dv_samp%d_seedS%d_seedB%d.set", NUM_SAMPLES, seed_sample, seed_borg);
//        sprintf(runtime, "./runtime/PortDPS_2dv_samp%d_seedS%d_seedB%d.runtime", NUM_SAMPLES, seed_sample, seed_borg);
//#endif
//        BORG_Algorithm_output_runtime(runtime);
//
//        BORG_Random_seed(seed_borg);
//        BORG_Archive result = BORG_Algorithm_ms_run(problem); // this actually runs the optimization
//
//        //If this is the master node, print out the final archive
//        if (result != NULL) {
//            output_file = fopen(output_filename, "w");
//            if (!output_file) {
//                BORG_Debug("Unable to open final output file\n");
//            }
//            BORG_Archive_print(result, output_file);
//            BORG_Archive_destroy(result);
//            fclose(output_file);
//        }
//
//        BORG_Problem_destroy(problem);
//
#endif
//
//        clock_t end = clock();
//        double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC / 60.0;
////        printf("%d\t%f\n", p, elapsed_secs);
//    }
//#endif
//
//#if (BORG_RUN_TYPE==2)
//    BORG_Algorithm_ms_shutdown();
//#endif

    MOEA_Init(NUM_OBJECTIVES, NUM_CONSTRAINTS);

      while (MOEA_Next_solution() == MOEA_SUCCESS) {
        MOEA_Read_doubles(NUM_DV, problem_dv);
        portfolioProblem(problem_dv, problem_objs, problem_constraints);
        MOEA_Write(problem_objs, problem_constraints);
    }
    MOEA_Terminate();

    return EXIT_SUCCESS;

}


