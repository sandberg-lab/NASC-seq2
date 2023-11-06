#include "arb_hypgeom.h"
#include "arb.h"
#include "arf.h"
#include "string.h"
#include "flint/fmpz.h"

void get_hyp1f1_first_1(arb_t res, const arb_t kon, const arb_t koff, const arb_t ksyn, const arb_t t, slong n, slong prec)
{
	arb_t a,b,z;
        arb_init(a);arb_init(b);arb_init(z);
        
	arb_neg(a, kon);
	arb_add_si(a,a,n,prec);

        arb_add_si(b,b,1,prec);
        arb_sub(b,b,kon,prec);
        arb_sub(b,b,koff,prec);
	arb_add_si(b,b,n,prec);

        arb_neg(z,t);
	arb_exp(z,z,prec);
        arb_mul(z, ksyn,z,prec);

        arb_hypgeom_1f1(res,a,b,z,0,prec);
	arb_clear(a);arb_clear(b);arb_clear(z);
}

void get_hyp1f1_first_2(arb_t res, const arb_t kon, const arb_t koff, const arb_t ksyn, slong n, slong prec)
{
	arb_t a,b,z;
        arb_init(a);arb_init(b);arb_init(z);
        arb_set(a, kon);
	arb_add_si(a,a,n,prec);

        arb_add(b,b,kon,prec);
        arb_add(b,b,koff,prec);
	arb_add_si(b,b,n,prec);

        arb_neg(z,ksyn);

        arb_hypgeom_1f1(res,a,b,z,0,prec);
	arb_clear(a);arb_clear(b);arb_clear(z);
}
void get_hyp1f1_second_1(arb_t res, const arb_t kon, const arb_t koff, const arb_t ksyn,const arb_t t, slong n, slong prec)
{
        arb_t a,b,z;
        arb_init(a);arb_init(b);arb_init(z);
        
	arb_set(a, koff);
	arb_add_si(a,a,n,prec);

	arb_set_si(b,1);
        arb_add(b,b,kon,prec);
	arb_add(b,b,koff,prec);
	arb_add_si(b,b,n,prec);

        arb_neg(z,t);
        arb_exp(z,z,prec);
        arb_mul(z, ksyn,z,prec);

        arb_hypgeom_1f1(res,a,b,z,0,prec);
        arb_clear(a);arb_clear(b);arb_clear(z);
}
void get_hyp1f1_second_2(arb_t res, const arb_t kon, const arb_t koff, const arb_t ksyn, slong n, slong prec)
{
        arb_t a,b,z;
        arb_init(a);arb_init(b);arb_init(z);

	arb_sub(a,a,koff,prec);
        arb_add_si(a,a,n,prec);

        arb_set_si(b,1);
        arb_sub(b,b,kon,prec);
        arb_sub(b,b,koff,prec);
        arb_add_si(b,b,n,prec);

        arb_neg(z,ksyn);

        arb_hypgeom_1f1(res,a,b,z,0,prec);
}

void get_coef_first(arb_t res, const arb_t ksyn, ulong n, slong prec)
{
	arb_t numer, denom;
	arb_init(numer);arb_init(denom);

	arb_pow_ui(numer, ksyn, n, prec);

	arb_fac_ui(denom, n, prec);

	arb_div(res, numer, denom, prec);

	arb_clear(numer);arb_clear(denom);
}

void get_factors_first(arb_t res1, arb_t res2, const arb_t kon, const arb_t koff, const arb_t t, ulong n, slong prec)
{
	arb_t numer1, numer2, denom1, denom2;
	arb_init(numer1);arb_init(numer2);arb_init(denom1);arb_init(denom2);

	
	arb_t numer1_part1, numer1_part2;
	arb_init(numer1_part1);arb_init(numer1_part2);

	arb_set(numer1_part1, kon);
	arb_neg(numer1_part1, numer1_part1);
	arb_rising_ui(numer1_part1, numer1_part1, n, prec);

	arb_set(numer1_part2, t);
	arb_mul_ui(numer1_part2, numer1_part2, n, prec);
	arb_neg(numer1_part2, numer1_part2);
	arb_exp(numer1_part2, numer1_part2, prec);

	arb_mul(numer1, numer1_part1, numer1_part2, prec);
	if ( n % 2 == 1)
	{
		arb_neg(numer1, numer1);
	}


	arb_add(denom1, kon, koff, prec);
	arb_neg(denom1, denom1);
	arb_add_ui(denom1, denom1, 1, prec);
	arb_rising_ui(denom1, denom1, n, prec);

	arb_div(res1, numer1, denom1, prec);
	arb_clear(numer1);arb_clear(denom1);

	arb_set(numer2, kon);
	arb_rising_ui(numer2, numer2, n, prec);

	arb_add(denom2, kon, koff,prec);
	arb_rising_ui(denom2, denom2, n, prec);

	arb_div(res2, numer2, denom2, prec);
	arb_clear(numer2);arb_clear(denom2);
}

void get_coef_second(arb_t res, const arb_t kon, const arb_t koff, const arb_t ksyn, const arb_t t, ulong n, slong prec)
{
	arb_t numer, denom;
	arb_init(numer);arb_init(denom);

	arb_t numer_part1, numer_part2;
	arb_init(numer_part1);arb_init(numer_part2);

	arb_set(numer_part1, ksyn);
	arb_pow_ui(numer_part1, numer_part1, n+1,prec);

	arb_set(numer_part2, kon);
	arb_add(numer_part2, numer_part2, koff, prec);
	arb_neg(numer_part2, numer_part2);
	arb_mul(numer_part2, numer_part2, t, prec);
	arb_exp(numer_part2, numer_part2, prec);

	arb_mul(numer, numer_part1, numer_part2, prec);
	arb_mul(numer,numer,kon,prec);

	arb_clear(numer_part1);arb_clear(numer_part2);

	arb_t denom_part1, denom_part2, denom_part3;
	arb_init(denom_part1);arb_init(denom_part2);arb_init(denom_part3);

	arb_set(denom_part1, kon);
	arb_add(denom_part1, denom_part1, koff, prec);

	arb_set(denom_part2, kon);
	arb_add(denom_part2, denom_part2, koff, prec);
	arb_neg(denom_part2, denom_part2);
	arb_add_ui(denom_part2, denom_part2, 1, prec);

	arb_fac_ui(denom_part3, n, prec);

	arb_mul(denom, denom_part1, denom_part2, prec);
	arb_mul(denom, denom, denom_part3, prec);

	arb_clear(denom_part1);arb_clear(denom_part2);arb_clear(denom_part3);

	arb_div(res, numer, denom, prec);
}
void get_factors_second(arb_t res1, arb_t res2, const arb_t kon, const arb_t koff, const arb_t t, ulong n, slong prec)
{
	arb_t numer1, numer2, denom1, denom2;
	arb_init(numer1);arb_init(numer2);arb_init(denom1);arb_init(denom2);

	arb_t numer1_part1, numer1_part2;
	
	arb_set(numer1_part1, koff);
	arb_rising_ui(numer1_part1, numer1_part1, n, prec);

	arb_set(numer1_part2, t);
        arb_mul_ui(numer1_part2, numer1_part2, n, prec);
        arb_neg(numer1_part2, numer1_part2);
        arb_exp(numer1_part2, numer1_part2, prec);

	arb_mul(numer1, numer1_part1, numer1_part2, prec);
        if ( n % 2 == 1)
        {
                arb_neg(numer1, numer1);
        }

	arb_clear(numer1_part1);arb_clear(numer1_part2);

	arb_add(denom1, kon, koff, prec);
	arb_add_ui(denom1, denom1, 1, prec);
	arb_rising_ui(denom1, denom1, n, prec);

	arb_div(res1, numer1, denom1, prec);
        arb_clear(numer1);arb_clear(denom1);

	arb_set(numer2, koff);
	arb_neg(numer2, numer2);
	arb_add_ui(numer2,numer2, 1, prec);
        arb_rising_ui(numer2, numer2, n, prec);

        arb_add(denom2, kon, koff, prec);
	arb_neg(denom2, denom2);
	arb_add_ui(denom2, denom2, 2, prec);
        arb_rising_ui(denom2, denom2, n, prec);

        arb_div(res2, numer2, denom2, prec);
        arb_clear(numer2);arb_clear(denom2);
}
void calculate_probs(double* res_array, double kon_double, double koff_double, double ksyn_double, double t_double, ulong nmax_ulong, ulong prec_ulong)
{
        // Setting up variables
	slong i;

	fmpz_t nmax_fmpz, prec_fmpz;
	fmpz_init_set_ui(nmax_fmpz, nmax_ulong);
	fmpz_init_set_ui(prec_fmpz, prec_ulong);

	slong nmax, prec;

	nmax = fmpz_get_si(nmax_fmpz);
	prec = fmpz_get_si(prec_fmpz);

	// Turn parameters of string type to arb variables
	arb_t kon,koff,ksyn,t;
        arb_init(kon);arb_init(koff);arb_init(ksyn),arb_init(t);
        arb_set_d(kon, kon_double);arb_set_d(koff, koff_double);arb_set_d(ksyn, ksyn_double);arb_set_d(t,t_double);
	
	arb_ptr res_first_1;
        res_first_1 = _arb_vec_init(nmax+1);

        arb_ptr res_first_2;
        res_first_2 = _arb_vec_init(nmax+1);

        arb_ptr res_second_1;
        res_second_1 = _arb_vec_init(nmax+1);

        arb_ptr res_coef_first;
        res_coef_first = _arb_vec_init(nmax+1);

        arb_ptr res_factors_first_r;
        res_factors_first_r = _arb_vec_init(nmax+1);

        arb_ptr res_factors_first_n_minus_r;
        res_factors_first_n_minus_r = _arb_vec_init(nmax+1);

        arb_ptr res_coef_second;
        res_coef_second = _arb_vec_init(nmax+1);

        arb_ptr res_factors_second_r;
        res_factors_second_r = _arb_vec_init(nmax+1);

        arb_ptr res_factors_second_n_minus_r;
        res_factors_second_n_minus_r = _arb_vec_init(nmax+1);

        arb_ptr res_second_2;
        res_second_2 = _arb_vec_init(nmax+2);

        for (i = 0; i <= nmax; i++){
                get_hyp1f1_first_1(&res_first_1[i], kon, koff, ksyn, t, i, prec);
                get_hyp1f1_first_2(&res_first_2[i], kon ,koff, ksyn, i, prec);
                get_hyp1f1_second_1(&res_second_1[i], kon, koff, ksyn, t, i, prec);
                get_coef_first(&res_coef_first[i], ksyn, i, prec);
                get_coef_second(&res_coef_second[i], kon, koff, ksyn, t, i, prec);
                get_factors_first(&res_factors_first_r[i], &res_factors_first_n_minus_r[i], kon, koff, t, i, prec);
                get_factors_second(&res_factors_second_r[i], &res_factors_second_n_minus_r[i], kon, koff, t, i, prec);
                get_hyp1f1_second_2(&res_second_2[i], kon,koff,ksyn,i, prec);
        }
        get_hyp1f1_second_2(&res_second_2[i], kon,koff,ksyn,i, prec);

        arb_ptr res_probs;
        res_probs = _arb_vec_init(nmax+1);

	
	arb_t res_term_1, res_term_2;
	arb_init(res_term_1);arb_init(res_term_2);

	arb_t temp_term_1, temp_term_2;
	arb_init(temp_term_1);arb_init(temp_term_2);
	
	// The formula is essentially based on multiplication and summations of the above results (and partly rely on values from n-1)
	slong n, r;
	for (n = 0; n <= nmax; n++){
		for (r = 0; r <= n; r++){
			arb_zero(temp_term_1);arb_zero(temp_term_2);

			arb_bin_uiui(temp_term_1, n, r, prec);
			arb_mul(temp_term_1, temp_term_1, &res_factors_first_r[r], prec);
			arb_mul(temp_term_1, temp_term_1, &res_factors_first_n_minus_r[n-r],prec);
			arb_mul(temp_term_1, temp_term_1, &res_first_1[r], prec);
			arb_mul(temp_term_1, temp_term_1, &res_first_2[n-r], prec);
			arb_add(res_term_1, res_term_1, temp_term_1, prec);

			arb_bin_uiui(temp_term_2, n, r, prec);
			arb_mul(temp_term_2, temp_term_2, &res_factors_second_r[r], prec);
			arb_mul(temp_term_2, temp_term_2, &res_factors_second_n_minus_r[n-r],prec);
			arb_mul(temp_term_2, temp_term_2, &res_second_1[r], prec);
			arb_mul(temp_term_2, temp_term_2, &res_second_2[n-r+1],prec);
			arb_add(res_term_2, res_term_2, temp_term_2, prec);
		}
		
		arb_mul(res_term_1, res_term_1, &res_coef_first[n], prec);
		arb_mul(res_term_2, res_term_2, &res_coef_second[n], prec);

		arb_add(&res_probs[n], &res_probs[n], res_term_1, prec);
		arb_add(&res_probs[n], &res_probs[n], res_term_2, prec);


		if  (n != nmax){
			arb_sub(&res_probs[n+1], &res_probs[n+1], res_term_2, prec);
		}	

		arb_zero(res_term_1);arb_zero(res_term_2);
	}
	arb_clear(res_term_1);arb_clear(res_term_2);
	arb_clear(temp_term_1);arb_clear(temp_term_2);
	
	
	// Populating the result array with mid point of arb ball

	for (n = 0; n <= nmax; n++){
//		arb_printd(&res_probs[n], 20);printf("\n");
//		printf("%f\n", arf_get_d(arb_midref(&res_probs[n]), MPFR_RNDF));
		res_array[n] = arf_get_d(arb_midref(&res_probs[n]), MPFR_RNDF);
	}

	// Clean up memory

	arb_clear(kon);arb_clear(koff);arb_clear(ksyn),arb_clear(t);

	_arb_vec_clear(res_probs, nmax+1);
	
	_arb_vec_clear(res_first_1,nmax+1);
	_arb_vec_clear(res_first_2,nmax+1);
	_arb_vec_clear(res_coef_first,nmax+1);
	_arb_vec_clear(res_factors_first_r,nmax+1);
	_arb_vec_clear(res_factors_first_n_minus_r,nmax+1);
	
	_arb_vec_clear(res_second_1,nmax+1);
	_arb_vec_clear(res_second_2,nmax+2);
	_arb_vec_clear(res_coef_second,nmax+1);
	_arb_vec_clear(res_factors_second_r,nmax+1);
	_arb_vec_clear(res_factors_second_n_minus_r,nmax+1);
}

void free_array(double* arr){
	free(arr);
	arr = NULL;
}
int main()
{
        ulong prec_ulong = 256;
        const ulong nmax_ulong = 20;

	double* res_array;
        res_array = malloc(sizeof(double)*(nmax_ulong+1));

        calculate_probs(res_array,1.17618015,25.0,60.0,1.0,nmax_ulong,prec_ulong);

        int i;

        for(i = 0; i <= 20; i++){
                printf("%d,%.10f\n",i, res_array[i]);
        }



}


