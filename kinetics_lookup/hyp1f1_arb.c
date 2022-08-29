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
        arb_clear(a);arb_clear(b);arb_clear(z);
}

double* calculate_hyp1f1(char* kon_str, char* koff_str, char* ksyn_str, char* t_str, ulong nmax_ulong, ulong prec_ulong)
{
	const slong size = (nmax_ulong+1)*4 + 1;

	double* res_array;
	res_array = malloc(sizeof(double)*size);
        slong i;

	fmpz_t nmax_fmpz, prec_fmpz;
	fmpz_init_set_ui(nmax_fmpz, nmax_ulong);
	fmpz_init_set_ui(prec_fmpz, prec_ulong);

	slong nmax, prec;

	nmax = fmpz_get_si(nmax_fmpz);
	prec = fmpz_get_si(prec_fmpz);

	arb_t kon,koff,ksyn,t;
        arb_init(kon);arb_init(koff);arb_init(ksyn),arb_init(t);
        arb_set_str(kon, kon_str, prec);arb_set_str(koff, koff_str, prec);arb_set_str(ksyn, ksyn_str, prec);arb_set_str(t,t_str,prec);

	arb_t res;
	arb_init(res);
	for (i = 0; i <= nmax; i++){
		get_hyp1f1_first_1(res, kon,koff,ksyn,t,i, prec);
		res_array[i] = arf_get_d(arb_midref(res),MPFR_RNDF);
		arb_zero(res);
	}
	slong offset;
	slong idx;
	offset = nmax+1;
	for (i = 0; i <= nmax; i++){
                get_hyp1f1_first_2(res, kon,koff,ksyn,i, prec);
		idx = i + offset;
                res_array[idx] = arf_get_d(arb_midref(res),MPFR_RNDF);
                arb_zero(res);
        }
	offset = 2*(nmax+1);
	for (i = 0; i <= nmax; i++){
                get_hyp1f1_second_1(res, kon,koff,ksyn,t,i, prec);
		idx = i + offset;
		res_array[idx] = arf_get_d(arb_midref(res),MPFR_RNDF);
                arb_zero(res);
        }
	offset = 3*(nmax+1);
	for (i = 0; i <= nmax; i++){
                get_hyp1f1_second_2(res, kon,koff,ksyn,i, prec);
		idx = i + offset;
		res_array[idx] = arf_get_d(arb_midref(res),MPFR_RNDF);
                arb_zero(res);
        }

	get_hyp1f1_second_2(res, kon,koff,ksyn,i, prec);
	idx += 1;
	res_array[idx] = arf_get_d(arb_midref(res),MPFR_RNDF);
	arb_clear(kon);arb_clear(koff);arb_clear(ksyn),arb_clear(t);arb_clear(res);	
	return res_array;
}




int main()
{
	ulong prec_ulong = 256;	
	const ulong nmax_ulong = 10;
	double* str;

	str = calculate_hyp1f1("1.17618015","25","60","1",nmax_ulong,prec_ulong);
	
	int i;
	
	for(i = 0; i < 45; i++){
		printf("%d,%.10f\n",i, str[i]);
	}



}

