#include <iostream>
#include <cmath>
#include <complex>
#include "fftw3.h"
using namespace std;
const int N=8,//空间点数(阶数)
NS=100;//时间步数
const double pi=atan(1)*4;
const double rb=-5,l=10,//计算域左边界，计算域长度
dt=0.1//时间步长
;
void init_guass(double* u0)//设置高斯函数初场
{
    for(int i=0;i<N;i++) 
	{
		u0[i]=exp(-pow((l*double(i)/(N-1)+rb),2));
	}
}
struct _tran//傅里叶变换的辅助类
{
	fftw_plan s2p,p2s;
	int n;
	_tran(int _n):n(_n)
	{
		double *up;
		fftw_complex *us;
		us  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*n);
		up = (double*) fftw_malloc(sizeof(double)*n);
    	s2p = fftw_plan_dft_c2r_1d(n,us,up,FFTW_MEASURE);
	    p2s = fftw_plan_dft_r2c_1d(n,up,us,FFTW_MEASURE);
		fftw_free(us);
    	fftw_free(up);
	}
	void tran(double* in,fftw_complex* out)//正变换：物理空间到谱空间
	{
		fftw_execute_dft_r2c(p2s,in,out);
	}
	void itran(fftw_complex* in,double* out)//逆变换：谱空间到物理空间
	{
		fftw_execute_dft_c2r(s2p,in,out);
		for(int i=0;i<n;i++) out[i]/=n;
	}
	~_tran()
	{
		fftw_destroy_plan(p2s);
		fftw_destroy_plan(s2p);
	}
}fft(N),fft2(2*N);
void prt(fftw_complex* in)//输出谱系数
{
	for(int i=0;i<N/2+1;i++) cout<<*((complex<double>*)(&in[i]))<<'\t';
    cout<<endl;
}
void prt(double* in)//输出网格点上的值
{
	for(int i=0;i<N;i++) cout<<in[i]<<'\t';
    cout<<endl;
}

void advance(fftw_complex *us)//单步时间推进
{
	complex<double> *a=(complex<double>*)(us);
	for(int i=0;i<N/2+1;i++)
	{
		a[i]=(complex<double>(1,0)-complex<double>(i*dt*(2*pi/l),0)*complex<double>(0,1))*a[i];
	}
}

/*
void advance(fftw_complex *us)//单步时间推进
{
	complex<double> *a=(complex<double>*)(us);
	for(int i=0;i<N/2+1;i++)
	{
        complex<double> t=(complex<double>(1,0)-complex<double>(i*dt*(2*pi/l),0)*complex<double>(0,1));
        t=t/complex<double>(abs(t),0);
		a[i]=t*a[i];
	}
}
*/
/*
void advance(fftw_complex *us)//单步时间推进
{
	complex<double> *a=(complex<double>*)(us);
    fftw_complex *ts=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(N/2+1));
    double *tp=(double*)fftw_malloc(sizeof(double)*N),
    *tp2=(double*)fftw_malloc(sizeof(double)*N);
    for(int i=0;i<N/2+1;i++) 
    {
        *(complex<double>*)(&ts[i])=complex<double>(i*dt*(2*pi/l),0)*complex<double>(0,1)*a[i];
    }
    fft.itran(ts,tp);
    fft.itran(us,tp2);
    for(int i=0;i<N;i++) tp[i]*=tp2[i];
    fft.tran(tp,ts);
    for(int i=0;i<N/2+1;i++)
	{
		a[i]-=(*(complex<double>*)(&ts[i]));
	}
    fftw_free(ts);
    fftw_free(tp);
    fftw_free(tp2);
}
*/
/*
void advance(fftw_complex *us)//单步时间推进
{
	complex<double> *a=(complex<double>*)(us);
    fftw_complex *ts=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(N/2+1));
    double *tp=(double*)fftw_malloc(sizeof(double)*N),
    *tp2=(double*)fftw_malloc(sizeof(double)*N);
    for(int i=0;i<N/2+1;i++) 
    {
        *(complex<double>*)(&ts[i])=complex<double>(i*dt*(2*pi/l),0)*complex<double>(0,1)*a[i];
    }
    fft.itran(ts,tp);
    fft.itran(us,tp2);
    for(int i=0;i<N;i++) tp[i]*=tp2[i];
    fft.tran(tp,ts);
    double nu=0.03;
    for(int i=0;i<N/2+1;i++)
	{
		a[i]-=(*(complex<double>*)(&ts[i]))+nu*complex<double>(pow(i*(2*pi/l),2)*dt,0)*a[i];
	}
    fftw_free(ts);
    fftw_free(tp);
    fftw_free(tp2);
}
*/
/*
void advance(fftw_complex *us)//单步时间推进
{
	complex<double> *a=(complex<double>*)(us);
    fftw_complex *ts=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(N+1)),
	*ts2=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(N+1));
    double *tp=(double*)fftw_malloc(sizeof(double)*N*2),
    *tp2=(double*)fftw_malloc(sizeof(double)*N*2);
    for(int i=0;i<N+1;i++) 
    {
        if(i<N/2+1) 
		{
			*(complex<double>*)(&ts[i])=complex<double>(i*(2*pi/l),0)*complex<double>(0,1)*a[i];
			*(complex<double>*)(&ts2[i])=a[i];
		}
		else
		{
			*(complex<double>*)(&ts[i])=complex<double>(0,0);
			*(complex<double>*)(&ts2[i])=complex<double>(0,0);
		}
	}
    fft2.itran(ts,tp);
    fft2.itran(ts2,tp2);
    for(int i=0;i<2*N;i++) tp[i]*=tp2[i];
    fft2.tran(tp,ts);
	double nu=0.02;
    for(int i=0;i<N/2+1;i++)
	{
		a[i]-=(*(complex<double>*)(&ts[i]))*dt+nu*complex<double>(pow(i*(2*pi/l),2)*dt,0)*a[i];
	}
    fftw_free(ts);
    fftw_free(tp);
    fftw_free(tp2);
}
*/
int main()
{
    fftw_complex *us;
    double *up;
	us  = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(N/2+1));
	up = (double*) fftw_malloc(sizeof(double)*N);
    init_guass(up);
	cout<<N<<'\t'<<NS<<'\t'<<rb<<'\t'<<l<<'\n';
	prt(up);
    fft.tran(up,us);
	//fft.itran(us,up);
	//prt(up);
	for(int i=0;i<NS;i++)
	{
		advance(us);
		fft.itran(us,up);
		prt(up);
	}
	fftw_free(us);
    fftw_free(up);
    return 0;
}
