#include <iostream>
#include <cmath>
#include <vector>
using namespace std;
/************************************************
多项式相关计算，用于生成基函数
************************************************/
struct poly
{
    vector<double> a;
    poly(int n):a(n+1,0){}
    poly(){}
    void set(int n)
    {
        a.resize(n);
        for(int i=0;i<n;i++) a[i]=0;
    }
    const int size()
    {
        return a.size();
    }
    const int size()const
    {
        return a.size();
    }
    double& operator[](int i)
    {
        return a[i];
    }
    const double& operator[](int i)const
    {
        return a[i];
    }
    double operator()(double x)
    {
        double r=0;
        for(int i=a.size()-1;i>0;i--)
        {
            r+=x*a[i];
        }
        return r+a[0];
    }
    poly integ()//不定积分
    {
        poly t;
        t.set(a.size()+1);
        for(int i=1;i<t.size();i++) t[i]=a[i-1]/i;
        return t;
    }
    double integ(double l,double u)//定积分
    {
        poly t;
        t.set(a.size()+1);
        for(int i=1;i<t.size();i++) t[i]=a[i-1]/i;
        return t(u)-t(l);
    }
    poly diff()//导函数
    {
        poly t;
        t.set(a.size()-1);
        for(int i=0;i<t.size();i++) t[i]=a[i+1]*(i+1);
        return t;
    }
    poly diff() const
    {
        poly t;
        t.set(a.size()-1);
        for(int i=0;i<t.size();i++) t[i]=a[i+1]*(i+1);
        return t;
    }
    double diff(double x)//x处导数值
    {
        poly t;
        t.set(a.size()-1);
        for(int i=0;i<t.size();i++) t[i]=a[i+1]*i;
        return t(x);
    }
};
poly operator*(const double k,const poly p)
{
    poly t=p;
    for(int i=0;i<p.size();i++) t[i]=k*p[i];
    return t;
}
poly operator*(const poly p,const double k)
{
    poly t=p;
    for(int i=0;i<p.size();i++) t[i]=k*p[i];
    return t;
}
poly operator/(const double k,const poly p)
{
    poly t=p;
    for(int i=0;i<p.size();i++) t[i]=k/p[i];
    return t;
}
poly operator/(const poly p,const double k)
{
    poly t=p;
    for(int i=0;i<p.size();i++) t[i]=k/p[i];
    return t;
}
poly operator*(const poly p1,const poly p2)
{
    poly t;
    t.set(p1.size()+p2.size()-1);
    for(int i=0;i<p1.size();i++)
    {
        for(int j=0;j<p2.size();j++) t[i+j]+=p1[i]*p2[j];
    }
    int i=t.size(),j=i-1;
    for(;j>=0;j--)
    {
        if(t[j]!=0) break;
    }
    t.a.resize(j+1);
    return t;
}
ostream& operator<<(ostream& out,const poly p)
{
    for(int i=p.size()-1;i>0;i--) out<<p[i]<<'\t';
    out<<p[0];
    return out;
}
vector<poly> creat_base(int n)//生成拉格朗日插值基函数
{
    vector<poly> r(n+1,poly(0));
    for(int i=0;i<n+1;i++)
    {
        r[i][0]=1;
        for(int j=0;j<n+1;j++)
        {
            if(j==i) continue;
            poly t(1);
            t[0]=-j/(i-j);
            t[1]=n/(i-j);
            r[i]=r[i]*t;
        }
    }
    return r;
}
/************************************************
矩阵相关计算
************************************************/
struct matrix
{
    vector<vector<double>> a;
    matrix(){}
    matrix(int r,int c):a(r,vector<double>(c,0)){}
    void set(int r,int c)
    {
        a.resize(r);
        for(int i=0;i<r;i++) 
        {
            a[i].resize(c);
            for(int j=0;j<c;j++) a[i][j]=0;
        }
    }
    vector<double>& operator[](int i)
    {
        return a[i];
    }
    const vector<double>& operator[](int i)const
    {
        return a[i];
    }
    const int row()
    {
        return a.size();
    }
    const int col()
    {
        if(a.size()==0) return 0;
        return a[0].size();
    }
    const int row() const
    {
        return a.size();
    }
    const int col() const
    {
        if(a.size()==0) return 0;
        return a[0].size();
    }
    matrix operator-()
    {
        matrix A(row(),col());
        for(int i=0;i<row();i++)
        {
            for(int j=0;j<col();j++)
            {
                A[i][j]=-a[i][j];
            }
        }
        return A;
    }
};
matrix operator+(const matrix& A,const matrix& B)
{
    if(A.row()!=B.row()||A.col()!=B.col())
    {
        cout<<"matrix operator+:维度不匹配"<<endl;
        return matrix();
    }
    matrix C(A.row(),A.col());
    for(int i=0;i<C.row();i++)
    {
        for(int j=0;j<C.col();j++)
        {
            C[i][j]=A[i][j]+B[i][j];
        }
    }
    return C;
}
matrix operator-(const matrix& A,const matrix& B)
{
    if(A.row()!=B.row()||A.col()!=B.col())
    {
        cout<<"matrix operator+:维度不匹配"<<endl;
        return matrix();
    }
    matrix C(A.row(),A.col());
    for(int i=0;i<C.row();i++)
    {
        for(int j=0;j<C.col();j++)
        {
            C[i][j]=A[i][j]-B[i][j];
        }
    }
    return C;
}
matrix operator*(const double k,const matrix& A)
{
    matrix C(A.row(),A.col());
    for(int i=0;i<C.row();i++)
    {
        for(int j=0;j<C.col();j++)
        {
            C[i][j]=A[i][j]*k;
        }
    }
    return C;
}
vector<double> operator*(const matrix& A,const vector<double>& x)
{
    if(A.row()!=x.size())
    {
        cout<<"operator *:维度不匹配"<<endl;
        return vector<double>();
    }
    vector<double> r(x.size(),0);
    for(int i=0;i<A.row();i++)
    {
        for(int j=0;j<A.col();j++)
        {
            r[i]+=A[i][j]*x[j];
        }
    }
    return r;
}
matrix operator*(const matrix& A,const double k)
{
    matrix C(A.row(),A.col());
    for(int i=0;i<C.row();i++)
    {
        for(int j=0;j<C.col();j++)
        {
            C[i][j]=A[i][j]*k;
        }
    }
    return C;
}
matrix operator/(const matrix& A,const double k)
{
    matrix C(A.row(),A.col());
    for(int i=0;i<C.row();i++)
    {
        for(int j=0;j<C.col();j++)
        {
            C[i][j]=A[i][j]/k;
        }
    }
    return C;
}
ostream& operator<<(ostream& out,const matrix& A)
{
    for(int i=0;i<A.row();i++)
    {
        for(int j=0;j<A.col()-1;j++)
        {
            out<<A[i][j]<<'\t';
        }
        out<<A[i][A.col()-1]<<'\n';
    }
    return out;
}
/************************************************
矢量和矩阵相关计算
************************************************/
ostream& operator<<(ostream& out,const vector<double>& A)
{
        for(int j=0;j<A.size()-1;j++)
        {
            out<<A[j]<<'\t';
        }
        out<<A[A.size()-1];
    return out;
}
double norm(const vector<double>& x)//二范数、模
{
    double r=0;
    for(int i=0;i<x.size();i++) r+=x[i]*x[i];
    return sqrt(r);
}
vector<double>& operator/(vector<double>& A,const double k)
{
    for(int i=0;i<A.size();i++)
    {
        A[i]=A[i]/k;
    }
    return A;
}
double operator*(const vector<double>& a,const vector<double>& b)
{
    if(a.size()!=b.size())
    {
        cout<<"点积维度不匹配"<<endl;
        return 0;
    }
    double r=0;
    for(int i=0;i<a.size();i++) r+=a[i]*b[i];
    return r;
}
double power_method(const matrix& A)//幂法求谱半径
{
    vector<double> t(A.row(),1),t0(t);
    double beta=0,beta0=1;
    while(abs(beta-beta0)>1e-7)
    {
        beta0=beta;
        t=t/norm(t);
        t0=t;
        t=A*t;
        beta=t0*t;
    }
    return beta;
}
void op1(matrix& A,int i,int j)
{
    for(int k=0;k<A[0].size();k++)
    {
        swap(A[i][k],A[j][k]);
    }
}
void op2(matrix& A,int i,double k)
{
    for(int j=0;j<A[0].size();j++)
    {
        A[i][j]*=k;
    }
}
void op3(matrix& A,int i,int j,double k)
{
    for(int l=0;l<A[0].size();l++)
    {
        A[j][l]+=k*A[i][l];
    }
}
int maxr(const matrix &A,int j)
{
    int Mr=j;
    double M=A[j][j];
    for(int i=j+1;i<A.a.size();i++)
    {
        if(abs(A[i][j])>M)
        {
            M=abs(A[i][j]);
            Mr=i;
        }
    }
    return Mr;
}
void inv_mul(matrix A,matrix& B)//inv(A)*B
{
    int rank=A.a.size();
    for(int i=0;i<A.a.size();i++)
    {
        int m=maxr(A,i);
        op1(A,i,m);
        op1(B,i,m);
        op2(B,i,1/A[i][i]);
        op2(A,i,1/A[i][i]);
        for(int j=i+1;j<A.a.size();j++)
        {
            op3(B,i,j,-A[j][i]);
            op3(A,i,j,-A[j][i]);
        }
    }
    for(int i=A.a.size()-1;i>=0;i--)
    {
        for(int k=i-1;k>=0;k--)
        {
            op3(B,i,k,-A[k][i]);
        } 
    }
}

matrix creat_At(int n,const vector<poly>& N)//获得时间项在标准单元上的系数矩阵
{
    if(N.size()!=n+1)
    {
        cout<<"creat_At:阶数不匹配"<<endl;
        return matrix();
    }
    matrix _At(n+1,n+1);
    for(int i=0;i<n+1;i++)
    {
        for(int j=0;j<n+1;j++)
        {
            _At[i][j]=(N[i]*N[j]).integ(0,1);
        }
    }
    return _At;
}
matrix creat_Ax(int n,const vector<poly>& N)//获得一阶空间导数项在标准单元上的系数矩阵
{
    if(N.size()!=n+1)
    {
        cout<<"creat_Ax:阶数不匹配"<<endl;
        return matrix();
    }
    matrix _Ax(n+1,n+1);
    for(int i=0;i<n+1;i++)
    {
        for(int j=0;j<n+1;j++)
        {
            poly t=N[j].diff();
            _Ax[i][j]=-(N[i]*t).integ(0,1);
        }
    }
    return _Ax;
}
matrix creat_Axx(int n,const vector<poly>& N)//获得二阶空间导数项在标准单元上的系数矩阵
{
    if(N.size()!=n+1)
    {
        cout<<"creat_Ax:阶数不匹配"<<endl;
        return matrix();
    }
    matrix _Ax(n+1,n+1);
    for(int i=0;i<n+1;i++)
    {
        for(int j=0;j<n+1;j++)
        {
            poly t1=N[i].diff(),t2=N[j].diff();
            _Ax[i][j]=(t1*t2).integ(0,1);
        }
    }
    return _Ax;
}
matrix get_A(int n,int o,const matrix& _A)//合并所有单元(该函数要求所有单元长度相等，如果不等需要额外修改合并方式)
{
    matrix A(o*n+1,o*n+1);
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<o+1;j++)
        {
            for(int k=0;k<o+1;k++)
            {
                A[i*o+j][i*o+k]+=_A[j][k];
            }
        }
    }
    return A;
}
const int order=1;
const int NE=32,//空间点数(阶数)
NS=100;//时间步数
const double rb=-5,l=10,//计算域左边界，计算域长度
dt=0.1,//时间步长
dx=l/NE;
void init_guass(vector<double> &u0)//设置高斯函数初场
{
    for(int i=0;i<u0.size();i++) 
	{
		u0[i]=exp(-pow((l*double(i)/(NE-1)+rb),2));
	}
}
int main()
{
    vector<double> u0(NE+1);
    init_guass(u0);
    vector<poly> N=creat_base(order);
    matrix _At=creat_At(order,N),_Ax=creat_Ax(order,N),_Axx=creat_Axx(order,N),
    At=get_A(NE,order,_At),Ax=get_A(NE,order,_Ax),Axx=get_A(NE,order,_Axx),A(NE+1,NE+1);
    double nu=0.5;
    A=At+Ax*(dt/dx)-nu*Axx*(dt/dx*dx);
    inv_mul(At,A);
    A[0][A.col()-1]=1;
    for(int i=0;i<A.col()-1;i++) A[0][i]=0;
    cout<<NE+1<<'\t'<<NS<<'\t'<<rb<<'\t'<<l<<'\n';
    cout<<u0<<'\n';
    for(int i=0;i<NS;i++)
    {
        u0=A*u0;
        cout<<u0<<'\n';
    }
    return 0;
}