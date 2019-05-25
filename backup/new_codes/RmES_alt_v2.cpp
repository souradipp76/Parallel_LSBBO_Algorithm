#include <bits/stdc++.h>
#include <iostream>
#include <random>
using namespace std;

#define N 200
#define M 2
#define max_iter 100000

double func(vector<double> &x,int fnum){
    double res=0;
switch(fnum){
case 1:    
    ///fEll(x)////
    for(int n = 0; n < N; n++)
    {
        res+=(pow(10,6.0*(n/float(N)))*x[n]*x[n]);
    }
    ////
break;
case 2:
    // ////fCig(x)/////
    res+=x[0]*x[0];
    for(int n = 0; n < N-1; n++)
     {
         res+=1000000*(x[n+1]*x[n+1]);
     }
    // /////
break;

case 3:
     /////fCtb(x)////
     res+=x[0]*x[0];
     for(int n = 0; n < N-2; n++)
     {
         res+=10000*(x[n+1]*x[n+1]);
     }
     res+=1000000*x[N-1]*x[N-1];
    // /////
break;
case 4:

    // /////fTab(x)////
     res+=1000000*x[0]*x[0];
     for(int n = 0; n < N-1; n++)
     {
         res+=(x[n+1]*x[n+1]);
     }
    // /////
break;

case 5:
    // /////fSch(x)////
     for(int n = 0; n < N; n++)
     {
      double temp=0;
      for(int j=0;j<n;j++)
          temp+=x[n]*x[n];
         res+=temp*temp;
     }
     ////
break;
case 6:

    /////fRos(x)////
     for(int n = 0; n < N-1; n++)
     {
         res+=(100*(x[n+1]-x[n]*x[n])*(x[n+1]-x[n]*x[n])+(x[n]-1)*(x[n]-1));
     }
    ////
break;
}
    return res;
}
int partition(vector<vector<double>> &a,int l,int u,int num)
{
    int i = (l - 1); 
    double fpivot = func(a[u],num);  

    for (int j = l; j <= u-1; j++) 
    { 
        if (func(a[j],num) <= fpivot) 
        { 
            i++; 
            a[i].swap(a[j]); 
        } 
    } 
    a[i + 1].swap(a[u]); 
    return (i + 1); 
}

void quick_sort(vector<vector<double>> &a,int l,int u,int num)
{
    int j;
    if(l<u)
    {
        j=partition(a,l,u,num);
        quick_sort(a,l,j-1,num);
        quick_sort(a,j+1,u,num);
    }
}


void find_rank(vector<vector<double>> F,vector<vector<double>> &R,int t,int mu){
    set<double> Funion;
    for(int i=0;i<mu;i++){
        Funion.insert(F[t&1][i]);
        Funion.insert(F[(t+1)&1][i]);
    }

    for(int i = 0;i<mu;i++){
        for(auto it = Funion.begin();it!=Funion.end();it++){
            if(F[t&1][i] == *it){
                R[t&1][i] = distance(Funion.begin(),it);
                break;
            }
        }
        for(auto it = Funion.begin();it!=Funion.end();it++){
            if(F[(t+1)&1][i] == *it){
                R[(t+1)&1][i] = distance(Funion.begin(),it);
                break;
            }
        }
    }

}

void update(vector<vector<double>> &P,vector<int> &tcap,int T,int m,int iter,vector<vector<double>>&p)
{
    int Tmin=INT_MAX,index;
    for(int i=0;i<m-1;i++)
    {
        if(tcap[i+1]-tcap[i]<Tmin)
        {
            Tmin=tcap[i+1]-tcap[i];
            index=i;
        }
        
    }
    if(Tmin>T||iter<m)
    {
        for(int i=0;i<m-1;i++)
        {
            for(int j=0;j<N;j++)
                P[i][j]=P[i+1][j];
            tcap[i]=tcap[i+1];
        }
    }
    else
    {
        for(int i=index;i<m-1;i++)
        {
            for(int j=0;j<N;j++)
                P[i][j]=P[i+1][j];
            tcap[i]=tcap[i+1];
        }   
    }
    tcap[m-1]=iter;
    for(int j=0;j<N;j++)
        P[m-1][j]=p[iter&1][j];
}


vector<double> RmES(int num){
    
    double sigma;
    vector<vector<double>> m(2,vector<double>(N,0));
    vector<vector<double>> p(2,vector<double>(N,0));
    vector<vector<double>> P(M,vector<double>(N,0));
    vector<int> tcap(M,0);
    vector<double> xbest(N,0);
    double s;

    sigma=20/3.0;
    srand (time(NULL));
    for(int n=0;n<N;n++)
    {
        double fn = ((float)rand())/RAND_MAX;
        float fr = 2*(fn-0.5);
        m[0][n]= 10*fr;
        xbest[n]=m[0][n];
        //cout<<xbest[n]<<" ";
    }
    //cout<<endl;

    int lambda = 4 + floor(3*log(N));//cout<<log(N)<<endl;
    int mu = floor(lambda/2);
    vector<double> w(mu,0);
    double mueffdem = 0;
    for (int i = 1; i <= mu; i++){
        double wnum = log(mu+1) - log(i);
        double wdem = mu*log(mu+1);
        for(int j = 1; j <= mu; j++){
            wdem = wdem - log(j);
        }
        w[i-1] = wnum/wdem;
        mueffdem = mueffdem + w[i-1]*w[i-1];
    } 
    double mueff = 1.0/mueffdem;
    double ccov =  1.0/(3*sqrt(N) + 5);
    double cc = 2.0/(N+7);
    double qstar = 0.3;
    double q = 0;
    double dsigma = 1000;
    double cs = 0.3;

    vector<vector<double>> x (lambda,vector<double>(N,0));
    vector<vector<double>> R (2,vector<double>(mu,0));
    vector<vector<double>> F (2,vector<double>(mu,0));
    double fm = func(m[0],num);
    for(int i = 0;i<mu;i++){
        F[0][i] = fm;
    }

    int t=0;
    vector<double> z(N);
    double fbest = func(xbest,num);
    while(t < max_iter && fbest>0.001){
        for(int i = 0; i<lambda; i++){
            std::random_device rd{};
            std::mt19937 gen{rd()};
            default_random_engine generator;
            normal_distribution<double> distribution(0,1.0);
            //double r = distribution(generator);
            //cout<<"r:"<<r<<endl;
            for(int n=0; n <N;n++){
                //z[n] = distribution(generator);
                z[n]=distribution(gen);//cout<<z[n];
                int val = 0;
                for(int j=0;j<M;j++){
                    double r  = distribution(gen);
                    val+=  pow(sqrt(1-ccov),M-(j+1))*r*P[j][n];
                }
                x[i][n] = m[t&1][n] + sigma*(pow(sqrt(1-ccov),M)*z[n] + sqrt(ccov)*val);
            }
            //cout<<endl;
            if(func(x[i],num) < fbest){
                for(int n=0; n<N;n++){
                    xbest[n] = x[i][n];
                    //cout<<xbest[n];
                }
                fbest = func(xbest,num);
                //cout<<endl;
            }

        }
    
        //cout<<"Iteration "<<t;//<<": x =(";
        //for(int n=0;n<N;n++)
        //    cout<<xbest[n]<<" ";
        //cout<<")"<<endl;
        //cout<<"fval:"<<fbest<<endl;
        //if(t%250==0)
        //    cout<<t<<","<<fbest<<endl;

        quick_sort(x,0,lambda-1,num);

        for(int i = 0;i<mu;i++){
            F[(t+1)&1][i] = func(x[i],num);
        }

    
        for(int n=0;n<N;n++){
            m[(t+1)&1][n]=0;
            for(int i=0;i<mu;i++){
                m[(t+1)&1][n] += w[i]*x[i][n];
            } 
            p[(t+1)&1][n] = (1-cc)*p[t&1][n] + sqrt(cc*(2-cc)*mueff)*(m[(t+1)&1][n] - m[t&1][n])/sigma;
        }

	    int T=N;
        update(P,tcap,T,M,t,p);

        find_rank(p,R,t,mu);
        q=0;
        for(int i=0;i<mu;i++){
            q += w[i]*(R[t&1][i]-R[(t+1)&1][i]);
        }
        q /= mu;

        s = (1-cs)*s +cs*(q-qstar);
        sigma = sigma*exp(s/dsigma);//cout<<sigma[t+1]<<endl;
        t++;
    }
    return xbest;
}

int main()
{
    int i;
    double T[2];
    for(i=1;i<=6;i++){
        for(int j = 0;j<2;j++){
        clock_t t1=clock(),t2;
        vector<double> v = RmES(i);
        t2=clock();
        T[j]=(double)(t2-t1)/CLOCKS_PER_SEC;
        cout<<i<<","<<j<<","<<T[j]<<endl;
        }
        cout<<"Avg="<<(T[0]+T[1])/2.0<<endl;
        //cout<<"Xopt = (";
        //for(int n=0;n<N;n++)
        //    cout<<v[n]<<", ";
        //cout<<")"<<endl;
        //cout<<"Function Value: "<<func(v,i)<<endl;
    }
    return 0;
}

