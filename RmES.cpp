#include <bits/stdc++.h>
#include <iostream>
#include <random>
using namespace std;

#define N 2
#define M 2
#define T 2
#define max_iter 10000

double func(vector<double> x){
    //return (x[0]*x[0] + 100*x[1]*x[1]);
    //return (100*(x[1]-x[0]*x[0])*(x[1]-x[0]*x[0])+(x[0]-1)*(x[0]-1));
    return (x[0]-1)*(x[0]-1)+ 2*(2*x[1]*x[1]-x[0])*(2*x[1]*x[1]-x[0]);
}

bool comp(vector<double> a, vector<double> b){
    return (func(a) <= func(b));
}

void find_rank(vector<vector<double>> F,vector<vector<double>> &R,int t,int mu){
    set<double> Funion;
    for(int i=0;i<mu;i++){
        Funion.insert(F[t][i]);
        Funion.insert(F[t+1][i]);
    }

    for(int i = 0;i<mu;i++){
        for(auto it = Funion.begin();it!=Funion.end();it++){
            if(F[t][i] == *it){
                R[t][i] = distance(Funion.begin(),it);
                break;
            }
        }
        for(auto it = Funion.begin();it!=Funion.end();it++){
            if(F[t+1][i] == *it){
                R[t+1][i] = distance(Funion.begin(),it);
                break;
            }
        }
    }

}

void update(vector<vector<double>> &P,vector<int> &tcap,int T,int m,int iter,vector<vector<double>>p)
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
    if(Tmin>T||t<m)
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
        P[m-1][j]=p[iter][j];
}


vector<double> RmES(){
    
    vector<double> sigma(max_iter+1,0);
    vector<vector<double>> m(max_iter+1,vector<double>(N,0));
    vector<vector<double>> p(max_iter+1,vector<double>(N,0));
    vector<vector<double>> P(M,vector<double>(N,0));
    vector<int> tcap(M,0);
    vector<double> xbest(N,0);
    vector<double> s(max_iter+1,0);

    sigma[0]=20/3.0;
    srand (time(NULL));
    for(int n=0;n<N;n++)
    {
    double fn = ((float)rand())/RAND_MAX;
    float fr = 2*(fn-0.5);
    m[0][n]= 10*fr;
    xbest[n]=m[0][n];
    cout<<xbest[n]<<" ";
    }
    cout<<endl;

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
    double dsigma = 1;
    double cs = 0.3;

    vector<vector<double>> x (lambda,vector<double>(N,0));
    vector<vector<double>> R (max_iter+1,vector<double>(mu,0));
    vector<vector<double>> F (max_iter+1,vector<double>(mu,0));
    double fm = func(m[0]);
    for(int i = 0;i<mu;i++){
        F[0][i] = fm;
    }

    for (int t = 0; t < max_iter; t++){
        for(int i = 0; i<lambda; i++){
            std::random_device rd{};
            std::mt19937 gen{rd()};
            default_random_engine generator;
            normal_distribution<double> distribution(0,1.0);
            //double r = distribution(generator);
            //cout<<"r:"<<r<<endl;
            vector<double> z(N);
            for(int n=0; n <N;n++){
                //z[n] = distribution(generator);
                z[n]=distribution(gen);//cout<<z[n];
                int val = 0;
                for(int j=0;j<M;j++){
                    double r  = distribution(gen);
                    val+=  pow(sqrt(1-ccov),M-(j+1))*r*P[j][n];
                }
                x[i][n] = m[t][n] + sigma[t]*(pow(sqrt(1-ccov),M)*z[n] + sqrt(ccov)*val);
            }
        //cout<<endl;
            if(func(x[i]) < func(xbest)){
                for(int n=0; n<N;n++){
                    xbest[n] = x[i][n];
            //cout<<xbest[n];
                }
        //cout<<endl;
            }

        }
    

    cout<<"fval:"<<func(xbest)<<endl;
    for(int i=0;i<lambda;i++)
    {
        cout<<"("<<x[i][0]<<","<<x[i][1]<<"):"<<func(x[i])<<" ";
    }
    cout<<endl;
        sort(x.begin(),x.end(),comp);
    for(int i=0;i<lambda;i++)
    {
        cout<<"("<<x[i][0]<<","<<x[i][1]<<"):"<<func(x[i])<<" ";
    }
    cout<<endl;
        for(int i = 0;i<mu;i++){
            F[t+1][i] = func(x[i]);
        }

    
        for(int n=0;n<N;n++){
            for(int i=0;i<mu;i++){
                m[t+1][n] += w[i]*x[i][n];
            } 
            p[t+1][n] = (1-cc)*p[t][n] + sqrt(cc*(2-cc)*mueff)*(m[t+1][n] - m[t][n])/dsigma;
        }
    /*cout<<"m=";
    for(int n=0;n<N;n++)
        cout<<m[t+1][n]<<" ";
    cout<<endl<<"p=";
    for(int n=0;n<N;n++)
        cout<<p[t+1][n]<<" ";
    cout<<endl;*/


        update(P,tcap,T,M,t,p);

        find_rank(p,R,t,mu);

        for(int i=0;i<mu;i++){
            q += w[i]*(R[t][i]-R[t+1][i]);
        }
        q /= mu;

        s[t+1] = (1-cs)*s[t] +cs*(q-qstar);
        sigma[t+1] = sigma[t]*exp(s[t+1]/dsigma);cout<<sigma[t+1]<<endl;
    }
    return xbest;
}

int main()
{
    vector<double> v = RmES();
    cout<<v[0]<<" "<<v[1]<<endl;
    return 0;
}

