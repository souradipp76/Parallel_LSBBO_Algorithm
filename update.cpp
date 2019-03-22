#include <bists/stdc++.h>

using namespace std;
#define N 2

int max_iter=10000;

void update(vector<vector<double>> &P,vector<int> &tcap,int T,int m,int iter,vector<vector<doubl>>p)
{
	int Tmin=INT_MAX,index;
	for(int i=0;i<m-1;i++)
	{
		if(tcap[i+1]-tcap[i])
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
	tcap[m-1]=t;
	for(int j=0;j<N;j++)
		P[m-1][j]=p[iter][j];
}
