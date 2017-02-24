#include <bits/stdc++.h>
#define ll long long
using namespace std;

bool cmp(node a, node b){
	if(a.x != b.x){
		return (a.x < b.x); 
	}
	else{
		return (a.y > b.y);
	}
}

ll divide(ll a[], ll x, ll n){
	ll temp = 0;
	for(int i = 0; i < n; ++i){
		temp = temp *10 + a[i];
		if(temp >= x){
			temp %= x;
		}
	}
	return temp;
}

ll multiply(ll a, ll b, ll m){
	ll result = 0;
	while(b>0){
		if(b % 2 == 1){
			result = (result % m) + (a % m);
			result %= m;
		}
		a = (a % m) + (a % m);
		a %= m;
		b=b/2;
	}
	return result;
}

ll binexp1(ll a, ll b, ll m){
	ll result = 1;
	while(b>0){
		if(b % 2 == 1){
			result = (result % m) * (a % m);
			result %= m;
		}
		a = (a % m) * (a % m);
		a %= m;
		b = b/2;
	}
	return result;
}

ll binexp2(ll a, ll b, ll m){
	ll result = 1;
	while(b>0){
		if(b % 2 == 1){
			result = multiply(result, a, m);
			result %= m;
		}
		a = multiply(a, a, m);
		a %= m;
		b = b/2;
	}
	return result;
}

void sieve(){
	for(int i = 2; i < 100001; ++i){
		if(isprime[i] == 0){
			for(int j = 2; j < 100001; ++j){
				isprime[i*j] = 1;
			}
		}
	}
}

void etf(){
	for(int i=0; i < 100001; ++i)
		phi[i] = i;
	for(int i=2; i < 100001; ++i){
		if(phi[i] == i ){
			for(int j=1; i*j < max; ++j){
				phi[i*j] /= i;
				phi[i*j] *= (i-1);
			}
		}
	}
}