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
			for(int j = 2; i*j < 100001; ++j){
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


// Miller rabin 


bool miller(ll n){
	if(n <=1 || n % 2 == 0){
		if(n != 2){
			return false;
		}
	}
	if(n == 2 || n == 3){
		return true;
	}
	ll d = n-1;
	while(d % 2 == 0){
		d /= 2;
	}
	ll a[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
	for(int i = 0; i < 12 && a[i] < n; ++i){
		ll temp = d;
		ll mod = binexp2(a[i], temp, n);
		if( mod == 1 || mod == n-1){
			continue;
		}
		while(temp != n-1 && mod != n-1){
			mod = multiply(mod, mod, n);
			temp *= 2;
		}
		if(mod != n-1){
			return false;
		}  
	}
	return true;
}


// Exteneded Euclid 


struct node{
	ll x, y;
};
struct node exteuclid(ll a, ll b){
	if(a % b == 0){
		node ret;
		ret.x = 0;
		ret.y = 1;
		return ret;
	}
	node ans = exteuclid(b, a % b);
	node ret;
	ret.x = ans.y;
	ret.y = ans.x - (a / b) * ans.y;
	return ret;
}
