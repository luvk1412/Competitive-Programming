#include <bits/stdc++.h>
#define ll long long
#define pb push_back
#define mp make_pair
#define fi first
#define se second
#define p(x) printf("%d\n", x)
#define pl(x) printf("%lld\n", x)
#define p2(x, y) printf("%d %d\n", x, y)
#define pl2(x, y) printf("%lld %lld\n", x, y)
#define pf(x) printf("%lf\n", x)
#define s(x) scanf("%d", &x)
#define sl(x) scanf("%lld", &x)
#define sf(x) scanf("%lf", &x)
#define INF 1e18+9 
using namespace std;

void fs(ll &x)
{
    register ll c = getchar_unlocked();
    x = 0;
    for(; (c<48 || c>57); c = getchar_unlocked());
    for(; c>47 && c<58; c = getchar_unlocked()){
    	x = (x<<1) + (x<<3) + c - 48;
    }
}

bool cmp(node a, node b){
	if(a.x != b.x){
		return (a.x < b.x); 
	}
	else{
		return (a.y > b.y);
	}
}

ll divide(string s, ll x){
	vector <ll> a;
	for(int i = 0; i < s.size(); ++i){
		a.pb((ll)s[i] - '0');
	}
	ll temp = 0;
	for(int i = 0; i < a.size(); ++i){
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
		if(b & 1){
			result = result + a;
			result %= m;
		}
		a = a << 1;
		a %= m;
		b = b >> 1;
	}
	return result;
}

ll binexp1(ll a, ll b, ll m){
	ll result = 1;
	while(b>0){
		if(b & 1){
			result = result * a;
			result %= m;
		}
		a = a * a;
		a %= m;
		b = b >> 1;
	}
	return result;
}

ll binexp2(ll a, ll b, ll m){
	ll result = 1;
	while(b>0){
		if(b & 1){
			result = multiply(result, a, m);
			result %= m;
		}
		a = multiply(a, a, m);
		a %= m;
		b = b >> 1;
	}
	return result;
}

void sieve(){
	for(int i = 2; i < 10000001; ++i){
		if(isprime[i] == 0){
			for(int j = 2; i*j < 10000001; ++j){
				isprime[i*j] = 1;
			}
		}
	}
}

void etf(){
	for(int i=0; i < 10000001; ++i)
		phi[i] = i;
	for(int i=2; i < 10000001; ++i){
		if(phi[i] == i ){
			for(int j=1; i*j < 10000001; ++j){
				phi[i*j] /= i;
				phi[i*j] *= (i-1);
			}
		}
	}
}


// miller rabin

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
        if(mod == 1){
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

//   extendid euclid

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

// binary indexed tree

void update(ll i, ll x, ll bit[], ll n){
	while(i <= n){
		bit[i] += x;
		i = i + (i & (-i));
	}
}
ll sum(ll i, ll bit[]){
	ll ans = 0;
	while(i != 0 && i > 0){
		ans += bit[i];
		i = (i & (i-1));
	}
	return ans;
}

Nth-fibonaaci-in-log-n

#define long long long
const long M = 1000000007; // modulo
map<long, long> F;
long f(long n) {
	if (F.count(n)) return F[n];
	long k=n/2;
	if (n%2==0) { // n=2*k
		return F[n] = (f(k)*f(k) + f(k-1)*f(k-1)) % M;
	} else { // n=2*k+1
		return F[n] = (f(k)*f(k+1) + f(k-1)*f(k)) % M;
	}
}

main(){
	long n;
	F[0]=F[1]=1;
	while (cin >> n)
	cout << (n==0 ? 0 : f(n-1)) << endl;
}
