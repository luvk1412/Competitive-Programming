#include <bits/stdc++.h>
#define ll long long
#define ull unsigned long long int
#define pb push_back
#define mp make_pair
#define fi first
#define se second
#define p(x) printf("%d\n", x)
#define pl(x) printf("%lld\n", x)
#define s(x) scanf("%d", &x)
#define sl(x) scanf("%lld", &x)
#define sf(x) scanf("%lf", &x)
#define INF 1e18+9
#define endl '\n'
#define FIO ios_base::sync_with_stdio(false)
using namespace std;
#define M 1000000007
#define MAX 1000001

time spent(){ 
	clock_t begin = clock();
	// CODE
	clock_t end = clock();
	double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
}

void ss(string &i){
	int temp=getchar_unlocked();
	while((temp < 'a' || temp > 'z') && (temp < 'A' || temp > 'Z'))
		temp=getchar_unlocked();
	while((temp >= 'a' && temp <= 'z') || (temp >= 'A' && temp <= 'Z')){
		i.pb((char)temp);
		temp = getchar_unlocked();
	}
}

void fs(ll &x){
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
 d	return result;
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

// segmet tree

void build(int a[], int v, int tl, int tr){
	if(tl == tr){
		t[v] = a[tl];
	}
	else{
		int tm = (tl+tr) >> 1;
		build(a, v << 1, tl, tm);
		build(a, (v << 1)+1, tm+1, tr);
		t[v] = t[v << 1] + t[(v << 1) + 1];
	}
}
ll sum(int v, int tl, int tr, int l, int r){
	if(l > r){
		return 0;
	}
	if(tl == l && tr == r){
		return t[v];
	}
	int tm = (tl + tr) >> 1;
	return sum(v << 1, tl, tm, l, min(tm, r)) + sum((v << 1)+1, tm + 1, tr, max(l, tm+1), r);
}
void update(int v, int tl, int tr, int pos, int val){
	if(tl == tr){
		t[v] = val;
	}
	else{
		int tm = (tl+tr) >> 1;
		if(pos <= tm){
			update(v << 1, tl, tm, pos, val);
		}
		else{
			update((v << 1) + 1, tm+1, tr, pos, val);
		}
		t[v] = t[v << 1] + t[(v << 1)+1];
	}
}

// DSU

int parent[MAX], rnk[MAX];
void make(int v){
	parent[v] = v;
	rnk[v] = 0;
}

int find(int v){
	if(parent[v] != v)
		parent[v] = find(parent[v]);
	return parent[v];
}
void Union(int a, int b){
	a = find(a);
	b = find(b);
	if(a != b){
		if(rnk[a] < rnk[b])
			swap(a, b);
		parent[b] = a;
		if(rnk[a] == rnk[b])
			rnk[a]++;
	}
}
