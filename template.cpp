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
#define bitcount __builtin_popcountll
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
		}
		a = multiply(a, a, m);
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

//    MST
void merge(int v, int vl, int vr){
	int i = 0, j = 0;
	for(i = 0, j = 0; i < t[vl].size() && j < t[vr].size();){
		if(t[vl][i] <= t[vr][j]){
			t[v].pb(t[vl][i]);
			++i;
		}
		else{
			t[v].pb(t[vr][j]);
			++j;
		}
	}
	while(i < t[vl].size()){
		t[v].pb(t[vl][i]);
		++i;
	}
	while(j < t[vr].size()){
		t[v].pb(t[vr][j]);
		++j;
	}
}
void build(int v, int tl, int tr){
	if(tl == tr){
		t[v].pb(a[tl]);
	}
	else{
		int tm = (tl+tr) >> 1;
		build(v << 1, tl, tm);
		build((v << 1)+1, tm+1, tr);
		merge(v, v << 1, (v << 1) + 1);
	}
}
int qry(int v, int tl, int tr, int l, int r, int k){
	if(l > r){
		return 0;
	}
	if(tl == l && tr == r){
		auto it = upper_bound(t[v].begin(), t[v].end(), k);
		return (t[v].end()-it);
	}
	int tm = (tl + tr) >> 1;
	return qry(v << 1, tl, tm, l, min(tm, r), k) + qry((v << 1)+1, tm + 1, tr, max(l, tm+1), r, k);
}

// LCA
#define MLOG 18
int parent[MAX][MLOG], level[MAX], tin[MAX], tout[MAX], timer, vis[MAX];
vector <int> g[MAX];
void dfs(int node){
	vis[node] = 1;
	tin[node] = ++timer;
	for(int i : g[node]){
		if(vis[i] == 0){
			parent[i][0] = node;
			level[i] = level[node]+1;
			dfs(i);
		}
	}
	tout[node] = ++timer;
}
void setparent(int n){
	for(int i = 1 ; i < MLOG ; ++i){
		for(int j = 1 ; j <= n; ++j){
			parent[j][i] = parent[parent[j][i-1]][i-1];
		}
	}
}
bool isanc(int top, int bot){
	return (tin[top] <= tin[bot]) && (tout[bot] <= tout[top]);
}
int lca(int a, int b){
	int diff = level[a] - level[b];
	if(diff < 0){
		swap(a, b);
		diff *= -1;
	}
	for(int i = 0; i < MLOG; ++i){
		if(diff & (1LL<<i)){
			a = parent[a][i];
		}
	}
	if(a != b){
		for(int i = MLOG-1; i >= 0; --i){
			if(parent[a][i] != parent[b][i]){
				a = parent[a][i];
				b = parent[b][i];
			}
		}
		a = parent[a][0];
	}
	return a;
} 
int dist(int a , int b){
	return level[a] + level[b] - 2 * level[lca(a , b)];
}

//  MATRIX EXPO

void matmul(ll a[3][3], ll b[3][3], ll c[3][3]){
	int i, j, k;
	for(i = 0; i < 3; ++i){
		for(j = 0; j < 3; ++j){
			c[i][j] = 0;
			for(k = 0; k < 3; ++k){
				c[i][j] += (a[i][k]*b[k][j];) % M;
				c[i][j] %= M;
			}
		}
	}
}

void matexp(ll a[3][3], ll b, ll ans[3][3]){
	int i, j;
	ll tmp[3][3];
	for(i = 0; i < 3; ++i){
		for(j = 0; j < 3; ++j){
			if(i == j){
				ans[i][j] = 1;
			}
			else{
				ans[i][j] = 0;
			}
		}
	}
	while(b > 0){
		if(b&1){
			matmul(ans, a, tmp);
			for(i = 0; i < 3; ++i)
				for(j = 0; j < 3; ++j)
					ans[i][j] = tmp[i][j];
		}
		matmul(a, a, tmp);
		for(i = 0; i < 3; ++i)
			for(j = 0; j < 3; ++j)
				a[i][j] = tmp[i][j];
		b >>= 1;
	}
}
// KMP
void LPS(string &p){
	lps[0] = 0;
	int i = 1, j = 0;
	while(i < p.size()){
		if(p[i] == p[j]){
			j++;
			lps[i] = j;
			i++;
		}
		else if(j > 0){
			j = lps[j-1];
		}
		else{
			lps[i] = 0;
			i++;
		}
	}
}
// TRIE
struct TRIE{
	int next[MAX][13];
	int end[MAX];
	int sz;
	void clear(){
		memset(next, -1, sizeof(next));
		memset(end, 0, sizeof(end));
		sz = 0;
	}
	int insert(string &s){
		int v = 0, i;
		for(i = 0; i < s.size(); ++i){
			if(next[v][s[i] - 'a'] == -1){
				next[v][s[i] - 'a'] = ++sz;
			}
			v = next[v][s[i] - 'a'];
		}
		end[v]++;
	}
	bool search(string &s){
		int v = 0, i;
		for(i = 0; i < s.size(); ++i){
			if(next[v][s[i] - 'a'] == -1)
				return false;
			v = next[v][s[i] - 'a'];
		}
		return (end[v] > 0);
	}
}tr;
