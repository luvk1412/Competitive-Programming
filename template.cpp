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
#define trace1(x)                cerr<<#x<<": "<<x<<endl
#define trace2(x, y)             cerr<<#x<<": "<<x<<" | "<<#y<<": "<<y<<endl
#define trace3(x, y, z)          cerr<<#x<<":" <<x<<" | "<<#y<<": "<<y<<" |\
									 "<<#z<<": "<<z<<endl
#define trace4(a, b, c, d)       cerr<<#a<<": "<<a<<" | "<<#b<<": "<<b<<" |\
									 "<<#c<<": "<<c<<" | "<<#d<<": "<<d<<endl
#define trace5(a, b, c, d, e)    cerr<<#a<<": "<<a<<" | "<<#b<<": "<<b<<" |\
									 "<<#c<<": "<<c<<" | "<<#d<<": "<<d<<" | "<<#e<< ": "<<e<<endl
#define trace6(a, b, c, d, e, f) cerr<<#a<<": "<<a<<" | "<<#b<<": "<<b<<" | "<<#c<<": "<< c<<" |\
								 "<<#d<<": "<<d<<" | "<<#e<< ": "<<e<<" | "<<#f<<": "<<f<<endl




//DEBUG

void __print(int x) {cerr << x;}
void __print(long x) {cerr << x;}
void __print(long long x) {cerr << x;}
void __print(unsigned x) {cerr << x;}
void __print(unsigned long x) {cerr << x;}
void __print(unsigned long long x) {cerr << x;}
void __print(float x) {cerr << x;}
void __print(double x) {cerr << x;}
void __print(long double x) {cerr << x;}
void __print(char x) {cerr << '\'' << x << '\'';}
void __print(const char *x) {cerr << '\"' << x << '\"';}
void __print(const string &x) {cerr << '\"' << x << '\"';}
void __print(bool x) {cerr << (x ? "true" : "false");}

template<typename T, typename V>
void __print(const pair<T, V> &x) {cerr << '{'; __print(x.first); cerr << ','; __print(x.second); cerr << '}';}
template<typename T>
void __print(const T &x) {int f = 0; cerr << '{'; for (auto &i: x) cerr << (f++ ? "," : ""), __print(i); cerr << "}";}
void _print() {cerr << "]\n";}
template <typename T, typename... V>
void _print(T t, V... v) {__print(t); if (sizeof...(v)) cerr << ", "; _print(v...);}
#ifndef ONLINE_JUDGE
#define debug(x...) cerr << "[" << #x << "] = ["; _print(x)
#else
#define debug(x...)
#endif


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
void write(int &x){
	register char buffor[35];
	register int i=0;
	do{
		buffor[i++]=(x%10)+'0';
		x/=10;
	} while(x);
	i--;
	while(i>=0) putchar_unlocked(buffor[i--]);
	putchar_unlocked('\n');
}
void fs(int &x){
	register int c = getchar_unlocked();
	x = 0;
	for(; (c<48 || c>57); c = getchar_unlocked());
	for(; c>47 && c<58; c = getchar_unlocked()){
		x = (x<<1) + (x<<3) + c - 48;
	}
}
void fs(long long &x){
	register long long c = getchar_unlocked();
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



int binexp1(int a, int b, int m){
	int result = 1;
	while(b>0){
		if(b & 1){
			result = (result  * 1LL * a) % m;
		}
		a = (a * 1LL * a) % m;
		b >>= 1;
	}
	return result;
}

long long multiply(long long a, long long b, long long m){
	long long result = 0;
	while(b>0){
		if(b & 1){
			result = result + a;
			result %= m;
		}
		a <<= 1;
		a %= m;
		b >>= 1;
	}
	return result;
}

long long binexp2(long long a, long long b, long long m){
	long long result = 1;
	while(b>0){
		if(b & 1){
			result = multiply(result, a, m);
		}
		a = multiply(a, a, m);
		b >>= 1;
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

bool miller(long long n){
	if(n <=1 || n % 2 == 0){
		if(n != 2){
			return false;
		}
	}
	if(n == 2 || n == 3){
		return true;
	}
	long long d = n-1;
    while(d % 2 == 0){
        d /= 2;
    }
	long long a[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
	for(int i = 0; i < 12 && a[i] < n; ++i){
		long long temp = d;
        long long mod = binexp2(a[i], temp, n);
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

void update(int i, int x){
	while(i <= N){
		bit[i] += x;
		i = i + (i & (-i));
	}
}

int sum(int i){
	int ans = 0;
	while(i != 0 && i > 0){
		ans += bit[i];
		i = (i & (i-1));
	}
	return ans;
}

void rupdate(int a, int b, int val){
	update(a, val, bit1);
	update(b + 1, -val, bit1);
	update(a, val * (a-1), bit2);
	update(b + 1, -val * b, bit2);
}

ll rsum(int i){
	return sum(i, bit1) * i - sum(i, bit2);
}

void update(int i, int j, int val){
	int y1;
	while (i <= N){
		y1 = j;
		while (y1 <= N){
			bit[i][y1] += val;
			y1 += (y1 & (-y1));
		}
		i += (i & -i);
	}
}

int sum(int i,int j){
	int sum= 0;         
	while(i){
		int y1 = j;
		while(y1){
			sum += bit[i][y1];
			y1 -= (y1 & (-y1));
		}
		i -= (i & (-i));
	}
	return sum;
}

// segment tree

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
//lazy propagation

void build(int v, int tl, int tr){
	lazy[v] = 0;
	if(tl == tr){
		t[v] = 1LL<<color[mp[tl]];
	}
	else{
		int mid = (tl + tr) >> 1;
		build(v<<1, tl, mid);
		build((v<<1) + 1, mid+1, tr);
		t[v] = t[v<<1] | t[(v<<1) + 1];
	}
}

void update(int v, int tl, int tr, int l, int r, long long val){
	if(lazy[v]){
		if(tl == tr){
			t[v] = lazy[v];
		}
		else{
			t[v] = lazy[v];
			lazy[v<<1] = lazy[v];
			lazy[(v<<1) + 1] = lazy[v];
		}
		lazy[v] = 0;
	}
	if(tl > r || tr < l)
		return;
	if(l <= tl && r >= tr){
		if(tl == tr){
			t[v] = (1LL << val);
		}
		else{
			t[v] = (1LL << val);
			lazy[v<<1] = (1LL<<val);
			lazy[(v<<1) + 1] = (1LL << val);
		}
		return;
	}
	int mid = (tl + tr) >> 1;
	update(v<<1, tl, mid, l, r, val);
	update((v<<1) + 1, mid + 1, tr, l, r, val);
	t[v] = t[v<<1] | t[(v<<1) + 1];
}

long long qry(int v, int tl, int tr, int l, int r){
	if(lazy[v]){
		if(tl == tr)
			t[v] = lazy[v];
		else{
			t[v] = lazy[v];
			lazy[v<<1] = lazy[v];
			lazy[(v<<1) + 1] = lazy[v];
		}
		lazy[v] = 0;
	}
	if(tr < l or tl > r)
		return 0LL;
	if(tl>=l and tr<=r){
		return t[v];
	}
	int mid = (tl + tr) >> 1;
	return qry(v << 1, tl, mid, l, r) + qry((v<<1) + 1, mid + 1, tr, l, r);
}


// DSU

int parent[N], rnk[N];

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

// Z - Algorithm
void getZarr(string &str){
	int n = str.length();
	int L, R, k; 
	L = R = 0;
	for(int i = 1; i < n; ++i){
		if (i > R){
			L = R = i;
			while (R<n && str[R-L] == str[R])
				R++;
			Z[i] = R-L;
			R--;
		}
		else{
			k = i-L;
			if (Z[k] < R-i+1)
				Z[i] = Z[k];
			else{
				L = i;
				while (R<n && str[R-L] == str[R])
					R++;
				Z[i] = R-L;
				R--;
			}
		}
	}
}

// TRIE
struct TRIE{
	int next[MAX][13];
	int end[MAX];
	int sz;
	void clear(){
		memset(next, -1, sizeof(next));`
		memset(end, 0, sizeof(end));
		sz = 0;
	}
	void insert(string &s){
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

//Dinics max flow ALGO
struct edge{
	int from, to, cap, flow;
};
vector <int> g[MAX];
vector <edge> e;
int source, sink, level[MAX], ptr[MAX];

void add_edge(int from, int to, int cap){
	edge e1 = {from, to, cap, 0};
	edge e2 = {to, from, 0, 0};
	g[from].push_back(e.size());
	e.push_back(e1);
	g[to].push_back(e.size());
	e.push_back(e2);
}

int bfs(){
	memset(level, -1, sizeof(level));
	level[source] = 0;
	queue <int> q;
	q.push(source);
	int v;
	while(!q.empty()){
		v = q.front();
		q.pop();
		if(v == sink){
			return 1;
		}
		for(auto i : g[v]){
		int to = e[i].to;
			if(level[to] == -1 && e[i].flow < e[i].cap){
				level[to] = level[v] + 1;
				q.push(to);
			}
		}
	}
	return 0;
}

int dfs(int v, int flow){
	if(flow <= 0)
		return 0;
	if(v == sink)
		return flow;
	for(; ptr[v] < g[v].size(); ++ptr[v]){
		int id = g[v][ptr[v]], to = e[id].to;
		if(level[to] != level[v] + 1)
			continue;
		int pushed = dfs(to, min(flow, e[id].cap - e[id].flow));
		if(pushed){
			e[id].flow += pushed;
			e[id^1].flow -= pushed;
			return pushed;
		}
	}
	return 0;
}
int dinic(){
	int flow = 0, val;
	while(bfs()){
		memset(ptr, 0, sizeof(ptr));
		while(val = dfs(source, INF)){
			flow += val;
		}
	}
	return flow;
}

// MO's
int block, n;

struct node{
	int l, r, idx;
	long long ans;
}qry[N];
bool cmp(node a, node b){
	if((a.l / block) != b.l / block)
		return (a.l / block) < (b.l/block);
	else
		return (a.r < b.r);
}

bool cmp1(node a, node b){
	return (a.idx < b.idx);
}
int ans;
void add(int val){
	//update
}
void remove(int val){
	//deupdate
}
int main(){
	sort(qry, qry + q, cmp);
	for(int i = qry[0].l; i <= qry[0].r; ++i){
		add(a[i]);
	}
	int curr, curl;
	qry[0].ans = ans;
	curl = qry[0].l;
	curr = qry[0].r;
	for(int i = 1; i < q; ++i){
		while(curl < qry[i].l){
			remove(a[curl++]);
		}
		while(curr < qry[i].r){
			add(a[++curr]);
		}
		while(curl > qry[i].l){
			add(a[--curl]);
		}
		while(curr > qry[i].r){
			remove(a[curr--]);
		}
		qry[i].ans = ans;
	}
	sort(qry, qry + q, cmp1);
	for(int i = 0; i < q; ++i){
		printf("%lld\n", qry[i].ans);
	}
}

//ordered Set | policy based data structure

#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
#include <bits/stdc++.h>
 
using namespace __gnu_pbds;
using namespace std;
 
typedef tree<int, null_type, less<int>, rb_tree_tag, tree_order_statistics_node_update> ordered_set;


int main(){
	ordered_set X;
	X.insert(1);
	X.insert(2);
	X.insert(4);
	X.insert(8);
	X.insert(16);
	cout<<*X.find_by_order(1)<<endl; // 2
	cout<<X.order_of_key(4)<<endl;   // 2
}

If we want to get map but not the set, as the second argument type must be used mapped type.
Apparently, the tree supports the same operations as the set (at least I haven't any problems with them before),
but also there are two new features — it is find_by_order() and order_of_key(). 
The first returns an iterator to the k-th largest element (counting from zero), the second — the
number of items in a set that are strictly smaller than our item. Example of use:


Bigint : https://ideone.com/TiMsbK
