#include <bits/stdc++.h>
#include <omp.h>

using namespace std;
#define int long long

int add(int a, int b, int MOD) {
    return (a + b) % MOD;
}

int sub(int a, int b, int MOD) {
    return (a - b + MOD) % MOD;
}

int mul(int a, int b, int MOD) {
    return (1LL * a * b) % MOD;
}

int powmod(int a, int b, int MOD) {
    int res = 1;
    while (b) {
        if (b & 1) res = (res * a) % MOD;
        a = (a * a) % MOD;
        b >>= 1;
    }
    return res;
}

int divi(int a,int b, int MOD)
{
    int ans=powmod(b,MOD-2, MOD);
    ans=mul(a,ans, MOD);
    return ans;
}

bool isPrime(int n)
{
    if (n <= 1) return false;
    if (n <= 3) return true;
    if (n%2 == 0 || n%3 == 0) return false;
    for (int i=5; i*i<=n; i=i+6)
        if (n%i == 0 || n%(i+2) == 0)
            return false;
    return true;
}

// Utility function to store prime factors of a number
void findPrimefactors(unordered_set<int> &s, int n)
{
    while (n%2 == 0)
    {
        s.insert(2);
        n = n/2;
    }
    // n must be odd at this point. So we can skip one element (Note i = i +2)
    for (int i = 3; i <= sqrt(n); i = i+2)
    {
        // While i divides n, print i and divide n
        while (n%i == 0)
        {
            s.insert(i);
            n = n/i;
        }
    }
    // This condition is to handle the case when n is a prime number greater than 2
    if (n > 2)
        s.insert(n);
}

// Function to find smallest primitive root of n
int findPrimitive(int MOD)
{
    unordered_set<int> s;

    // Check if MOD is prime or not
    if (isPrime(MOD)==false)
        return -2;

    // Find value of Euler Totient function of MOD. Since MOD is a prime number, the value of Euler totient function is MOD-1 as there are MOD-1 relatively prime numbers.
    int phi = MOD-1;

    // Find prime factors of phi and store in a set
    findPrimefactors(s, phi);

    // Check for every number from 2 to phi
    for (int r=2; r<=phi; r++)
    {
        // Iterate through all prime factors of phi and check if we found a power with value 1
        bool flag = false;
        for (auto it = s.begin(); it != s.end(); it++)
        {
            // Check if r^((phi)/primefactors) mod MOD is 1 or not
            if (powmod(r, phi/(*it), MOD) == 1)
            {
                flag = true;
                break;
            }
        }

        // If there was no power with value 1.
        if (flag == false)
        return r;
    }

    // If no primitive root found
    return -1;
}

vector<int> ntt(vector<int>poly, int root, int MOD) {
    int n = poly.size();
    if (n == 1) return poly;
    vector<int>a0(n/2),a1(n/2);
    #pragma omp parallel for
    for (int i = 0, j = 0; i < n; i += 2, ++j) {
        a0[j] = poly[i];
        a1[j] = poly[i + 1];
    }
    #pragma omp task
    vector<int> y_even = ntt(a0, powmod(root, 2, MOD), MOD);
    #pragma omp task
    vector<int> y_odd = ntt(a1, powmod(root, 2, MOD), MOD);
    #pragma omp taskwait
    vector<int> y(n);
    int w = 1;
    for (int i = 0; i < n / 2; i++) {
        int u = y_even[i];
        int v = y_odd[i];
        int k= mul(w,v,MOD);
        y[i] = add(u,k,MOD);
        y[i + n / 2] = sub(u,k,MOD);
        w = (w * root) % MOD;
    }
    return y;
}
vector<int> forward_ntt(vector<int>poly, int MOD, int root) {
    poly=ntt(poly, root, MOD);
    return poly;
}

vector<int> inverse_ntt(vector<int>poly, int MOD, int invroot) {
    poly=ntt(poly, invroot, MOD);
    int sz=poly.size();
    for(int i=0;i<sz;i++)
    {
        poly[i]=divi(poly[i],sz, MOD);
    }
    return poly;
}

// Multiply two polynomials using NTT
vector<int> multiply_polynomials(vector<int>poly1, vector<int>poly2,int MOD) {
    int max_degree = poly1.size() + poly2.size() - 1;
    int n = 1 << (32 - __builtin_clz(max_degree));  
    vector<int> padded_poly1(n), padded_poly2(n);
    for (int i = 0; i < poly1.size(); i++) padded_poly1[i] = poly1[i];
    for (int i = 0; i < poly2.size(); i++) padded_poly2[i] = poly2[i];

    // finding primitive root
    int prim_root=findPrimitive(MOD);
    int root=powmod(prim_root,(MOD-1)/n, MOD);
    int invroot=powmod(prim_root,MOD-1-((MOD-1)/n), MOD);

    // cout<<prim_root<<" "<<root<<" "<<invroot<<endl;
    // int ans=1;
    // for(int i=0;i<=MOD;i++)
    // {
    //     cout<<ans<<" ";
    //     ans=mul(ans,prim_root,MOD);
    // }
    // cout<<endl;

    // forward ntt
    vector<int> ntt1 = forward_ntt(padded_poly1,MOD,root);
    vector<int> ntt2 = forward_ntt(padded_poly2,MOD,root);

    // Pointwise multiplication and inverse NTT
    vector<int> product(n);
    for (int i = 0; i < n; i++) {
        product[i] = mul(ntt1[i],ntt2[i],MOD);
    }
    
    // for (int x : ntt1) cout << x << " ";
    // cout<<endl;
    // for (int x : ntt2) cout << x << " ";
    // cout<<endl;
    // for (int x : product) cout << x << " ";
    // cout<<endl;

    // inverse ntt
    vector<int> result = inverse_ntt(product,MOD,invroot);
    return vector<int>(result.begin(), result.begin() + max_degree);
}
int mod_inverse(int a, int m) {
    a = a % m;
    for (int x = 1; x < m; x++) {
        if ((a * x) % m == 1) {
            return x;
        }
    }
    return 1;
}
int crt(int a,int b,int MOD1,int MOD2)
{
    int M=(MOD1*MOD2);
    int ans=0;
    ans+=(a*(M/MOD1)*mod_inverse(M/MOD1,MOD1));
    ans+=(b*(M/MOD2)*mod_inverse(M/MOD2,MOD2));
    ans%=M;
    return ans;
}
vector<int>ntt_crt(vector<int>result1,vector<int>result2,int MOD1,int MOD2)
{
    int sz=result1.size();
    vector<int>result(sz);
    for(int i=0;i<sz;i++)
    {
        result[i]=crt(result1[i],result2[i],MOD1,MOD2);
    }
    return result;
}
signed main() {

    // Finding modulo

    // int ans=pow(2LL,8);
    // for(int i=1;i<=1000;i++)
    // {        
    //     if(isPrime(i*ans+1))
    //     {
    //         cout<<i*ans+1<<endl;
    //     }
    // }
    // 8 -> 257,769
    const int MOD1 = 257;
    const int MOD2 = 769;

    // input as numbers

    // vector<int> a,b;
    // int n,m;
    // cin>>n>>m;
    // while(n)
    // {
    //     a.push_back(n%10);
    //     n=n/10;
    // }
    // while(m)
    // {
    //     b.push_back(m%10);
    //     m=m/10;
    // }
    // reverse(a.begin(),a.end());
    // reverse(b.begin(),b.end());
    // int sz1=a.size();
    // int sz2=b.size();
    // vector<int>result = multiply_polynomials(a, b);
    // int ans=0;
    // int sz=result.size();
    // for(int i=0;i<(sz1+sz2-1);i++)
    // {
    //     ans=(ans*10)+result[i];
    // }    
    // cout<<ans<<endl;

    // input as vector of polynomials

    int n,m;
    cin>>n>>m;
    vector<int> poly1(n),poly2(m);
    for(int i=0;i<n;i++)
    {
        cin>>poly1[i];
    }
    for(int i=0;i<m;i++)
    {
        cin>>poly2[i];
    }
    vector<int> result1 = multiply_polynomials(poly1,poly2,MOD1);
    vector<int> result2 = multiply_polynomials(poly1,poly2,MOD2);
    vector<int>result=ntt_crt(result1,result2,MOD1,MOD2);
    for (int x : result1) cout << x << " ";
    cout<<endl;
    for (int x : result2) cout << x << " ";
    cout<<endl;
    cout << "Resultant polynomial coefficients:\n";
    for (int x : result) cout << x << " ";
    cout<<endl;

    // input as vectors

    // vector<int> poly1,poly2;
    // string s1,s2;
    // cin>>s1>>s2;
    // for(auto it:s1)
    // {
    //     poly1.push_back(it-'0');
    // }
    // for(auto it:s2)
    // {
    //     poly2.push_back(it-'0');
    // }
    // vector<int> result1 = multiply_polynomials(poly1,poly2,MOD1);
    // vector<int> result2 = multiply_polynomials(poly1,poly2,MOD2);
    // vector<int>result=ntt_crt(result1,result2,MOD1,MOD2);
    // cout << "Resultant polynomial coefficients:\n";
    // for (int x : result) cout << x << " ";
    // cout<<endl;

    // cout<<crt(2,3,3,5)<<endl;

    return 0;
}


/*

input:
4 4
111 111 111 111
111 111 111 111
output:
12321 24642 36963 49284 36963 24642 12321

*/