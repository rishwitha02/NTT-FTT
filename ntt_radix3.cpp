#include <bits/stdc++.h>
#include <omp.h>

using namespace std;
#define int long long
int root=362597;
int invroot=249134;
const int MOD = 472393;

int add(int a, int b) {
    return (a + b) % MOD;
}

int add(int a, int b,int c) {
    return add(add(a,b),c);
}

int sub(int a, int b) {
    return (a - b + MOD) % MOD;
}

int mul(int a, int b) {
    return (1LL * a * b) % MOD;
}

int mul(int a, int b, int c) {
    return mul(mul(a,b),c);
}

int powmod(int a, int b) {
    int res = 1;
    while (b) {
        if (b & 1) res = (res * a) % MOD;
        a = (a * a) % MOD;
        b >>= 1;
    }
    return res;
}

int divi(int a,int b)
{
    int ans=powmod(b,MOD-2);
    ans=mul(a,ans);
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
            if (powmod(r, phi/(*it)) == 1)
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

vector<int> ntt(vector<int>poly, int root) {
    int n = poly.size();
    if (n == 1) return poly;
    vector<int>a0(n/3),a1(n/3),a2(n/3);
    #pragma omp parallel for
    for (int i = 0, j = 0; i < n; i += 3, ++j) {
        a0[j] = poly[i];
        a1[j] = poly[i + 1];
        a2[j] = poly[i + 2];
    }
    int root3 = powmod(root, 3);
    #pragma omp task
    vector<int> y_0 = ntt(a0, root3);
    #pragma omp task
    vector<int> y_1 = ntt(a1, root3);
    #pragma omp task
    vector<int> y_2 = ntt(a2, root3);
    #pragma omp taskwait
    vector<int> y(n);
    int w = 1;
    for (int i = 0; i < n / 3; i++) {
        int u = y_0[i];
        int v = y_1[i];
        int x = y_2[i];
        int w2=mul(w,powmod(root,n/3));
        int w3=mul(w2,powmod(root,n/3));
        y[i] = add(u,mul(w,v),mul(w,w,x));
        y[(n/3)+i] = add(u,mul(w2,v),mul(w2,w2,x));
        y[(2*n/3)+i] = add(u,mul(w3,v),mul(w3,w3,x));
        w=mul(w,root);
    }
    return y;
}
vector<int> forward_ntt(vector<int>poly) {
    poly=ntt(poly, root);
    return poly;
}

vector<int> inverse_ntt(vector<int> poly) {
    poly = ntt(poly, invroot);
    int sz = poly.size();
    for (int i = 0; i < sz; i++) {
        poly[i] = divi(poly[i], sz);
    }
    return poly;
}

// Multiply two polynomials using NTT
vector<int> multiply_polynomials(vector<int>poly1, vector<int>poly2) {
    int max_degree = poly1.size() + poly2.size() - 1;
    int n = 1;
    while(n<max_degree)
    {
        n=n*3;
    }
    vector<int> padded_poly1(n), padded_poly2(n);
    for (int i = 0; i < poly1.size(); i++) padded_poly1[i] = poly1[i];
    for (int i = 0; i < poly2.size(); i++) padded_poly2[i] = poly2[i];
        
    // finding primitive root
    int prim_root=findPrimitive(MOD);
    root=powmod(prim_root,(MOD-1)/n);
    invroot=powmod(prim_root,MOD-1-((MOD-1)/n));

    // forward ntt
    vector<int> ntt1 = forward_ntt(padded_poly1);
    vector<int> ntt2 = forward_ntt(padded_poly2);

    // Pointwise multiplication and inverse NTT
    vector<int> product(n);
    for (int i = 0; i < n; i++) {
        product[i] = mul(ntt1[i],ntt2[i]);
    }
    
    // for (int x : ntt1) cout << x << " ";
    // cout<<endl;
    // for (int x : ntt2) cout << x << " ";
    // cout<<endl;
    // for (int x : product) cout << x << " ";
    // cout<<endl;

    // inverse ntt
    vector<int> result = inverse_ntt(product);
    return vector<int>(result.begin(), result.begin() + max_degree);
}

signed main() {
    // finding modulus
    
    // int ans=pow(3LL,9);
    // for(int i=1;i<=27;i++)
    // {
    //     if(isPrime(ans+1)>0)
    //     {
    //         cout<<ans<<endl;
    //     }
    //     ans=ans*i;
    // }

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
    // cout << "Resultant:\n";
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
    vector<int> result = multiply_polynomials(poly1,poly2);
    cout << "Resultant polynomial coefficients:\n";
    for (int x : result) cout << x << " ";
    cout<<endl;

    // input number as string

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
    // vector<int> result = multiply_polynomials(poly1,poly2);
    // cout << "Resultant polynomial coefficients:\n";
    // for (int x : result) cout << x << " ";
    // cout<<endl;

    return 0;
}