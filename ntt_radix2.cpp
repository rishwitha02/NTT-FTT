#include <bits/stdc++.h>
#include <omp.h>

using namespace std;
#define int long long
int root;
int invroot;
const int MOD = 998244353;

int add(int a, int b) {
    return (a + b) % MOD;
}

int sub(int a, int b) {
    return (a - b + MOD) % MOD;
}

int mul(int a, int b) {
    return (1LL * a * b) % MOD;
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
    vector<int>a0(n/2),a1(n/2);
    #pragma omp parallel for
    for (int i = 0, j = 0; i < n; i += 2, ++j) {
        a0[j] = poly[i];
        a1[j] = poly[i + 1];
    }
    #pragma omp task
    vector<int> y_even = ntt(a0, powmod(root, 2));
    #pragma omp task
    vector<int> y_odd = ntt(a1, powmod(root, 2));
    #pragma omp taskwait
    vector<int> y(n);
    int w = 1;
    for (int i = 0; i < n / 2; i++) {
        int u = y_even[i];
        int v = y_odd[i];
        int k= mul(w,v);
        y[i] = add(u,k);
        y[i + n / 2] = sub(u,k);
        w = (w * root) % MOD;
    }
    return y;
}
vector<int> forward_ntt(vector<int>poly) {
    poly=ntt(poly, root);
    return poly;
}

vector<int> inverse_ntt(vector<int>poly) {
    poly=ntt(poly, invroot);
    int sz=poly.size();
    for(int i=0;i<sz;i++)
    {
        poly[i]=divi(poly[i],sz);
    }
    return poly;
}

// Multiply two polynomials using NTT
vector<int> multiply_polynomials(vector<int>poly1, vector<int>poly2) {
    int max_degree = poly1.size() + poly2.size() - 1;
    int n = 1 << (32 - __builtin_clz(max_degree));  
    vector<int> padded_poly1(n), padded_poly2(n);
    for (int i = 0; i < poly1.size(); i++) padded_poly1[i] = poly1[i];
    for (int i = 0; i < poly2.size(); i++) padded_poly2[i] = poly2[i];

    // finding modulo

        
    // finding primitive root
    int prim_root=findPrimitive(MOD);
    root=powmod(prim_root,(MOD-1)/n);
    invroot=powmod(prim_root,MOD-1-((MOD-1)/n));
    // cout<<root<<" "<<invroot<<endl;

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
    // vector<int> result = multiply_polynomials(poly1,poly2);
    // cout << "Resultant polynomial coefficients:\n";
    // for (int x : result) cout << x << " ";
    // cout<<endl;

    return 0;
}


/*

input:
123456789987654321123456789987654321123456789987654321123456789987654321123456789987654321123456789987654321123456789987654321123456789987654321
123456789987654321123456789987654321123456789987654321123456789987654321123456789987654321123456789987654321123456789987654321123456789987654321

output:
Resultant polynomial coefficients:
1 4 10 20 35 56 84 120 165 218 276 336 395 450 498 536 561 570 563 544 518 490 465 448 444 458 495 556 636 728 825 920 1006 1076 1123 1140 1125 1084 1026 960 895 840 804 796 825 894 996 1120 1255 1390 1514 1616 1685 1710 1687 1624 1534 1430 1325 1232 1164 1134 1155 1232 1356 1512 1685 1860 2022 2156 2247 2280 2249 2164 2042 1900 1755 1624 1524 1472 1485 1570 1716 1904 2115 2330 2530 2696 2809 2850 2811 2704 2550 2370 2185 2016 1884 1810 1815 1908 2076 2296 2545 2800 3038 3236 3371 3420 3373 3244 3058 2840 2615 2408 2244 2148 2145 2246 2436 2688 2975 3270 3546 3776 3933 3990 3935 3784 3566 3310 3045 2800 2604 2486 2475 2584 2796 3080 3405 3740 4054 4316 4495 4560 4495 4316 4054 3740 3405 3080 2796 2584 2475 2486 2604 2800 3045 3310 3566 3784 3935 3990 3933 3776 3546 3270 2975 2688 2436 2246 2145 2148 2244 2408 2615 2840 3058 3244 3373 3420 3371 3236 3038 2800 2545 2296 2076 1908 1815 1810 1884 2016 2185 2370 2550 2704 2811 2850 2809 2696 2530 2330 2115 1904 1716 1570 1485 1472 1524 1624 1755 1900 2042 2164 2249 2280 2247 2156 2022 1860 1685 1512 1356 1232 1155 1134 1164 1232 1325 1430 1534 1624 1687 1710 1685 1616 1514 1390 1255 1120 996 894 825 796 804 840 895 960 1026 1084 1125 1140 1123 1076 1006 920 825 728 636 556 495 458 444 448 465 490 518 544 563 570 561 536 498 450 395 336 276 218 165 120 84 56 35 20 10 4 1 


*/

/*

input:
1234
5678

output:
Resultant polynomial coefficients:
5 16 34 60 61 52 32 

*/