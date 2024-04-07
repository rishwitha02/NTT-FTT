#include <iostream>
#include <vector>

using namespace std;

long long mod_inverse(long long a, long long m) {
    a = a % m;
    for (int x = 1; x < m; x++) {
        if ((a * x) % m == 1) {
            return x;
        }
    }
    return 1;
}

long long chinese_remainder_theorem(int n, vector<long long>& a, vector<long long>& mod) {
    long long M = 1;
    for (int i = 0; i < n; i++) {
        M *= mod[i];
    }

    long long x = 0;
    for (int i = 0; i < n; i++) {
        long long Mi = M / mod[i];
        long long Mi_inv = mod_inverse(Mi, mod[i]);
        x += a[i] * Mi * Mi_inv;
    }

    return x % M;
}

int main() {
    int n;
    cin >> n;
    vector<long long> a(n);
    vector<long long> mod(n);
    for (int i = 0; i < n; i++) {
        cin >> a[i];
        cin >> mod[i];
    }
    long long x = chinese_remainder_theorem(n, a, mod);
    cout << "The value of x that satisfies all the equations is: " << x << endl;

    return 0;
}

/*

Test #1:
In:
7
1 2
1 3
2 5
1 7
4 11
6 13
8 17
Out:
300007

Test #2:
3
1 2
1 3
3 17
Out:
37

Test #3:
1
100 106
Out:
100

Test #4:
2
2 3
3 5
Out:
8
*/
