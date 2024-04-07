#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

using namespace std;

typedef complex<double> Complex;
double PI = acos(-1);
void fft(vector<Complex>& a, bool invert) {
    int n = a.size();
    if (n <= 1) return;

    vector<Complex> a0(n / 2), a1(n / 2);
    for (int i = 0, j = 0; i < n; i += 2, ++j) {
        a0[j] = a[i];
        a1[j] = a[i + 1];
    }

    fft(a0, invert);
    fft(a1, invert);

    double angle = (invert ? -1 : 1) * 2 * PI / n;
    Complex w(1), wn(cos(angle), sin(angle));

    for (int i = 0; i < n / 2; ++i) {
        Complex t = w * a1[i];
        a[i] = a0[i] + t;
        a[i + n / 2] = a0[i] - t;
        w *= wn;
    }
}

void fft(vector<Complex>& a) {
    fft(a, false);
}

void ifft(vector<Complex>& a) {
    fft(a, true);
    for (Complex& x : a) {
        x /= a.size();
    }
}

vector<int> multiply(const vector<int>& a, const vector<int>& b) {
    // Convert integers to Complex numbers for FFT
    vector<Complex> fa(a.begin(), a.end()), fb(b.begin(), b.end());
    int n = 1;
    while (n < (a.size()+b.size())) n <<= 1;
    fa.resize(n);
    fb.resize(n);

    // Perform FFT
    fft(fa);
    fft(fb);

    // Pointwise multiplication in frequency domain
    for (int i = 0; i < n; ++i) {
        fa[i] *= fb[i];
    }

    // Perform inverse FFT
    ifft(fa);

    // Extract real parts and round to integers
    vector<int> result(n);
    for (int i = 0; i < n; ++i) {
        result[i] = round(fa[i].real());
    }

    return result;
}

int main() {
    int n,m;
    cin>>n>>m;
    vector<int>a(n),b(m);
    for(int i=0;i<n;i++)
    {
        cin>>a[i];
    }
    for(int i=0;i<m;i++)
    {
        cin>>b[i];
    }
    vector<int> result = multiply(a, b);
    int carry = 0;
    int sz=n+m-1;
    for (int i = (sz-1); i>=0; --i) {
        result[i] += carry;
        carry = result[i] / 10;
        result[i] %= 10;
    }
    cout << "Multiplication Result: ";
    if(carry)
    {
        cout<<1;
    }
    for(int i = 0; i<sz; i++)
    {
        cout<<result[i];
    }
    return 0;
}
