#include<iostream>
#include<string>
using namespace std;
int makeEqualLength(string &str1, string &str2)
{
	int len1 = str1.size();
	int len2 = str2.size();
	if (len1 < len2)
	{
		for (int i = 0 ; i < (len2 - len1) ; i++)
			str1 = '0' + str1;
		return len2;
	}
	else if (len1 > len2)
	{
		for (int i = 0 ; i < len1 - len2 ; i++)
			str2 = '0' + str2;
	}
	return len1; 
}
string addStrings(string first,string second)
{
	string result=""; 
	int len = makeEqualLength(first, second);
	int carry = 0;
	for (int i = len-1 ; i >= 0 ; i--)
	{
		int firstnum = first[i] - '0';
		int secondnum = second[i] - '0';
		int sum = (firstnum + secondnum + carry);
		carry=sum/10;
		sum=(sum%10);
		char ch=(sum+'0');
		result = ch + result;
	}
	if (carry) result = '1' + result;
	return result;
}
string subStrings(string first,string second)
{
	string result=""; 
	int len = makeEqualLength(first, second);
	for (int i = len-1 ; i >= 0 ; i--)
	{
		int firstnum = first[i] - '0';
		int secondnum = second[i] - '0';
		int diff = (firstnum - secondnum);
		char ch=(diff+'0');
		if(diff<0)
		{
			int j=i-1;
			while(j>=0)
			{
				first[j]=(((first[j]-'0')+10-1)%10)+'0';
				if(first[j]!='9')break;
				else j--;
			}
			ch=(10+diff)+'0';
		}
		result = ch + result;
	}
	return result;
}
string multiply(string x, string y)
{
	int n = makeEqualLength(x, y);
	if (n == 0) return 0;
	if (n == 1) return to_string((x[0] - '0')*(y[0] - '0'));
	int fh = n/2; 
	int sh = (n-fh);

	string a = x.substr(0, fh);
	string b = x.substr(fh, sh);
	string c = y.substr(0, fh);
	string d = y.substr(fh, sh);

	string ac = multiply(a, c);
	string bd = multiply(b, d);
	string e = multiply(addStrings(a, b), addStrings(c, d));
	string f = subStrings(e,addStrings(ac,bd));

	for (int i = 0; i < 2*sh; i++)
		ac.append("0");
	for (int i = 0; i < sh; i++)
		f.append("0");

	string ans = addStrings(addStrings(ac,f),bd);
    int st=0,sz=ans.size();
    for(int i=0;i<sz;i++)
    {
        if(ans[i]=='0')st++;
		else break;
    }
	return ans.substr(st,sz-st);
}
int main()
{
	int t;
	cout<<"Number of test cases: ";
	cin>>t;
	string ch;
	getline(cin, ch);
	while(t--)
	{
		string num1,num2;
		cout<<"First number: ";
		getline(cin, num1);
		cout<<"Second number: ";
		getline(cin, num2);
		// cout<<stoi(num1)*stoi(num2)<<endl;
		string ans=multiply(num1,num2);
		if(ans.empty())ans="0";
		cout<<ans<<endl;
	}
}
