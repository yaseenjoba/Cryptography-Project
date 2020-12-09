#include<iostream>
#include<vector>
#include<bits/stdc++.h>
#include<bitset>
#include<algorithm>
#include<math.h>
#include<cmath>
#include <stdio.h>
#include "Cephes.h"
#define debug(x )cout<<'['<<#x<<" is "<<x<<"]"<<endl;
using namespace std;
double normalCDF(double x) // Phi(-∞, x) aka N(x)
{
	return erfc(-x / std::sqrt(2)) / 2;
}
void Copy(vector<int>& v1, vector<int>& v2) {
	v2 = v1;
}
//this function for test 10.(additional test)
int BerlekampMassey(vector<int>& s)
{
	int L, N, m, d;
	int n = s.size();
	vector<int>c(n), b(n), t(n);

	//Initialization
	b[0] = c[0] = 1;
	N = L = 0;
	m = -1;

	//Algorithm core
	while (N < n)
	{
		d = s[N];
		for (int i = 1; i <= L; i++)
			d ^= c[i] & s[N - i];            //(d+=c[i]*s[N-i] mod 2)
		if (d == 1)
		{
			Copy(c, t);    //T(D)<-C(D)
			for (int i = 0; (i + N - m) < n; i++)
				c[i + N - m] ^= b[i];
			if (L <= (N >> 1))
			{
				L = N + 1 - L;
				m = N;
				Copy(t, b);    //B(D)<-T(D)
			}
		}
		N++;
	}
	return L;
}


int main() {
	vector<int> key = { 65,66,67,68,69,90,150,200,14/*,23,78,233,44*/ };
	vector<int>s(256);
	vector<int>t(256);

	for (int i = 0; i < 256; i++) {
		s[i] = i;
		t[i] = key[i % key.size()];
	}
	int temp = 0;
	for (int i = 0; i < 256; i++) {
		temp = (temp + s[i] + t[i]) % 256;
		swap(s[i], s[temp]);
	}

	int i = 0, j = 0;
	vector<int>gkey;
	while (gkey.size() != 125000) {
		i = (i + 1) % 256;

		j = (j + s[i]) % 256;

		swap(s[i], s[j]);

		int t = (s[i] + s[j]) % 256;
		gkey.push_back(s[t]);
	}

	//convert the key from bytes to bits.
	string tempForConvert = "";
	for (int i = 0; i < gkey.size(); i++)
	{
		tempForConvert += bitset<8>(gkey[i]).to_string();
	}
	vector<int>finalKey(tempForConvert.size());
	for (int i = 0; i < tempForConvert.size(); i++)
	{
		finalKey[i] = tempForConvert[i] - '0';
	}

	// Test number 1 : Frequency (Monobit) Test.
	vector<int>forTest1 = finalKey;
	//initialize Sn value.
	int Sn = 0;
	for (int i = 0; i < forTest1.size(); i++)
	{
		forTest1[i] = 2 * forTest1[i] - 1;
		Sn += forTest1[i];
	}
	//compute Sobs value.

	double Sobs = abs(Sn) / sqrt(forTest1.size());
	long double P_value = erfc(Sobs / sqrt(2));
	if (P_value < 0.01) {
		cout << "Test 1 (Frequency (Monobit) Test)  : Failed " << endl;
	}
	else cout << "Test 1 (Frequency (Monobit) Test)  : Passed" << endl;

	//Test number 2  : Frequency Test within a Block.
	// set the length of each block (M) = 100;
	int M = 10000;
	vector<int>forTest2 = finalKey;
	vector<double>PiOfI;
	int N = forTest2.size() / M;
	for (int i = 1; i < forTest2.size(); i += M)
	{
		double pi = 0;
		for (int j = 1; j < i + M; j++)
		{
			pi += ((forTest2[(i % M) - 1] * M) + j);
		}
		pi /= M;
		PiOfI.push_back(pi);
	}
	//cout << PiOfI.size() << endl;

	//find x^2(obs).
	long double Xops = 0;
	for (int i = 0; i < N; i++)
	{
		Xops += pow(PiOfI[i] - 0.5, 2);
	}

	Xops *= (4 * M);
	P_value = Cephes::cephes_igamc(N / 2, Xops / 2);

	if (P_value < 0.01) {
		cout << "Test 2 (Frequency Test within a Block)  : Failed " << endl;
	}
	else cout << "Test 2 (Frequency Test within a Block)  : Passed" << endl;
	//Test number 3 : Runs Test.
	//calculate Pi.
	vector<int>forTest3 = finalKey;
	double pi = 0;
	for (int i = 0; i < forTest3.size(); i++)
	{
		pi += forTest3[i];
	}
	pi /= forTest3.size();
	//Determine if the prerequisite Frequency test is passed(and its passed).
	//calculate Vnobs..
	int Vnobs = 0;
	for (int i = 0; i < forTest3.size() - 1; i++)
	{
		Vnobs += (forTest3[i] == forTest3[i + 1]) ? 0 : 1;
	}
	Vnobs++;
	double n = forTest3.size();

	P_value = erfc(abs(Vnobs - 2 * n * pi * (1 - pi)) / (2 * sqrt(2 * n) * pi * (1 - pi)));

	if (P_value < 0.01) {
		cout << "Test 3 (Runs Test)  : Failed " << endl;
	}
	else cout << "Test 3 (Runs Test)  : Passed" << endl;

	// Test number 13(Cumulative Sums (Cusum) Test) using mode 0
	vector<int>forTest13 = finalKey;
	for (int i = 0; i < forTest13.size(); i++)
	{
		forTest13[i] = 2 * forTest13[i] - 1;
	}
	for (int i = 1; i < forTest13.size(); i++) {
		forTest13[i] += forTest13[i - 1];
	}
	int z = INT_MIN;
	for (int i = 0; i < (int)forTest13.size(); i++) {
		z = max(z, abs(forTest13[i]));
	}
	P_value = 1;

	int k1 = int(((-1 * n / z) + 1) / 4);
	int k2 = ((n / z) - 1) / 4;
	n = forTest13.size();
	double temp1 = 0;
	for (int i = k1; i <= k2; i++)
	{
		temp1 += normalCDF(((4 * i + 1) * z) / sqrt(n)) - normalCDF(((4 * i - 1) * z) / sqrt(n));

	}
	double temp2 = 0;
	k1 = ((-1 * n / z) - 3) / 4;
	for (int i = k1; i <= k2; i++)
	{
		temp2 += normalCDF(((4 * i + 3) * z) / sqrt(n)) - normalCDF(((4 * i + 1) * z) / sqrt(n));

	}
	P_value -= temp1;
	P_value += temp2;
	if (P_value < 0.01) {
		cout << "Test 13 (Cumulative Sums (Cusum) Test)  : Failed " << endl;
	}
	else cout << "Test 13 (Cumulative Sums (Cusum) Test)  : Passed" << endl;

	//Test 15 :  Random Excursions Variant Test.
	vector<int>Test = finalKey;
	for (int i = 0; i < Test.size(); i++)
	{
		int t = (2 * Test[i] - 1);
		Test[i] = t;

	}

	for (int i = 1; i < Test.size(); i++)
	{
		Test[i] += Test[i - 1];
	}
	Test.push_back(0);
	vector<int>forTest_15 = {};
	forTest_15.push_back(0);
	for (int i = 0; i < Test.size(); i++)
	{
		forTest_15.push_back(Test[i]);
	}
	//ceil(forTest_15.size() / 4.0)
	int J = 1490;
	map<int, int>freq;//ξ (x)
	for (int i = 0; i < forTest_15.size(); i++)
	{
		freq[forTest_15[i]]++;
	}
	map<int, bool>isCalc;
	vector<double>P;
	for (int i = 0; i < forTest_15.size(); i++)
	{
		if (forTest_15[i] <= 9 && forTest_15[i] >= -9 && !isCalc[forTest_15[i]]) {
			isCalc[forTest_15[i]] = true;
			double p = erfc(abs(freq[forTest_15[i]] - J) / sqrt(2 * J * (4 * abs(forTest_15[i]) - 2)));
			P.push_back(p);
		}
	}
	bool isRandom = true;
	for (int i = 0; i < P.size(); ++i) {
		if (P[i] < 0.01) {
			isRandom = false;
			break;
		}
	}
	if (!isRandom) {
		cout << "Test 15 (Random Excursions Variant Test)  : Failed " << endl;
	}
	else cout << "Test 15 (Random Excursions Variant Test)  : Passed" << endl;
	//(additional test)
	//Test number 10  :   Linear Complexity Test (with  m =10000 and usign Berlekamp-Massey algorithm).
	vector<int>forTest10 = finalKey;
	M = 10;
	double meanValue = (M / 2) + ((9 + pow(-1, M + 1)) / 36) + (((M / 3) + (2 / 9)) / pow(2, M));
	vector<int>temp13;
	vector<double>T;
	for (int i = 0; i < forTest13.size(); i += M)
	{
		temp13.clear();
		for (int j = i; j < i + M; j++)
		{
			temp13.push_back(forTest13[j]);
		}

		T.push_back(pow(-1, M) * (BerlekampMassey(temp13) - meanValue) + 2.0 / 9.0);
	}
	vector<int>v(7);
	v.assign(7, 0);
	for (int i = 0; i < T.size(); i++)
	{
		if (T[i] <= -2.5) {
			v[0]++;
		}
		else if (T[i] > -2.5 && T[i] <= -1.5) {
			v[1]++;
		}
		else if (T[i] > -1.5 && T[i] <= -0.5) {
			v[2]++;
		}
		else if (T[i] > -0.5 && T[i] <= 0.5) {
			v[3]++;
		}
		else if (T[i] > 0.5 && T[i] <= 1.5) {
			v[4]++;
		}
		else if (T[i] > 1.5 && T[i] <= 2.5) {
			v[5]++;
		}
		else v[6]++;
	}
	vector<double>tempPi(7);
	tempPi[0] = 0.010417;
	tempPi[1] = 0.03125;
	tempPi[2] = 0.125;
	tempPi[3] = 0.5;
	tempPi[4] = 0.25;
	tempPi[5] = 0.0625;
	tempPi[6] = 0.020833;
	Xops = 0;
	int k = 6;
	N = forTest13.size() / M;
	for (int i = 0; i <= k; i++)
	{
		Xops += (pow(v[i] - (N * tempPi[i]), 2)) / (N * tempPi[i]);
	}

	P_value = Cephes::cephes_igamc(Xops / 2, double(k / 2.0));

	P_value *= tgamma(double(k / 2.0));
	P_value /= tgamma(Xops / 2);

	if (P_value < 0.01) {
		cout << "Test 10 (Linear Complexity Test)  : Faild" << endl;
	}
	else cout << "Test 10 (Linear Complexity Test)  : Passed" << endl;


	return 0;
}

