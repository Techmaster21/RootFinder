/* Everything works in the usual cases, and the input is totally variable. 

////
x^7−14*x^5+49*x^3−36*x - pasted input doesn't work due to long -'s 
////
////
equation x^11 - 
2x^10-10010.5370709241x^9 + 20023.0741418482x^8 + 105382.651717972x^7 - 
230788.3775777922x^6 - 119407.103358567x^5 + 469602.584294926x^4 - 
176669.430017213x^3 - 116263.7242605x^2 + 58131.86213025x
Comes up with only 9 real roots when in fact there are 11. -- due to zero of 
multiplicity two, sturm's theorem is not supposed to work properly with zeros 
of higher multiplicity than one --- find a fix.
Also, the root that is to be zero is not exactly zero -- in fact, its not 
terribly close. It's only e-15. Might this be an issue for your nonzero - 
leading term checker for remainder cleanup as well? How close to zero is zero?
////
////
Currently there are no memory leaks. However, it would be nice for further 
implementation if it could be made so that the user doesn't have to free every 
call to getDeriv.
////
Investigate why if condition is now neccessary in getRemainder function(where 
rTerms will equal zero), for cases such as 2x + 3. It didn't used to be, to my 
knowledge.


Thinks that x^5 has a root, doesn't come up with exactly 0.

Root bounds in the Sturm paper may prove useful to pick starting points 
for Newton's method.

You may not be able to work with non-integer degrees due to synthetic 
division and (maybe) sturm's theorem.

test if deriv function is awesome.

Infinity check is now DIAMONDS (And by that I mean symbolic so BRING ON GIANT 
ROOTS)

Newton's method seems to converge perfectly to a quadratic with a root of 
-10000. Maybe it's because it's a quadratic that it doesn't have trouble.
*/
#include <iostream>
#include <cmath>
#include <string>
#include <cstdio>
#include <cstdlib> 
using namespace std;

typedef struct term {
	double degree;
	double coeff;
} term;

term* func();
/* This function is responsible for getting the equation input from the user and
   for filling in the term array properly, so that it can be used by the other
   functions */
double getVal(double x, term* func);
/* This function returns the value of the function func at x */
term* getDeriv(term* func, int n);
/* This function returns the function for the nth derivative of the function
   func using simple symbolic polynomial differentiation */
double getRoot(term* equation);
/* This function returns a root of the function equation */
term* divideOut(double a, term* equation);
/* This function divides out (x - a) from function equation using simple 
   synthetic division and returns the result */
int rootNum(term* equation);
/* This function determines how many real roots the polynomial has by utilizing
   a Sturm chain */
term* returnRemainder(term* p, term* q);
/* This function is the magnum opus. It is required for the rootNum() function
   and determines the quotient and remainder of the result of polynomial divison
   of p/q using an expanded synthetic division algorithm for non-monic divisors.
   It then returns the remainder function. */
int compareRoots(const void*, const void*);

int main()
{
	term* equation;
	int i, rNum;
	double root;
//	double x;
	srand(time(0));
	equation = func();
//	for (i = 0; i < highest; ++i)
//		cout << equation[i].coeff << ' ' << equation[i].degree << ' '
//		     << endl; //for debugging
//	cout << "Please enter the value of x: ";
//	cin >> x;
//	cout << "The value of the function at x is " << getVal(x, equation) << endl;
//	cout << "The value of the derivative of the function at x is "
//	     << getVal(x, getDeriv(equation, 1)) << endl;
	cout << "There should be " << (rNum = rootNum(equation)) << " real roots"
	     << endl;
	cout << "The roots are: ";
	for (i = 0; i < rNum; ++i) {
		cout << (root = getRoot(equation)) << ((i == rNum-1) ? "" : ", ");
		equation = divideOut(root, equation);
	}
	cout << endl;
	delete[] equation;
	return 0;
}

term* func()
{
	string cur, curCpy, operators;
	int i, n, eaten = 0, deg = 0, num, p, coeffLen = 0, highest = 0;
	double a;
	char b, c;
	char aStr[128], degStr[128];
	term *tempEquation;
	term *equation;
	signed short int *negative;
	
	cout << "Enter function:" << endl;
	getline(cin, cur);
	
	for (n = 0; n < cur.length(); ++n) { 
		while (cur[n] == ' ')
			cur.erase(n, 1);
	} // Super duper space remover
	for (n = 0; n < cur.length(); ++n) { 
		while (cur[n] == '*')
			cur.erase(n, 1);
	} //Super duper asterisk remover
	// Makes it appear like it understands things like 2*x^2 + 5*x but it really
	// just removes the asterisks. Still will not understand 2*3x properly (it
	// will interpret it as 23x)
	curCpy = cur;
	cur = "";
	for (i = 0;; ++i) {
		if (curCpy[0] == '+' || curCpy[0] == '-') {
			operators.push_back(curCpy[0]);
			curCpy.erase(0, 1);
		}
		p = sscanf(curCpy.c_str(),"%lf%n%[x]%[ ^]%i%n", &a, &coeffLen, &b, 
	                                                    &c, &deg, &eaten);
		if (p == 0) {
			curCpy.insert(0, "1");
			p = sscanf(curCpy.c_str(),"%lf%n%[x]%[ ^]%i%n", &a, &coeffLen, &b, 
		                                                    &c, &deg, &eaten);
		}
		if (p == 1) {
			curCpy.insert(coeffLen, "x^0");
			p = sscanf(curCpy.c_str(),"%lf%n%[x]%[ ^]%i%n", &a, &coeffLen, &b, 
		                                                    &c, &deg, &eaten);
		}
		if (p == 2) {
			curCpy.insert(coeffLen + 1, "^1");
			p = sscanf(curCpy.c_str(),"%lf%n%[x]%[ ^]%i%n", &a, &coeffLen, &b, 
		                                                    &c, &deg, &eaten);
		}
		sprintf(aStr, "%.30g", a);
		sprintf(degStr, "%i", deg);
		cur = cur + aStr + "x" + "^" + degStr + " ";
		curCpy.erase(0,eaten);
		if (deg > highest)
			highest = deg;
		if (curCpy == "") {
			cur.erase(cur.length() - 1, cur.length());
			num = i + 1;
			break;
		}
	} /* This for loop is used to determine the number of terms and to make the 
	     input understandable (as in the case of '2x') before actually processing 
	     it into an array */
	++highest;
	equation = new term[highest];
	tempEquation = new term[num];
	for (i = 0; i < num; ++i) {
		sscanf(cur.c_str(),"%lf%*c%*c%lf%n", &(tempEquation[i].coeff), 
			   &(tempEquation[i].degree), &eaten);
		cur.erase(0,eaten);
		// The above now works magically with equations like x^2 + 3x - 5 
		// because of the shennanigans performed above, but now something like 
		// x^2 - -3x - 5 no longer works and in fact breaks it some times for a 
		// as of yet unknown reason.
		// Also I'm not 100% sure that the erasing is necessary
	}
	negative = new signed short int[num];
	negative[0] = 1;
	for((operators.length() == num ? i = 0 : i = 1), n = 0; i < num;
	    ++i, ++n ) {
		if (operators[n] == '-')
			negative[i] = -1;
		else
			negative[i] = 1;
	}
	for ((operators.length() == num ? i = 0 : i = 1); i < num; ++i) {
		tempEquation[i].coeff *= negative[i];
	}
	for (i = 0; i < highest; ++i) {
		equation[i].degree = (highest-1)-i;
		equation[i].coeff = 0;
	}
	for (i = 0; i < highest; ++i) {
		for (n = 0; n < num; ++n) {
			if (equation[i].degree == tempEquation[n].degree) {
				equation[i].coeff = tempEquation[n].coeff;
				// Should be able to add break here but last time I think it 
				// broke it.
			}
		}
	}
	delete[] tempEquation;
	delete[] negative;
	return equation;
}

double getVal(double x, term* a)
{
	int i;
	double val = 0;
	for (i = 0; i < a[0].degree + 1; ++i)
		val += a[i].coeff*pow(x,a[i].degree);
	return val;
}

term* getDeriv(term *a, int order)
{
	term *deriv;
	term *derivClean;
	int limit = a[0].degree + 1;
	int i, n, q = 0;
	deriv = new term[limit];
	for (i = 0; i < limit; ++i)
		deriv[i] = a[i];	
	for (i = 0; i < order; ++i) {
		for (n = 0; n < limit; ++n) {
			deriv[n].coeff *= deriv[n].degree;
			--deriv[n].degree;
			if (deriv[n].degree == -1)
				deriv[n].degree = 0;
		}
	}
	for (i = 0; i < limit; ++i)
		if (deriv[i].degree == 0)
			++q;
	derivClean = new term[limit - (q-1)];
	for (i = 0; i < limit - (q-1); ++i) {
		derivClean[i].coeff = deriv[i].coeff;
		derivClean[i].degree = deriv[i].degree;
	}
	delete[] deriv;
	return derivClean;
}

double getRoot(term *equation)
{
	double x = 1; // need to develop a way to choose a better starting point.
	int i;
	if (getVal(0, equation) == 0)
		return 0;
	term* deriv = getDeriv(equation, 1);
	// isfinite only works up to long double. In case you ever want to make more 
	// precise. I'd imagine signbit() has a similar requirement
	while (!isfinite(getVal(x, equation)/getVal(x, deriv)))
		x += (static_cast<double>(rand())/RAND_MAX);
	for (i = 0; i < 99; ++i) {
		x = x - (getVal(x, equation)/getVal(x, deriv));
	}
	delete[] deriv;
	return x;
}

term* divideOut(double x, term* equation)
{
	int i, highest = equation[0].degree + 1;
	term* cleanEq;
	for (i = 1; i < highest; i++)
		equation[i].coeff = x * equation[i - 1].coeff + equation[i].coeff; 
	for (i = 0; i < highest; i++) {
		--equation[i].degree;
	}
	if (equation[highest-1].degree == -1) {
		cleanEq = new term[highest-1];
		for(i = 0; i < highest - 1; ++i)
			cleanEq[i] = equation[i];
	}
	delete[] equation;
	return cleanEq;
}

int rootNum(term* equation)
{
	int i, n, c = 0, j = 0, k = 0, limit;
	limit = equation[0].degree + 2;
	term *f[limit];
	f[0] = equation;
	f[1] = getDeriv(f[0], 1);
	for(i = 2; i < equation[0].degree + 2; ++i) {
		f[i] = returnRemainder(f[i-2],f[i-1]);
		for (n = 0; n < f[i][0].degree + 1; ++n)
			f[i][n].coeff = -f[i][n].coeff;		
		if (f[i][0].degree == 0){
			if (f[i][0].coeff == 0) {
				delete[] f[i];
				c = i;
				break;
			}
			c = i + 1;
			break;
		}
	}
	// c tells number of equations in f	
	bool sign[c];
	
	// Positive Infinity test
	for (i = 0; i < c; ++i) {
		sign[i] = signbit(f[i][0].coeff);
	}
	for (i = 1; i < c; ++i) {
		if (sign[i-1] != sign[i])
			++j;
	}
	
	// Negative Infinity Test
	for (i = 0; i < c; ++i) {
		if (static_cast<int>(f[i][0].degree) % 2 == 0)
			sign[i] = signbit(f[i][0].coeff);
		else if (static_cast<int>(f[i][0].degree) % 2 == 1)
			sign[i] = signbit(-(f[i][0].coeff));
		else
			cout << "Something asploded." << endl;
	}
	for (i = 1; i < c; ++i) {
		if (sign[i-1] != sign[i])
			++k;
	}
	
	
	for(i = 1; i < c; ++i)
		delete[] f[i];

	return (abs(j - k));
}

term* returnRemainder(term* a, term* b)
{
	// make certain sent equation is fully fleshed out, zero coeff and all
	int highestDegreeA, highestDegreeB, i, n, qTerms;
	highestDegreeA = a[0].degree + 1;
	highestDegreeB = b[0].degree + 1;
	qTerms = (highestDegreeA - (b[0].degree));
	double coeffB[highestDegreeB-1], coeffA[highestDegreeA];
	double leadCoeff = b[0].coeff;
	term *quotient = new term[qTerms]; 
	// If the dividend is of any degree, the bar separating quotient and 
	// remainder goes at dividend degree + 1 - divisor degree

	// Store the negative of the coefficients of b into coeffB, not including 
	// the first
	for (i = 1; i < highestDegreeB; ++i) {
		coeffB[i - 1] = -(b[i].coeff);
	}
	for (i = 0; i < highestDegreeA; ++i) {
		coeffA[i] = a[i].coeff;
	}
	for (i = qTerms - 1, n = 0; i >= 0; --i, ++n)
		quotient[i].degree = n;
	quotient[0].coeff = coeffA[0]/leadCoeff;
	for (i = 1; i < qTerms; ++i) {
		quotient[i].coeff = coeffA[i];
		for (n = 0; n < i; ++n) {
			quotient[i].coeff += quotient[i-1].coeff*coeffB[n];
			// I think valgrind believes the above illegally tries to read past 
			// the array in the case of 5x^5 + 5x^1 + 5x^0 as the equation. 
			// Testing is required.
		}
		quotient[i].coeff /= leadCoeff;
	}
	if (b[0].degree != 0) { //Ugliest workaround ever. This used to work though, 
	// without having to do this. Why does it break now without it?
	term *remainder;
	int rTerms = b[0].degree;
	remainder = new term[rTerms];
	int q, z, m;
	for (i = rTerms - 1, n = 0, z = 0; i >= 0; --i, ++n)
		remainder[i].degree = n;
	// start loop at qTerms for coeffA
	for (i = qTerms,  z = 0, q = rTerms; i < highestDegreeA; ++i, ++z, --q) {
		remainder[z].coeff = coeffA[i];
		for (n = 0, m = 1; n < q && m > 0; ++n) {
			m = qTerms - (n+1);
			remainder[z].coeff += quotient[m].coeff*coeffB[n+z];
		}
	}
	// Clean out leading near-zero entries. Not certain about what bounds should 
	// be.
	// Maybe cycle through and replace near-zero with zero and then just check 
	// for zero instead? Maybe.
	for (i = 0; remainder[i].coeff < 1e-100 && remainder[i].coeff > -1e-100 && 
	     i < rTerms - 1; ++i); // I think this is the issue. It's trying to 
	     // to read something it ain't sure even exists.
	     // Except even then, what does it return? What is in remainder, if 
	     // anything?
	term *cleanRemainder = new term[rTerms - i];
	for (n = 0; n < (rTerms - i); ++n)
		cleanRemainder[n] = remainder[n+i];
	delete[] remainder;
	return cleanRemainder;
	}
	else {
		term *zero = new term[1];
		zero[0].coeff = 0;
		zero[0].degree = 0;
		delete[] quotient;
		return zero;
	}
}

int compareRoots(const void *a, const void *b)
{
	// This will sort roots. It can sort in ascending order right now, but I am 
	// unable to determine an exact way to put the negative valued root first, 
	// when they exist in a pair, due to them not being exactly equal to the -.
	return ceil((abs(*(double*)a) - abs(*(double*)b)));
}
