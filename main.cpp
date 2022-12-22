// Written by: Pashshoev Bakhtierzhon KMBO-04-19
// Contains functions to solve ______ equations:
    // QUADRATIC ---> quadraticEqSolve(a, b, c, &roots)
    // CUBIC -------> cubicEqSolve(a, b, c, d, &roots)
    // QUARTIC -----> quarticEqSolve(a, b, c, d, e, &roots)

#include <iostream>
#include <cmath>
#include <vector>
#include "excerptMod2.h"
#include <complex>

#define PRINT true // for printing results of tests

using namespace std;

//typedef long double fp_t;
//typedef double fp_t;
typedef float fp_t;

//______________________            HELPER FUNCTIONS        ______________________

// Flexible comparison of 2 floating point variables
template<typename fp_t>
inline bool isEqual(const fp_t  &A, const fp_t  &B)
{
    return abs(A - B) < numeric_limits<fp_t>::epsilon() * (abs(B) + abs(A)) * 0.5;
}

// Flexible comparison to zero
template<typename  fp_t>
inline bool isZero( const fp_t &x )
{ return FP_ZERO == std::fpclassify(x); }

// Suppression cubic member of quartic
// x^4+b*x^3+c*x^2+dx+e  ---> y^4 + px^2 + qx + r, x = y - b/4
template<typename fp_t>
inline void preProcessing(fp_t  b, fp_t  c, fp_t d, fp_t  e, fp_t & p, fp_t &q,fp_t &r ){
    // p = -3/8 * b^2 + c
    p = fma<fp_t>(-b * static_cast<fp_t>(0.375L) ,b, c);
    // q = b^3 / 8 - b*c/2 + d = b/2 (b^2 * 0.25 - c) + d
    q = fma<fp_t>(static_cast<fp_t>(0.5L) * b,
                  fma<fp_t>(static_cast<fp_t>(0.25L) * b, b, -c),d);
    // r = -3 * pow(b, 4)/256 + b*b*c / 16 - b*d / 4 + e;
    r = fma<fp_t>(static_cast<fp_t>(0.25L) * b, fma<fp_t>(static_cast<fp_t>(0.25L) * b,
                                                          fma<fp_t>(static_cast<fp_t>(-0.1875L)*b, b, c),-d), e);
}

// Flexible suppression imaginary part of complex number
template<typename fp_t>
inline complex<fp_t> epsilonComplex(const complex<fp_t> &x)
{
    return abs(x) * numeric_limits<fp_t>::epsilon() > abs(x.imag()) ? complex<fp_t>(x.real(), 0) : x;
}

//______________________            MAIN FUNCTIONS        ______________________

// QUADRATIC: ax^2 + bx + c = 0
// roots - a container to put roots. Size must be >= 2
// returns: number of real roots
// Reference:
    // The Solutions of the Quadratic Equation Obtained by the Aid of the Trigonometry
    // Author(s): H. T. R. Aude
    // Source: National Mathematics Magazine, Vol. 13, No. 3 (Dec., 1938), pp. 118-121
    // Published by: Mathematical Association of America
    // Stable URL: http://www.jstor.org/stable/3028750
template<typename fp_t>
int quadraticEqSolve(fp_t a, fp_t b, fp_t c, vector<fp_t> &roots){
    // Normalizing
    if(isinf(b /= a))  return 0;
    if(isinf(c /= a)) return 0;
    a = 1;

    // constants
    const fp_t oneHalf=static_cast<fp_t>(0.5L);
    const fp_t pi = static_cast<fp_t>(numbers::pi);

    // temp variables
    fp_t sqrC, tang, tmp;
    complex<fp_t> alfa;

    if (c < 0){
        sqrC = sqrt(-c);
        tmp = -b/(2 * sqrC);
        if(isinf(tmp)) return 0;// check if there is no inf
        alfa = epsilonComplex<fp_t>(-atan<fp_t>(tmp));

        if (! alfa.imag()) {
            tang = tan(fma<fp_t>(pi, oneHalf, alfa.real()) * oneHalf);
            roots[0] = -sqrC * tang;
            roots[1]= sqrC / tang;

            if(isinf(roots[1])) return 1; // check if there is no inf
            return 2;
        }
        return 0; // No real roots
    }
    else if (abs(b) >= 2 * sqrt(c)){
        sqrC = sqrt(c);
        tmp = -2 * sqrC / b;
        if(isinf(tmp)) return 0; // check if there is no inf
        alfa = epsilonComplex(asin<fp_t>(tmp));
        tang = tan(alfa.real() * oneHalf);

        if(!alfa.imag()){
            roots[0] = sqrC * tang;
            roots[1] = sqrC / tang;

            if(isinf(roots[1])) return 1; // check if there is no inf
            return 2;
        }
        return 0; // No real roots
    }
    return 0; // D < 0 => No real roots
}


// CUBIC: ax^3 + bx^2 + cx + d = 0
// roots - a container to put roots. Size must be >= 3
// returns: number of real roots
// Reference:
    // Author: Not found
    // File Name: Algorithm to find the real roots of the cubic equation-computation_flow_chart.pdf
template<typename fp_t>
int cubicEqSolve(fp_t a, fp_t b, fp_t c, fp_t d, vector<fp_t> &roots){

    // Normalizing
    if(isinf(b /= a)) return 0;
    if(isinf(c /= a)) return 0;
    if(isinf(d /= a)) return 0;
    a = 1;

    fp_t e,f,g, h, i, absI, j, tmp;// temp variables that will need for computation several times (>=3)
    int numOfRoots = 0; // total number of roots
    // constants
    const fp_t oneThird = static_cast<fp_t>(1.0L/3.0L); // temp var, will use it 7 times
    const fp_t pi = static_cast<fp_t>(numbers::pi);
    const fp_t piThird = pi * oneThird;

    // temp computations
    e = b * oneThird;
    f = fma<fp_t>(-b,e,c);
    g = fma<fp_t>(-2,pow<fp_t>(e, 3),fma<fp_t>(c,e, -d));
    h = sqrt(abs(f) * 4 * oneThird);

    if (isZero<fp_t>(f)){// CASE 0: Triple root
        roots[0] = roots[1] = roots[2] = - e;
        return 3;
    }

    i = 4 * g / pow<fp_t>(h, 3);
    absI = abs(i);
    if (isinf(i)) return 0; // checking for inf

    //CASE 1: Only one real root
    if (f > 0){
        roots[0] = fma<fp_t>(h,sinh(asinh(i) * oneThird),-e);
        numOfRoots = 1;
    }
    else if (absI > 1){
        tmp = h * i/absI;
        if(isinf(tmp)) return 0;
        roots[0] = fma<fp_t>(tmp,cosh(acosh(absI) * oneThird),-e);
        numOfRoots = 1;
    }
    else // CASE 2: Three real roots
    {
        j = acos(i) * oneThird;
        roots[0] = fma<fp_t>(h, cos(j), -e);
        roots[1] = fma<fp_t>(h, cos(fma<fp_t>(2, piThird, j)), -e);
        if (isZero(j))
            roots[2] = roots[1];
        else if(isEqual(j, piThird))
            roots[2] = roots[0];
        else
            roots[2] = fma<fp_t>(h, cos(fma<fp_t>(2, piThird, -j)), -e);
        numOfRoots = 3;
    }
    return numOfRoots;
}


// QUARTIC: ax^4 + b*x^3 + cx^2 + dx + e = 0
// roots - a container to put roots. Size must be >= 4
// returns: number of real roots
// Reference:
    // On the Solution of the Real Quartic
    // Author(s): William F. Carpenter,
    // Source: Mathematics Magazine, Vol. 39, No. 1 (Jan., 1966), pp. 28-30
    // Published by: Mathematical Association of America
    // Stable URL: http://www.jstor.org/stable/2688990
template<typename fp_t>
int quarticEqSolve(fp_t a, fp_t b, fp_t c, fp_t d, fp_t e, vector<fp_t> &roots){

    // Normalizing
    if(isinf(b /= a))  return 0;
    if(isinf(c /= a)) return 0;
    if(isinf(d /= a)) return 0;
    if(isinf(e /= a)) return 0;
    a = 1;
    int numOfRoots = 0; // total number of found roots

    // x^4+b*x^3+c*x^2+dx+e  ---> y^4 + cx^2 + dx + e, x = y - b/4
    preProcessing(b,c,d,e, c, d, e);

    // getting roots from specific cubic equation
    vector<fp_t> cubicRoots(3);
    const int cnumOfRoots = cubicEqSolve<fp_t>(a,2*c,pr_product_difference<fp_t>(c,c,4, e),fma<fp_t>(-d, d, 0), cubicRoots);

    // getting a positive root from cubic solution
    fp_t r;
    const fp_t zero = static_cast<fp_t>(0.0L);
    const fp_t oneHalf = static_cast<fp_t>(0.5L), minusOneFourth = static_cast<fp_t>(-0.25L);

    if(cnumOfRoots && cubicRoots [0] > zero) r = cubicRoots[0];
    else if(cnumOfRoots > 1 && cubicRoots[1] > zero) r = cubicRoots[1];
    else if(cnumOfRoots > 1 && cubicRoots[2] > zero) r = cubicRoots[2];
    else return 0;

    // variables
    const fp_t r1 = sqrt(r);
    const fp_t tmp1 = pr_product_difference<fp_t>(r, minusOneFourth, c, oneHalf);
    const fp_t tmp2 = 1/(2 * r1);

    if(isinf(tmp2)) return 0;

    const complex<fp_t> t1 = epsilonComplex<fp_t>(sqrt<fp_t>(fma<fp_t>(-d, tmp2,tmp1)));
    const complex<fp_t> t2 = epsilonComplex<fp_t>(sqrt<fp_t>(fma<fp_t>(d, tmp2,tmp1)));

    if (isZero(t1.imag())){
        roots[0] = fma<fp_t>( b , minusOneFourth,fma<fp_t>(r1,oneHalf,t1.real()));
        roots[1] = fma<fp_t>( b , minusOneFourth,fma<fp_t>(r1,oneHalf,-t1.real()));
        numOfRoots +=2;
    }
    if (isZero(t2.imag())){
        roots[numOfRoots] = fma<fp_t>( b , minusOneFourth,fma<fp_t>(-r1,oneHalf,t2.real()));
        roots[numOfRoots + 1] = fma<fp_t>( b , minusOneFourth,fma<fp_t>(-r1,oneHalf,-t2.real()));
        numOfRoots +=2;
    }
    return numOfRoots;
}



//______________________            TESTER FUNCTIONS        ______________________


// Function to test quadratic solution:
// testCount - total count of tests
// dist - maximum distance between roots
template<typename fp_t>
void testQuadraticAdv(const int testCount, const fp_t dist){
    int const P=2; // power
    fp_t const low=-1, high=1; // [low, high]
    fp_t absMaxError, relMaxError; // variables for each test Errors
    int numOfFoundRoots, cantFind = 0;

    fp_t maxAbsAllofTest = -1, maxRelAllofTest = -1; // maximum from maxAbsoluteError and maxRelError from all [testCount] tests

    long double absErrors = 0;
    long double relError = 0;

    vector<fp_t> coefficients(P+1);
    vector<fp_t> trueRoots(P);

    for(int i=0; i < testCount; ++i){
        vector<fp_t> foundRoots(P);
        generate_polynomial<fp_t>(P, 0, 2, 0, dist,
                                  low, high, trueRoots, coefficients);

        numOfFoundRoots = quadraticEqSolve<fp_t>(coefficients[2],coefficients[1],coefficients[0], foundRoots);

        if(!numOfFoundRoots){
            cantFind++;
        }
        else{
            compare_roots<fp_t>(numOfFoundRoots, 2, foundRoots, trueRoots, absMaxError, relMaxError);
            maxAbsAllofTest = absMaxError > maxAbsAllofTest ? absMaxError : maxAbsAllofTest;
            absErrors += absMaxError;
            maxRelAllofTest = relMaxError > maxRelAllofTest ? relMaxError : maxRelAllofTest;
            relError += relMaxError;
        }
    }
    if(PRINT){
        cout<<"\n\n\t\t\tQUADRATIC TEST RESULTS\n\n";
        cout<<"Max distance: "<< dist << endl;
        cout<<"Total count of tests: "<<testCount<<endl;
        cout<<"Couldn't find roots: " << cantFind <<" times "<<endl;
        cout<<"Mean absMaxError = "<< absErrors / (testCount - cantFind) << endl;
        cout<<"Max {absMaxError_i | i = 0, ..., 1e6} from all of the tests: "<<maxAbsAllofTest<<endl;
        cout<<"Mean RelMaxError = "<< relError / (testCount - cantFind)  << endl;
        cout<<"Max {RelMaxError_i | i = 0, ..., 1e6} all of the tests: "<<maxRelAllofTest<<endl;
    }
}

// Function to test cubic solution:
// testCount - total count of tests
// dist - maximum distance between roots
template <typename fp_t>
void testCubicAdv(const int testCount, const fp_t dist){

    int const P = 3; // power, total number of tests
    int numOfFoundRoots, cantFind = 0, tripleRoot=0;

    fp_t const low=-1, high=1; // [low, high], max distance between clustered roots
    fp_t absMaxError, relMaxError; // variables for each test Errors
    fp_t maxAbsAllofTest = -1, maxRelAllofTest = -1; // maximum from maxAbsoluteError and maxRelError from all [testCount] tests

    long double absErrors = 0, relError;

    vector<fp_t> coefficients(P+1);
    vector<fp_t> trueRoots(P);

    for(int i=0; i < testCount; ++i){
        vector<fp_t> foundRoots(P);

        generate_polynomial<fp_t>(P, 0, P, 0, dist,
                                  low, high, trueRoots, coefficients);
        numOfFoundRoots = cubicEqSolve<fp_t> (coefficients[3],coefficients[2],coefficients[1],coefficients[0], foundRoots);

        if(numOfFoundRoots == 1){
            cantFind++;
        }
        else{
            compare_roots<fp_t>(numOfFoundRoots, 3, foundRoots, trueRoots, absMaxError, relMaxError);
            maxAbsAllofTest = absMaxError > maxAbsAllofTest ? absMaxError : maxAbsAllofTest;
            absErrors += absMaxError;
            maxRelAllofTest = relMaxError > maxRelAllofTest ? relMaxError : maxRelAllofTest;
            relError += relMaxError;
        }
    }

    if (PRINT) {
        cout << "\n\n\t\t\tCUBIC TEST RESULTS\n\n";
        cout << "Max distance: " << dist << endl;
        cout << "Total count of tests: " << testCount << endl;
        cout << "ONLY ONE ROOT: " << cantFind << " times " << endl;
        cout << "Mean absMaxError = " << absErrors / (testCount - cantFind) << endl;
        cout << "Max {absMaxError_i | i = 0, ..., 1e6} from all of the tests: " << maxAbsAllofTest << endl;
        cout << "Mean RelMaxError = " << relError / (testCount - cantFind) << endl;
        cout << "Max {RelMaxError_i | i = 0, ..., 1e6} all of the tests: " << maxRelAllofTest << endl;
    }
    ////////////////////////////////////////
}

// Function to test quartic solution:
// testCount - total count of tests
// dist - maximum distance between roots
template<typename fp_t>
void testQuarticAdv(const int testCount, const fp_t dist){

    const int P = 4; // power
    const fp_t low=-1, high=1; // [low, high]
    fp_t absMaxError, relMaxError; // variables for each test Errors
    int numOfFoundRoots, cantFind = 0, twoRoots = 0;
    fp_t maxAbsAllofTest = -1, maxRelAllofTest = -1; // maximum from maxAbsoluteError and maxRelError from all [testCount] tests

    long double absErrors = 0;
    long double relError = 0;

    vector<fp_t> coefficients(P+1);
    vector<fp_t> trueRoots(P);

    for(int i=0; i < testCount; ++i){
        vector<fp_t> foundRoots(P);
        generate_polynomial<fp_t>(P, 0, P, 0, dist,
                                  low, high, trueRoots, coefficients);
        numOfFoundRoots = quarticEqSolve<fp_t> (coefficients[4], coefficients[3],coefficients[2],coefficients[1],coefficients[0], foundRoots);

        if(!numOfFoundRoots){
            cantFind++;
        }
        else{
            compare_roots<fp_t>(numOfFoundRoots, 4, foundRoots, trueRoots, absMaxError, relMaxError);
            maxAbsAllofTest = absMaxError > maxAbsAllofTest ? absMaxError : maxAbsAllofTest;
            absErrors += absMaxError;
            maxRelAllofTest = relMaxError > maxRelAllofTest ? relMaxError : maxRelAllofTest;
            relError += relMaxError;
        }
    }
    if(PRINT){
        cout<<"\n\n\t\t\tQUARTIC TEST RESULTS\n\n";
        cout<<"Max distance: "<< dist << endl;
        cout<<"Total count of tests: "<<testCount<<endl;
        cout<<"Couldn't find roots: " << cantFind <<" times "<<endl;
        cout<<"Mean absMaxError = "<< absErrors / (testCount - cantFind) << endl;
        cout<<"Max {absMaxError_i | i = 0, ..., 1e6} from all of the tests: "<<maxAbsAllofTest<<endl;
        cout<<"Mean RelMaxError = "<< relError / (testCount - cantFind)  << endl;
        cout<<"Max {RelMaxError_i | i = 0, ..., 1e6} all of the tests: "<<maxRelAllofTest<<endl;
    }
}



int main(){
    setlocale(LC_ALL, "ru");
    cout<<setprecision(12);
    const int testCount = 1000'000; // total number of tests
    const fp_t dist = 1e-5;  // maximum distance between roots


    testQuadraticAdv<fp_t>(testCount, dist);
    testCubicAdv<fp_t>(testCount, dist);
    testQuarticAdv<fp_t>(testCount, dist);
    return 0;
}
