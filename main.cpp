#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

//typedef double fp_t;
typedef float fp_t;

template <typename fp_t>
int sgn(fp_t val) {
    return (fp_t(0) < val) - (val < fp_t(0));
}

template <typename fp_t>
int discriminant(fp_t a, fp_t b, fp_t c, vector<fp_t> &roots){
    //a, b, c - coefficients of equation
    //roots - a container to put roots. Size must be >= 2

    //returns: number of real roots
    if (!a)
        return 0;
    b/=a;
    c/=a;
    a/=a;
    fp_t D= fma(b,b, -4*a*c);
    if(D < 0)
        return 0;
    fp_t chisl=abs(b)+ sqrt(D);
    fp_t znam=-2*a*sgn<fp_t>(b);
    roots[0] =chisl/znam;
    if(!roots[0])
        roots[1]=c/roots[0];
    else
        roots[1]=roots[0];
    return 2;
}

template<typename fp_t>
int completingTheSquare(fp_t a, fp_t b, fp_t c, vector<fp_t> &roots){
    //a, b, c - coefficients of equation
    //roots - a container to put roots. Size must be >= 2

    //returns: number of real roots
    if (!a)
        return 0;
    b/=a;
    c/=a;
    a/=a;
    fp_t x1,x2, h=-b/2, k=fma(h,h, -c); ///////////
    if (k < 0)
        return 0;
    roots[0]=h+sqrt(k);
    roots[1]=h-sqrt(k);

}

template<typename fp_t>
int quadTrigAude(fp_t a, fp_t b, fp_t c, vector<fp_t> &roots){
    // Trigonometric solution for quadratic equation
    //a, b, c - coefficients of equation
    //roots - a container to put roots. Size must be >= 2

    //returns: number of real roots
    if (!a || !b && c>0){
        return 0;
    }

    fp_t x1,x2, alfa, pi = static_cast<fp_t>(numbers::pi);
    b/=a; c/=a; a/=a;

    if (!c){
        roots[0]=0;
        roots[1]=-b;
        return 2;
    }

    if (c < 0){
        alfa=pi/2 - atan(-b/(2 * sqrt(-c)));
        x1=sqrt(-c/a)*1/tan(alfa/2);
        x2=-sqrt(-c/a)*tan(alfa/2);
    }
    else if (b >= 2*sqrt(c)){
        alfa=std::asin(-2*sqrt(c)/b);
        x1=sqrt(c/a)*1/tan(alfa/2);
        x2=sqrt(c/a)*tan(alfa/2);
    }
    else{
        return 0; // D < 0 => No real roots
    }

    roots[0]=x1;
    roots[1]=x2;
    return 2;
}

template<typename fp_t>
int cubicEqSolution(fp_t a, fp_t b, fp_t c, fp_t d, vector<fp_t> &roots){
    //a, b, c, d - coefficients of equation
    //roots - a container to put roots. Size must be >= 3

    //returns: number of real roots
    if (!a) // it's not a cubic
        return 0;

    fp_t e,f,g, tmp, h, i, j;// var's that will need for computation
    fp_t x1, x2, x3; // roots
    fp_t pi=static_cast<fp_t>(numbers::pi);

    e=b/(3*a);
    f=fma<fp_t>(-b,e,c)/a;
    tmp=fma<fp_t>(c,e, -d);
    g=fma<fp_t>(-2*a,pow<fp_t>(e,3),tmp)/a;

    // triple root
    fp_t oneThird=static_cast<fp_t>(1.0L/3);
    if (!f){
        x1=fma<fp_t>(pow<fp_t>(g,oneThird),1,-e);
        cout<<"Triple root\n";
        roots[0]=x1;
        roots[1]=x1;
        roots[2]=x1;
        return 3;
    }
    h = sqrt(4*abs(f)/3);
    i = 4*g/pow<fp_t>(h,3);

    //one real root
    if (f>0){
        x1=fma<fp_t>(h,sinh(asinh(i)/3),-e);
        roots[0]=x1;

        cout<<"Only one root1\n";
        return 1;
    }
    fp_t absI=abs(i);
    if (absI>1){
        x1=fma<fp_t>(h*i/absI,cosh(acosh(absI)/3),-e);
        roots[0]=x1;

        cout<<"Only one root2\n";
        return 1;
    }

    // 3 различные корни или существует корень кратности 2
    j=acos(i)/3;
    x1=fma<fp_t>(h, cos(j), -e);
    x2=fma<fp_t>(h, cos(2*pi/3+j), -e);
    roots[0]=x1;
    roots[1]=x2;
    if (j==0){
        cout<<"j="<<j<<endl;
        roots[2]=x2;
        cout<<"2 unique roots\t"<<endl;
    }
    else if(j==pi/3){
        cout<<"j="<<j<<endl;
        roots[2]=x1;
        cout<<"2 unique roots\t"<<endl;
    }
    else{
        roots[2]=fma<fp_t>(h, cos(2*pi/3-j), -e);
        cout<<"3 unique roots\t"<<endl;
    }
    return 3;
}

template<typename fp_t>
int quarticCruchagaGarver(fp_t a, fp_t b, fp_t c, fp_t d, fp_t e, vector<fp_t> &roots){
    // Method described by Garver, updated version of Cruchaga's
    // Solves equation like: ax^4 + cx^2 + dx + e (without cubic (b*x^3) member)
    //a, b, c, d, e - coefficients of equation
    //roots - a container to put roots. Size must be >= 3

    //returns: number of real roots
    if (!a || b){
        return 0;
    }
    b/=a; c/=a; d/=a; e/=a; a/=a;

    fp_t newC=fma<fp_t>(c,c,-4*e);

    fp_t k1, k2, k3;
    vector<fp_t> cubicRoots(3);
    int numOfRoots= cubicEqSolution<fp_t>(a,2*c,newC,-d*d, cubicRoots);

    if (numOfRoots==3){
        k1=sqrt(cubicRoots[0]);
        k2=sqrt(cubicRoots[1]);

        fp_t tmpRoot=sqrt(cubicRoots[2]);
        fp_t tmpProd=k1*k2*tmpRoot;
        if (tmpProd==-d){
            k3=tmpRoot;
        }else{
            k3=-tmpRoot;
        }
        roots[0]=(k1+k2+k3)/2;
        roots[1]=(k1-k2-k3)/2;
        roots[2]=(-k1+k2-k3)/2;
        roots[3]=(-k1-k2+k3)/2;
        return 4;
    }
    return 0;
}

template<typename fp_t>
int quarticCarpenter(fp_t a, fp_t b, fp_t c, fp_t d, fp_t e, vector<fp_t> &roots){
    //Method described by Carpenter
    //Solves equation like: ax^4 + cx^2 + dx + e (without cubic (b*x^3) member)
    //a, b, c, d, e - coefficients of equation
    //roots - a container to put roots. Size must be >= 3

    //returns: number of real roots
    if (!a || b){
        return 0;
    }
    b/=a; c/=a; d/=a; e/=a; a/=a;
    fp_t newC=fma<fp_t>(c,c,-4*e);

    vector<fp_t> cubicRoots(3);
    int numOfRoots = cubicEqSolution<fp_t>(a,2*c,newC,-d*d, cubicRoots);
    fp_t r;
    if(cubicRoots[0]>0) r=cubicRoots[0];
    else if(numOfRoots > 1 && cubicRoots[1]>0) r=cubicRoots[1];
    else if(numOfRoots > 1 &&  cubicRoots[2]>0) r=cubicRoots[2];
    else return 0;

    fp_t r1=sqrt(r);
    fp_t tmp1=fma<fp_t>(-r/4,1,-c/2);
    fp_t tmp2=sqrt(fma<fp_t>(tmp1, 1,-d/(2*r1)));
    roots[0]=fma<fp_t>(r1/2,1,tmp2);
    roots[1]=fma<fp_t>(r1/2,1,-tmp2);

    fp_t tmp3=sqrt(fma<fp_t>(d, 1/(2*r1),tmp1));
    roots[2]=fma<fp_t>(r1/2,1,tmp3);
    roots[3]=fma<fp_t>(r1/2,1,-tmp3);

    return 4;
}



void testQuadratic(){

    fp_t a=1.0,b=2.0,c=-4.0;
    vector<fp_t> res2(2);
    int num_of_roots = quadTrigAude<fp_t>(a,b,c,res2);

    if (num_of_roots){
        cout<<"x1 ="<<res2[0]<<" x2 = "<<res2[1];
    }
    else{
        cout<<"No real roots\n";
    }
}
void testCubic(){
    fp_t a=1.0,b=0.0,c=-4.0, d=-15.0; //  3
//    fp_t a=1.0,b=6.0,c=14.0, d=15.0; //  -3
//    fp_t a=1.0,b=103.0,c=305.0, d=500.0; //  -100
//
//    fp_t a=2.0,b=58.0,c=190.0, d=-250.0; // 1, -5, -25
//    fp_t a=-5,b=30,c=225, d=-250; // 1, -5, 10
//    fp_t a=1.0,b=54.0,c=195.0, d=-250.0; //  1, -5, -250
//
//    fp_t a=1.0,b=3.0,c=-9.0, d=5.0; // 1, 1, -5
//    fp_t a=1.0,b=9.0,c=15.0, d=-25.0; // 1, -5, -5
//    fp_t a=1.0,b=7.0,c=11.0, d=5.0; // -1, -1, -5
//
//    fp_t a=1.0,b=-9.0,c=27.0, d=-27.0; //  3, 3, 3
//    fp_t a=-125,b=750,c=-1500, d=1000; //  2, 2, 2
//    fp_t a=1,b=-18,c=108, d=-216; //  -6, -6, -6
//    fp_t a=1,b=-18.75,c=117.1875, d=-244.140625; //  -6.25, -6.25, -6.25

    vector<fp_t> res3(3);

    int num_of_roots= cubicEqSolution<fp_t>(a,b,c,d,res3);
    if (num_of_roots==3){
        cout<<"x1 = "<<res3[0]<<"; x2 = "<<res3[1]<<"; x3 = "<<res3[2]<<endl;
    }
    else if(num_of_roots){
        cout<<"x1 = "<<res3[0]<<endl;
    }
}




int main() {

//Quadratic
//    testQuadratic();

//Cubic equation
//    testCubic();
// Quartic
    return 0;
}
