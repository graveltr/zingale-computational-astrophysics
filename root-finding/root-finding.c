#include <stdio.h>
#include <math.h>
#include <float.h>

double f1(double x) {
  return x * x - 1.0;
}

double evalFunction(double (*fptr)(double), double x) {
  return fptr(x);
}

double evalDerivative(double (*func)(double), double x) {
  double deltax = 0.01;
  double rightVal = func(x + deltax);
  double leftVal = func(x - deltax);

  return (rightVal - leftVal) / (2.0 * deltax);
}

double findRoot(double (*func)(double), double x0) {
  double deltaxGoal = 0.5;
  double currDeltax = DBL_MAX;
  double currX = x0;

  int maxIter = 30; 
  int currIter = 0;

  while (deltaxGoal < fabs(currDeltax) && currIter < maxIter) {
    printf("deltaxgoal: %lf \n", deltaxGoal);
    printf("currDeltax: %lf \n", currDeltax);
    printf("is deltaxgoal less than currdeltax?: %d \n", deltaxGoal < currDeltax);

    printf("iter: %d \n", currIter);

    printf("currX: %lf \n", currX);

    double currY = func(currX);
    printf("currY: %lf \n", currY);

    double currYPrime = evalDerivative(func, currX);
    currDeltax = -1.0 * (currY / currYPrime);
    printf("currDeltax: %lf \n", currDeltax);

    currX = currX + currDeltax;
    currIter++;
  }

  return currX;
}

int main() {
  double (*f1ptr)(double) = f1;

  printf("root: %lf \n", findRoot(f1ptr, 3.0));
  
  return 0;
}
