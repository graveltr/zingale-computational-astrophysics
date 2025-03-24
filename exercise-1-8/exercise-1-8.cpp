#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

constexpr double PI = 3.14159265358979323846;

/*
 * Consider the orbit of Earth around the Sun. If ew work in the units of
 * astronomical units, years, and solar masses, then Newton's gravitational
 * constant and the solar mass together are simply GM = 4 \pi^2. We can write
 * the ODE system describing the motion of Earth as
 *
 * \dot{\vec{x}} = \vec{v}
 * \dot{\vec{v}} = - G M \vec{r} / r^3
 *
 * If we take the coordinate system such that the Sun is at the origin, then
 * \vec{x} = (x, y) is the position of the Earth and \vec{r} = x \hat{x} + y
 * \hat{y} is the radius vector.
 *
 * Take as the initial conditions the planet at perihelion:
 *
 * x_0 = 0, y_0 = a (1 - e)
 * (\vec{v} \cdot \hat{x})_0 = - \sqrt{(G M / a) (1 + e / 1 - e)}
 * (\vec{v} \cdot \hat{y})_0 = 0
 *
 * where a is the semi-major axis and e is the eccentricity of the orbit.
 * Integrate this system for a single orbital period with the first-order Euler
 * and RK4 method and measure the convergence by integrating at a number of
 * different \Delta t's. You'll need to define some measure of error, you can
 * consider a number of different metrics, e.g., the change in radius after a
 * single orbit.
 */

void dumpMatrixToFile(std::vector<std::vector<double>>& matrix) {
  std::string filename = "output.txt";
  std::ofstream outFile(filename);
  
  if (!outFile) {
      std::cerr << "Error opening file for writing: " << filename << std::endl;
      return;
  }
  
  for (const auto& row : matrix) {
      for (const auto& element : row) {
          outFile << element << " ";
      }
      outFile << std::endl;  // Newline after each row
  }
  
  outFile.close();
  
  if (!outFile) {
      std::cerr << "Error occurred while writing to file: " << filename << std::endl;
  }
}

std::vector<double> computeVecYNPlusOneEuler(std::vector<double>& vecYN, double deltaT) {
  std::vector<double> vecYNPlusOne(vecYN.size());
  
  double xN = vecYN[0];
  double yN = vecYN[1];
  double vxN = vecYN[2];
  double vyN = vecYN[3];

  double r = std::sqrt(xN * xN + yN * yN);

  double xNPlus1 = xN + vxN * deltaT;
  double yNPlus1 = yN + vyN * deltaT;
  double vxNPlus1 = vxN - deltaT * (4.0 * PI * PI * xN) / (r * r * r);
  double vyNPlus1 = vyN - deltaT * (4.0 * PI * PI * yN) / (r * r * r);

  vecYNPlusOne[0] = xNPlus1;
  vecYNPlusOne[1] = yNPlus1;
  vecYNPlusOne[2] = vxNPlus1;
  vecYNPlusOne[3] = vyNPlus1;

  return vecYNPlusOne;
}

std::vector<double> computeVecYNPlusOneRK(std::vector<double>& vecYN, double deltaT) {
  std::vector<double> vecYNPlusOne(vecYN.size());
  
  double xN = vecYN[0];
  double yN = vecYN[1];
  double vxN = vecYN[2];
  double vyN = vecYN[3];

  double r = std::sqrt(xN * xN + yN * yN);

  double xNPlus1 = xN + vxN * deltaT;
  double yNPlus1 = yN + vyN * deltaT;
  double vxNPlus1 = vxN - deltaT * (4.0 * PI * PI * xN) / (r * r * r);
  double vyNPlus1 = vyN - deltaT * (4.0 * PI * PI * yN) / (r * r * r);

  vecYNPlusOne[0] = xNPlus1;
  vecYNPlusOne[1] = yNPlus1;
  vecYNPlusOne[2] = vxNPlus1;
  vecYNPlusOne[3] = vyNPlus1;

  return vecYNPlusOne;
}

int main() {
  std::vector<std::vector<double>> vecY;

  double deltaT = 0.0001;
  double currT = 0.0;
  double tMax = 50.0;

  double a = 1.0;
  double e = 0.0;

  // Initial conditions
  double x0 = 0.0;
  double y0 = a * (1.0 - e);
  double vx0 = -1.0 * std::sqrt((1.0 / a) * (4.0 * PI * PI) * ((1.0 + e) / (1.0 - e)));
  double vy0 = 0.0;
  vecY.push_back({x0, y0, vx0, vy0});
  while (currT < tMax) {
    vecY.push_back(computeVecYNPlusOneEuler(vecY.back(), deltaT));
    currT += deltaT;
    std::cout << currT << std::endl;
  }

  dumpMatrixToFile(vecY);
  return 0;
}

/*
  for (const auto& row: vecY) {
    for (const auto& elt: row) {
      std::cout << elt << " ";
    }
    std::cout << std::endl;
  }
*/