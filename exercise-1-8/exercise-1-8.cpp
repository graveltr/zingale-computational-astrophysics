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

// Overload operator* to scale a vector by a double (scalar)
std::vector<double> operator*(const std::vector<double>& vec, double scalar) {
    std::vector<double> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = vec[i] * scalar;
    }
    return result;
}

// Overload operator* to scale a vector by a double (scalar) - other direction
std::vector<double> operator*(double scalar, const std::vector<double>& vec) {
    return vec * scalar;
}

// Function to add two vectors element-wise
std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vectors must be of the same size");
    }
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] + b[i];
    }
    return result;
}

void dumpMatrixToFile(std::vector<std::vector<double>>& matrix, std::string filename) {
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

std::vector<double> evaluateF(const std::vector<double>& vec) {
  double x = vec[0];
  double y = vec[1];
  double vx = vec[2];
  double vy = vec[3];

  double r = std::sqrt(x * x + y * y);

  std::vector<double> ret({
    vx,
    vy, 
    -1.0 * (4.0 * PI * PI * x) / (r * r * r),
    -1.0 * (4.0 * PI * PI * y) / (r * r * r)
  });

  return ret;
}

std::vector<double> computeVecYNPlusOneEuler(std::vector<double>& vecYN, double deltaT) {
  return vecYN + deltaT * evaluateF(vecYN);
}

std::vector<double> computeVecYNPlusOneRK(std::vector<double>& vecYN, double deltaT) {
  std::vector<double> k1, k2, k3, k4;

  k1 = evaluateF(vecYN);
  k2 = evaluateF(vecYN + ((deltaT / 2.0) * k1));
  k3 = evaluateF(vecYN + ((deltaT / 2.0) * k2));
  k4 = evaluateF(vecYN + deltaT * k3);

  return vecYN + (deltaT / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
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
  dumpMatrixToFile(vecY, "euler-output.txt");

  currT = 0.0;
  vecY.clear();
  vecY.push_back({x0, y0, vx0, vy0});
  while (currT < tMax) {
    vecY.push_back(computeVecYNPlusOneRK(vecY.back(), deltaT));
    currT += deltaT;
    std::cout << currT << std::endl;
  }
  dumpMatrixToFile(vecY, "rk-output.txt");

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