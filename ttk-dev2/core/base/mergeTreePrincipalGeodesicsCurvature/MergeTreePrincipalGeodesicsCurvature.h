/// \ingroup base
/// \class ttk::MergeTreePrincipalGeodesicsCurvature
/// \author XXX
/// \date 2022.
///
/// This module defines the %MergeTreePrincipalGeodesicsCurvature class that
/// computes TODO
///
/// \b Related \b publication: \n
/// TODO
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <MergeTreePrincipalGeodesicsBase.h>
#include <Triangulation.h>

namespace ttk {

  /**
   * The MergeTreePrincipalGeodesicsCurvature class provides methods to compute
   * TODO
   */
  class MergeTreePrincipalGeodesicsCurvature
    : virtual public Debug,
      public MergeTreePrincipalGeodesicsBase {

  protected:
    // distances[i][j] stores the right and down distance for the node [i][j]
    // in the grid
    std::vector<std::vector<std::vector<double>>> distances_;
    // diagDistances_[i][j] stores the diagonal distance (upper-left to
    // down-right) for the node [i][j] in the grid
    std::vector<std::vector<double>> diagDistances_;

  public:
    MergeTreePrincipalGeodesicsCurvature();

    template <class dataType>
    void computeSurfaceDistances(
      std::vector<ftm::MergeTree<dataType>> &surfaceTrees,
      bool computeDiagDistances = false) {
      // Preprocessing
      preprocessingTrees<dataType>(surfaceTrees);

      // Computation
      unsigned int dim = std::sqrt(surfaceTrees.size()); // k
      distances_ = std::vector<std::vector<std::vector<double>>>(
        dim,
        std::vector<std::vector<double>>(dim, std::vector<double>(2, 0.0)));
      diagDistances_
        = std::vector<std::vector<double>>(dim, std::vector<double>(dim, 0.0));
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif
      for(unsigned int i = 0; i < dim; ++i) {
        for(unsigned int j = 0; j < dim; ++j) {
          int index = i * dim + j;
          if(j != dim - 1) {
            int index2 = i * dim + j + 1;
            dataType distance;
            computeOneDistance<dataType>(
              surfaceTrees[index], surfaceTrees[index2], distance, true);
            distances_[i][j][0] = distance;
          }
          if(i != dim - 1) {
            int index2 = (i + 1) * dim + j;
            dataType distance;
            computeOneDistance<dataType>(
              surfaceTrees[index], surfaceTrees[index2], distance, true);
            distances_[i][j][1] = distance;
          }
          if(computeDiagDistances and j != dim - 1 and i != dim - 1) {
            int index2 = (i + 1) * dim + j + 1;
            dataType distance;
            computeOneDistance<dataType>(
              surfaceTrees[index], surfaceTrees[index2], distance, true);
            diagDistances_[i][j] = distance;
          }
        }
      }

      // Postprocessing
      for(unsigned int i = 0; i < surfaceTrees.size(); ++i)
        postprocessingPipeline<dataType>(&(surfaceTrees[i].tree));
    }

    template <class dataType>
    void computeSurfaceArea(std::vector<ftm::MergeTree<dataType>> &surfaceTrees,
                            std::vector<std::vector<double>> &surfacePoints,
                            std::vector<std::vector<double>> &pointsArea,
                            std::vector<std::vector<double>> &wassersteinArea,
                            std::vector<std::vector<double>> &ratioArea) {
      computeSurfaceDistances<dataType>(surfaceTrees, true);

      unsigned int dim = std::sqrt(surfacePoints.size()); // k
      pointsArea = std::vector<std::vector<double>>(
        dim - 1, std::vector<double>(dim - 1, 0.0));
      wassersteinArea = std::vector<std::vector<double>>(
        dim - 1, std::vector<double>(dim - 1, 0.0));
      ratioArea = std::vector<std::vector<double>>(
        dim - 1, std::vector<double>(dim - 1, 0.0));
      for(unsigned int i = 0; i < dim - 1; ++i) {
        for(unsigned int j = 0; j < dim - 1; ++j) {
          int i0 = i * dim + j;
          int i1 = i * dim + j + 1;
          int i2 = (i + 1) * dim + j;
          int i3 = (i + 1) * dim + j + 1;
          pointsArea[i][j]
            = triangleArea3D(
                surfacePoints[i0], surfacePoints[i1], surfacePoints[i3])
              + triangleArea3D(
                surfacePoints[i0], surfacePoints[i2], surfacePoints[i3]);

          wassersteinArea[i][j]
            = triangleAreaFromSides(distances_[i][j][0], diagDistances_[i][j],
                                    distances_[i][j + 1][1])
              + triangleAreaFromSides(distances_[i][j][1], diagDistances_[i][j],
                                      distances_[i + 1][j][0]);

          ratioArea[i][j] = wassersteinArea[i][j] / pointsArea[i][j];
        }
      }
    }

    void computeInputPoints(std::vector<std::vector<double>> &inputTs,
                            std::vector<std::vector<double>> &surfacePoints,
                            std::vector<std::vector<double>> &surfaceTs,
                            std::vector<std::vector<double>> &inputPoints) {
      unsigned int noPoints = inputTs.size();
      inputPoints = std::vector<std::vector<double>>(
        noPoints, std::vector<double>(3, 0.0));

      // Find quadrant of each point
      // 0 --- 2
      // |     |
      // 1 --- 3
      unsigned int dim = std::sqrt(surfaceTs.size()); // k
      std::vector<std::vector<int>> quadPoints(noPoints, std::vector<int>(4));
      for(unsigned int i = 0; i < noPoints; ++i) {
        int i0 = 0, i1 = 0;
        for(unsigned int j = 0; j < dim; ++j) {
          double t = 1.0 / (dim - 1) * j;
          if(t < inputTs[i][0])
            i0 = j;
          if(t < inputTs[i][1])
            i1 = j;
        }
        quadPoints[i][0] = i0 * dim + i1 + 1;
        quadPoints[i][1] = i0 * dim + i1;
        quadPoints[i][2] = (i0 + 1) * dim + i1 + 1;
        quadPoints[i][3] = (i0 + 1) * dim + i1;
      }

      // Iterate through each point
      std::vector<std::vector<double>> coef(
        noPoints, std::vector<double>(4, 0.0));
      for(unsigned int i = 0; i < noPoints; ++i) {
        std::vector<std::vector<double>> points(4, std::vector<double>(2));
        for(unsigned int j = 0; j < 4; ++j) {
          for(unsigned int k = 0; k < 2; ++k) {
            points[j][k] = surfaceTs[quadPoints[i][j]][k];
          }
        }

        // Find barycentric coordinates
        std::vector<std::vector<double>> trianglePoints(3);
        std::vector<int> indexes(3);
        std::vector<double> mid{(points[3][0] - points[1][0]) / 2.0,
                                (points[0][1] - points[1][1]) / 2.0};
        if(inputTs[i][0] > points[1][0] + mid[0]) {
          indexes[0] = 2;
          indexes[1] = 3;
          if(inputTs[i][1] > points[1][1] + mid[1])
            indexes[2] = 0;
          else
            indexes[2] = 1;
        } else {
          indexes[0] = 0;
          indexes[1] = 1;
          if(inputTs[i][1] > points[1][1] + mid[1])
            indexes[2] = 2;
          else
            indexes[2] = 3;
        }

        for(unsigned int j = 0; j < 3; ++j)
          trianglePoints[j] = points[indexes[j]];

        std::vector<double> area(3);
        coef[i][indexes[2]]
          = triangleArea(inputTs[i], trianglePoints[0], trianglePoints[1]);
        coef[i][indexes[1]]
          = triangleArea(inputTs[i], trianglePoints[0], trianglePoints[2]);
        coef[i][indexes[0]]
          = triangleArea(inputTs[i], trianglePoints[1], trianglePoints[2]);
        double sumArea
          = coef[i][indexes[2]] + coef[i][indexes[1]] + coef[i][indexes[0]];
        for(unsigned int j = 0; j < 3; ++j)
          coef[i][indexes[j]] /= sumArea;

        // Compute new point
        for(unsigned int j = 0; j < 3; ++j)
          for(unsigned int k = 0; k < 4; ++k)
            inputPoints[i][j]
              += coef[i][k] * surfacePoints[quadPoints[i][k]][j];
      }
    }

    //-------------------------------------------------------------------------
    // Utils
    //-------------------------------------------------------------------------
    double triangleArea(std::vector<double> &x1,
                        std::vector<double> &x2,
                        std::vector<double> &x3) {
      return std::abs((x1[0] * x2[1] + x2[0] * x3[1] + x3[0] * x1[1]
                       - x1[1] * x2[0] - x2[1] * x3[0] - x3[1] * x1[0])
                      / 2.0);
    }

    double triangleArea3D(std::vector<double> &x1,
                          std::vector<double> &x2,
                          std::vector<double> &x3) {
      std::vector<double> AB(3), AC(3);
      for(unsigned int i = 0; i < 3; ++i) {
        AB[i] = x2[i] - x1[i];
        AC[i] = x3[i] - x1[i];
      }
      return std::sqrt(std::pow(AB[1] * AC[2] - AB[2] * AC[1], 2)
                       + std::pow(AB[2] * AC[0] - AB[0] * AC[2], 2)
                       + std::pow(AB[0] * AC[1] - AB[1] * AC[0], 2))
             / 2.0;
    }

    double triangleAreaFromSides(double a, double b, double c) {
      double s = (a + b + c) / 2.0;
      return std::sqrt(s * (s - a) * (s - b) * (s - c));
    }

  }; // MergeTreePrincipalGeodesicsCurvature class

} // namespace ttk
