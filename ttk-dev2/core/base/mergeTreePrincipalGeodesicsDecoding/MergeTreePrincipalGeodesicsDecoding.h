/// \ingroup base
/// \class ttk::MergeTreePrincipalGeodesicsDecoding
/// \author XXX
/// \date 2021.
///
/// This module defines the %MergeTreePrincipalGeodesicsDecoding class that
/// computes
/// TODO
///
/// \b Related \b publication: \n
/// TODO
///

#pragma once

// ttk common includes
#include <Debug.h>
#include <MergeTreePrincipalGeodesicsBase.h>

namespace ttk {

  /**
   * The MergeTreePrincipalGeodesicsDecoding class provides methods to compute
   * TODO
   */
  class MergeTreePrincipalGeodesicsDecoding
    : virtual public Debug,
      public MergeTreePrincipalGeodesicsBase {

  protected:
    std::vector<std::vector<double>> tEllipses_, tRectangle_, tSurface_;
    std::vector<std::vector<std::vector<double>>> tGeodesics_;

    std::vector<double> geodesicsDistances_; // distance between extremities

    std::vector<bool> surfaceIsBoundary_;
    std::vector<int> surfaceBoundaryID_;

    bool computeReconstructionError_ = true;

  public:
    MergeTreePrincipalGeodesicsDecoding();

    //----------------------------------------------------------------------------
    // Utils
    //----------------------------------------------------------------------------
    template <class dataType>
    void preprocessBarycenter(ftm::MergeTree<dataType> &barycenter) {
      if(not isPersistenceDiagram_) {
        bool useMinMax = true;
        bool cleanTree = false;
        bool pt = 0.0;
        std::vector<int> nodeCorr;
        preprocessingPipeline<dataType>(barycenter, 0.0, 100.0, 100.0,
                                        branchDecomposition_, useMinMax,
                                        cleanTree, pt, nodeCorr, false);
      }
    }

    template <class dataType>
    void processInputTrees(
      std::vector<ttk::ftm::MergeTree<dataType>> &inputTrees) {
      treesNodeCorr_ = std::vector<std::vector<int>>(inputTrees.size());
      preprocessingTrees<dataType>(inputTrees, treesNodeCorr_);
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif
      for(unsigned int i = 0; i < inputTrees.size(); ++i)
        postprocessingPipeline<dataType>(&(inputTrees[i].tree));
    }

    // TODO Manage double input
    template <class dataType>
    void getGeodesicsMiddle(ftm::MergeTree<dataType> &barycenter,
                            std::vector<std::vector<std::vector<double>>> &vS,
                            std::vector<std::vector<std::vector<double>>> &v2s,
                            std::vector<double> &middle) {
      int cptDivide = 0;
      std::vector<double> alpha(2, 0.0);
      for(unsigned int i = 0; i < 2; ++i) {
        cptDivide = 0;
        for(unsigned int j = 0; j < vS[i].size(); ++j) {
          if(barycenter.tree.isNodeAlone(j))
            continue;
          for(unsigned int k = 0; k < 2; ++k) {
            if(std::abs(v2s[i][j][k]) < 1e-12)
              continue;
            alpha[i] += (vS[i][j][k] / v2s[i][j][k]);
            ++cptDivide;
          }
        }
        alpha[i] /= cptDivide;
      }

      middle = std::vector<double>(2, 0.0);
      for(unsigned int i = 0; i < 2; ++i)
        middle[i] = alpha[i] / (1 + alpha[i]);
    }

    template <class dataType>
    void computeGeodesicsDistance(
      std::vector<ttk::ftm::MergeTree<dataType>> &barycenters) {
      // Preprocessing
      preprocessBarycenter<dataType>(barycenters[0]);
      if(barycenters.size() > 1)
        preprocessBarycenter<dataType>(barycenters[1]);

      // Compute
      geodesicsDistances_ = std::vector<double>(vS_.size());
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif
      for(unsigned int i = 0; i < vS_.size(); ++i) {
        ftm::MergeTree<dataType> extremityV1, extremityV2;
        getInterpolation<dataType>(
          barycenters[0], vS_[i], v2s_[i], 0.0, extremityV1);
        getInterpolation<dataType>(
          barycenters[0], vS_[i], v2s_[i], 1.0, extremityV2);
        // Get distance
        dataType distance;
        computeOneDistance(
          extremityV1, extremityV2, distance, true, useDoubleInput_);
        if(barycenters.size() > 1) {
          ftm::MergeTree<dataType> extremityV1_2, extremityV2_2;
          getInterpolation<dataType>(
            barycenters[1], trees2Vs_[i], trees2V2s_[i], 0.0, extremityV1_2);
          getInterpolation<dataType>(
            barycenters[1], trees2Vs_[i], trees2V2s_[i], 1.0, extremityV2_2);
          // Get distance
          dataType distance2;
          computeOneDistance(extremityV1_2, extremityV2_2, distance2, true,
                             useDoubleInput_, false);
          distance = mixDistances(distance, distance2);
        }
        geodesicsDistances_[i] = distance;
      }

      // Postprocessing
      postprocessingPipeline<dataType>(&(barycenters[0].tree));
      if(barycenters.size() > 1)
        postprocessingPipeline<dataType>(&(barycenters[1].tree));
    }

    //----------------------------------------------------------------------------
    // Construct functions
    //----------------------------------------------------------------------------
    template <class dataType>
    void
      reconstruction(ftm::MergeTree<dataType> &barycenter,
                     std::vector<ftm::MergeTree<dataType>> &inputTrees,
                     std::vector<ftm::MergeTree<dataType>> &reconstructedTrees,
                     std::vector<double> &reconstructionErrors,
                     bool isSecondInput = false) {
      auto &vSToUse = (isSecondInput ? trees2Vs_ : vS_);
      auto &v2sToUse = (isSecondInput ? trees2V2s_ : v2s_);

      // Preprocessing
      preprocessBarycenter<dataType>(barycenter);
      if(inputTrees.size() != 0)
        preprocessingTrees<dataType>(inputTrees);

      // Reconstruction
      reconstructedTrees
        = std::vector<ftm::MergeTree<dataType>>(allTreesTs_.size());
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif
      for(unsigned int i = 0; i < reconstructedTrees.size(); ++i)
        getMultiInterpolation<dataType>(
          barycenter, vSToUse, v2sToUse, allTreesTs_[i], reconstructedTrees[i]);

      // Compute reconstruction error (if input trees are provided)
      if(inputTrees.size() != 0 and computeReconstructionError_) {
        auto reconstructionError = computeReconstructionError(
          barycenter, inputTrees, vSToUse, v2sToUse, allTreesTs_,
          reconstructionErrors);
        std::stringstream ss;
        ss << "Reconstruction Error = " << reconstructionError;
        printMsg(ss.str());
      }

      // Postprocessing
      postprocessingPipeline<dataType>(&(barycenter.tree));
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif
      for(unsigned int i = 0; i < reconstructedTrees.size(); ++i)
        postprocessingPipeline<dataType>(&(reconstructedTrees[i].tree));
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif
      for(unsigned int i = 0; i < inputTrees.size(); ++i)
        postprocessingPipeline<dataType>(&(inputTrees[i].tree));
    }

    // TODO Manage double input
    template <class dataType>
    void constructGeodesicsTrees(
      ftm::MergeTree<dataType> &barycenter,
      std::vector<std::vector<ftm::MergeTree<dataType>>> &geodesicsTrees) {

      // Preprocessing
      preprocessBarycenter<dataType>(barycenter);

      std::vector<double> middle;
      getGeodesicsMiddle<dataType>(barycenter, vS_, v2s_, middle);

      // Construct geodesics trees
      geodesicsTrees = std::vector<std::vector<ftm::MergeTree<dataType>>>(
        vS_.size(), std::vector<ftm::MergeTree<dataType>>(k_));
      tGeodesics_ = std::vector<std::vector<std::vector<double>>>(
        geodesicsTrees.size(),
        std::vector<std::vector<double>>(
          geodesicsTrees[0].size(), std::vector<double>(vS_.size(), 0.0)));
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif
      for(unsigned int i = 0; i < geodesicsTrees.size(); ++i)
        for(unsigned int j = 0; j < geodesicsTrees[i].size(); ++j) {
          double t = 1.0 / (k_ - 1) * j;
          getInterpolation<dataType>(
            barycenter, vS_[i], v2s_[i], t, geodesicsTrees[i][j]);
          tGeodesics_[i][j][i] = t;
          int i2 = (i + 1) % 2;
          tGeodesics_[i][j][i2] = middle[i2];
        }

      // Postprocessing
      postprocessingPipeline<dataType>(&(barycenter.tree));
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif
      for(unsigned int i = 0; i < geodesicsTrees.size(); ++i)
        for(unsigned int j = 0; j < geodesicsTrees[i].size(); ++j)
          postprocessingPipeline<dataType>(&(geodesicsTrees[i][j].tree));
    }

    // TODO Manage double input
    template <class dataType>
    void constructGeodesicsEllipses(
      ftm::MergeTree<dataType> &barycenter,
      std::vector<ftm::MergeTree<dataType>> &geodesicsEllipses) {

      // Preprocessing
      preprocessBarycenter<dataType>(barycenter);

      // Init
      unsigned int noSample = k_ * 2;
      geodesicsEllipses = std::vector<ftm::MergeTree<dataType>>(noSample);
      tEllipses_ = std::vector<std::vector<double>>(geodesicsEllipses.size());

      std::vector<std::vector<std::vector<double>>> vS(2), v2s(2);
      vS[0] = vS_[0];
      vS[1] = vS_[1];
      v2s[0] = v2s_[0];
      v2s[1] = v2s_[1];

      // Get middle of geodesic
      std::vector<double> middle;
      getGeodesicsMiddle<dataType>(barycenter, vS, v2s, middle);

      // Get ellipses
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif
      for(unsigned int i = 0; i < noSample; ++i) {
        double angle = 360.0 / noSample * i;
        double pi = 3.14159265359;
        double radius = 1.0;
        double x = -1 * radius * std::cos(-1 * angle * pi / 180);
        double y = -1 * radius * std::sin(-1 * angle * pi / 180);

        // 0: upper-left ; 1: upper-right ; 2: bottom-right ; 3: bottom-left
        int quadrant = (x < 0.0 ? (y > 0.0 ? 0 : 3) : (y > 0.0 ? 1 : 2));
        if(quadrant == 0 or quadrant == 3)
          x = (x + 1.0) * middle[0];
        if(quadrant == 1 or quadrant == 2)
          x = x * (1.0 - middle[0]) + middle[0];
        if(quadrant == 0 or quadrant == 1)
          y = y * (1.0 - middle[1]) + middle[1];
        if(quadrant == 2 or quadrant == 3)
          y = (y + 1.0) * middle[1];

        // Get interpolation
        if(x > 1.0 or y > 1.0 or x < 0.0 or y < 0.0)
          printErr("[constructGeodesicsEllipses] extrapolation.");
        std::vector<double> ts{x, y};
        getMultiInterpolation<dataType>(
          barycenter, vS, v2s, ts, geodesicsEllipses[i]);
        tEllipses_[i] = ts;
      }

      // Postprocessing
      postprocessingPipeline<dataType>(&(barycenter.tree));
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif
      for(unsigned int i = 0; i < geodesicsEllipses.size(); ++i)
        postprocessingPipeline<dataType>(&(geodesicsEllipses[i].tree));
    }

    unsigned int getNumberOfRectangles(unsigned int rectangleMultiplier = 1) {
      return k_ * 4 * rectangleMultiplier;
    }

    // TODO Manage double input
    template <class dataType>
    void constructGeodesicsRectangle(
      ftm::MergeTree<dataType> &barycenter,
      std::vector<ftm::MergeTree<dataType>> &geodesicsRectangle,
      unsigned int rectangleMultiplier = 1) {

      // Preprocessing
      preprocessBarycenter<dataType>(barycenter);

      // Init
      unsigned int noSample = getNumberOfRectangles(rectangleMultiplier);
      geodesicsRectangle = std::vector<ftm::MergeTree<dataType>>(noSample);
      tRectangle_ = std::vector<std::vector<double>>(geodesicsRectangle.size());

      std::vector<std::vector<std::vector<double>>> vS(2), v2s(2);
      vS[0] = vS_[0];
      vS[1] = vS_[1];
      v2s[0] = v2s_[0];
      v2s[1] = v2s_[1];

      // Get ellipses
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif
      for(unsigned int i = 0; i < noSample; ++i) {
        // 0: left ; 1: up; 2: right ; 3: bottom
        int quadrant = i / (noSample / 4);

        double x = 0.0, y = 0.0;
        double offset = (i % (noSample / 4)) / ((noSample / 4.0) - 1.0);
        switch(quadrant) {
          case 0: // Left
            x = 0.0;
            y = 0.0 + offset;
            break;
          case 1: // Up
            x = 0.0 + offset;
            y = 1.0;
            break;
          case 2: // Right
            x = 1.0;
            y = 1.0 - offset;
            break;
          case 3: // Bottom
          default:
            x = 1.0 - offset;
            y = 0.0;
        }

        // Get interpolation
        std::vector<double> ts{x, y};
        getMultiInterpolation<dataType>(
          barycenter, vS, v2s, ts, geodesicsRectangle[i]);
        tRectangle_[i] = ts;
      }

      // Postprocessing
      postprocessingPipeline<dataType>(&(barycenter.tree));
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif
      for(unsigned int i = 0; i < geodesicsRectangle.size(); ++i)
        postprocessingPipeline<dataType>(&(geodesicsRectangle[i].tree));
    }

    template <class dataType>
    void constructGeodesicsSurface(
      ftm::MergeTree<dataType> &barycenter,
      std::vector<ftm::MergeTree<dataType>> &geodesicsSurface,
      bool isSecondInput = false) {

      // Preprocessing
      preprocessBarycenter<dataType>(barycenter);

      // Init
      unsigned int noSample = k_ * k_;
      geodesicsSurface = std::vector<ftm::MergeTree<dataType>>(noSample);
      tSurface_ = std::vector<std::vector<double>>(geodesicsSurface.size());
      surfaceIsBoundary_ = std::vector<bool>(geodesicsSurface.size(), false);
      surfaceBoundaryID_ = std::vector<int>(geodesicsSurface.size(), -1);

      std::vector<std::vector<std::vector<double>>> vS(2), v2s(2);
      vS[0] = (isSecondInput ? trees2Vs_[0] : vS_[0]);
      vS[1] = (isSecondInput ? trees2Vs_[1] : vS_[1]);
      v2s[0] = (isSecondInput ? trees2V2s_[0] : v2s_[0]);
      v2s[1] = (isSecondInput ? trees2V2s_[1] : v2s_[1]);

      // Get surface
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif
      for(unsigned int i = 0; i < k_; ++i) {
        for(unsigned int j = 0; j < k_; ++j) {
          double x = 1.0 / (k_ - 1) * i;
          double y = 1.0 / (k_ - 1) * j;
          int index = i * k_ + j;
          // Get Boundary ID
          if(i == 0 or j == 0 or i == k_ - 1 or j == k_ - 1) {
            surfaceIsBoundary_[index] = true;
            int boundaryID = 0;
            if(i == 0)
              boundaryID = j;
            else if(j == k_ - 1)
              boundaryID = k_ + i;
            else if(i == k_ - 1)
              boundaryID = k_ * 2 + (k_ - j);
            else // if(j == 0)
              boundaryID = k_ * 3 + (k_ - i);
            surfaceBoundaryID_[index] = boundaryID;
          }
          // Get interpolation
          std::vector<double> ts{x, y};
          getMultiInterpolation<dataType>(
            barycenter, vS, v2s, ts, geodesicsSurface[index]);
          tSurface_[index] = ts;
        }
      }

      // Postprocessing
      postprocessingPipeline<dataType>(&(barycenter.tree));
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(this->threadNumber_)
#endif
      for(unsigned int i = 0; i < geodesicsSurface.size(); ++i)
        postprocessingPipeline<dataType>(&(geodesicsSurface[i].tree));
    }
  }; // MergeTreePrincipalGeodesicsDecoding class

} // namespace ttk
