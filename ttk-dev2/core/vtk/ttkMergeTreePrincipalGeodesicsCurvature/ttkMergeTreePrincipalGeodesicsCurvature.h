/// \ingroup vtk
/// \class ttkMergeTreePrincipalGeodesicsCurvature
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \date 2022.
///
/// \brief TTK VTK-filter that wraps the
/// ttk::MergeTreePrincipalGeodesicsCurvature module.
///
/// This VTK filter uses the ttk::MergeTreePrincipalGeodesicsCurvature module to
/// compute TODO
///
/// \param Input vtkMultiBlockDataSet Trees
/// \param Input vtkUnstructuredGrid Surface
/// \param Output vtkDataSet.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// See the corresponding standalone program for a usage example:
///   - standalone/MergeTreePrincipalGeodesicsCurvature/main.cpp
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::MergeTreePrincipalGeodesicsCurvature
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkMergeTreePrincipalGeodesicsCurvatureModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkUnstructuredGrid.h>

// TTK Base Includes
#include <MergeTreePrincipalGeodesicsCurvature.h>

class TTKMERGETREEPRINCIPALGEODESICSCURVATURE_EXPORT
  ttkMergeTreePrincipalGeodesicsCurvature
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::MergeTreePrincipalGeodesicsCurvature // and we inherit from
                                                        // the base class
{
private:
  /**
   * Add all filter parameters only as private member variables and
   * initialize them here.
   */
  bool ComputeGeodesicDistance = false;
  bool ComputeSurfaceArea = false;

  // ----------------------
  // Data for visualization
  // ----------------------
  // Trees
  std::vector<ttk::ftm::MergeTree<double>> inputMTrees;
  std::vector<vtkUnstructuredGrid *> inputTreesNodes;
  std::vector<vtkUnstructuredGrid *> inputTreesArcs;
  std::vector<vtkDataSet *> inputTreesSegmentation;
  // Trees infos
  std::vector<std::vector<double>> inputTs;
  std::vector<bool> inputIsSurface;
  int noInputIsSurface;
  std::vector<int> inputCorr;
  // Surface
  std::vector<std::vector<double>> surfacePoints;
  std::vector<std::vector<double>> surfaceTs;
  // Input Surface
  vtkUnstructuredGrid *geodSurface;
  vtkMTimeType geodSurfaceMTime;
  vtkMultiBlockDataSet *blockTrees;
  // Output
  std::vector<std::vector<double>> pointsArea, wassersteinArea, ratioArea;
  std::vector<std::vector<double>> inputPoints;

  void setDataVisualization(int numTrees) {
    // Trees
    inputMTrees = std::vector<ttk::ftm::MergeTree<double>>(numTrees);
    inputTreesNodes = std::vector<vtkUnstructuredGrid *>(numTrees);
    inputTreesArcs = std::vector<vtkUnstructuredGrid *>(numTrees);
    inputTreesSegmentation = std::vector<vtkDataSet *>(numTrees);
  }

  bool isDataVisualizationFilled() {
    return (ComputeGeodesicDistance == (distances_.size() != 0))
           and (ComputeSurfaceArea
                == (pointsArea.size() != 0 and wassersteinArea.size() != 0
                    and ratioArea.size() != 0))
           and inputPoints.size() != 0;
  }

  void resetDataVisualization() {
    setDataVisualization(0);
    distances_.clear();
    pointsArea.clear();
    wassersteinArea.clear();
    ratioArea.clear();
    inputPoints.clear();
  }

public:
  /**
   * Automatically generate getters and setters of filter
   * parameters via vtkMacros.
   */
  vtkSetMacro(ComputeGeodesicDistance, bool);
  vtkGetMacro(ComputeGeodesicDistance, bool);
  vtkSetMacro(ComputeSurfaceArea, bool);
  vtkGetMacro(ComputeSurfaceArea, bool);

  /**
   * This static method and the macro below are VTK conventions on how to
   * instantiate VTK objects. You don't have to modify this.
   */
  static ttkMergeTreePrincipalGeodesicsCurvature *New();
  vtkTypeMacro(ttkMergeTreePrincipalGeodesicsCurvature, ttkAlgorithm);

protected:
  /**
   * Implement the filter constructor and destructor
   * (see cpp file)
   */
  ttkMergeTreePrincipalGeodesicsCurvature();
  ~ttkMergeTreePrincipalGeodesicsCurvature() override;

  /**
   * Specify the input data type of each input port
   * (see cpp file)
   */
  int FillInputPortInformation(int port, vtkInformation *info) override;

  /**
   * Specify the data object type of each output port
   * (see cpp file)
   */
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  /**
   * Pass VTK data to the base code and convert base code output to VTK
   * (see cpp file)
   */
  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;

  template <class dataType>
  int run(vtkInformationVector *outputVector,
          std::vector<vtkMultiBlockDataSet *> &inputTrees);

  template <class dataType>
  int runCompute(vtkInformationVector *outputVector,
                 std::vector<vtkMultiBlockDataSet *> &inputTrees);

  template <class dataType>
  int runOutput(vtkInformationVector *outputVector,
                std::vector<vtkMultiBlockDataSet *> &inputTrees);
};
