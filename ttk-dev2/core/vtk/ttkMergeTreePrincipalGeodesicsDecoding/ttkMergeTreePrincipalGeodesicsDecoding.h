/// \ingroup vtk
/// \class ttkMergeTreePrincipalGeodesicsDecoding
/// \author Mathieu Pont <mathieu.pont@lip6.fr>
/// \date 2022.
///
/// \brief TTK VTK-filter that wraps the
/// ttk::MergeTreePrincipalGeodesicsDecoding module.
///
/// This VTK filter uses the ttk::MergeTreePrincipalGeodesicsDecoding module to
/// compute
/// TODO
///
/// \param Input vtkMultiBlockDataSet Barycenter
/// \param Input vtkMultiBlockDataSet Geodesics
/// \param Input vtkTable Coefficients
/// \param Input vtkTable Branches Correlation
/// \param Output vtkMultiBlockDataSet Trees
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutputDataObject()).
///
/// See the corresponding standalone program for a usage example:
///   - standalone/MergeTreePrincipalGeodesicsDecoding/main.cpp
///
/// See the related ParaView example state files for usage examples within a
/// VTK pipeline.
///
/// \sa ttk::MergeTreePrincipalGeodesicsDecoding
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkMergeTreePrincipalGeodesicsDecodingModule.h>

// VTK Includes
#include <ttkAlgorithm.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkTable.h>
#include <vtkUnstructuredGrid.h>

// TTK Base Includes
#include <MergeTreePrincipalGeodesicsDecoding.h>
#include <ttkMergeTreeVisualization.h>

class TTKMERGETREEPRINCIPALGEODESICSDECODING_EXPORT
  ttkMergeTreePrincipalGeodesicsDecoding
  : public ttkAlgorithm // we inherit from the generic ttkAlgorithm class
  ,
    protected ttk::MergeTreePrincipalGeodesicsDecoding // and we inherit from
                                                       // the base class
{
private:
  /**
   * Add all filter parameters only as private member variables and
   * initialize them here.
   */

  // Output options
  bool OutputInputTrees = false;
  bool OutputBarycenter = true;
  bool ReconstructInputTrees = true;
  bool ComputeReconstructionError = true;
  bool ConstructGeodesicsTrees = false;
  bool ConstructEllipses = false;
  bool ConstructRectangle = false;
  unsigned int RectangleMultiplier = 1;
  bool ConstructSurface = false;
  bool ProcessSecondInput = false;

  // ----------------------
  // Data for visualization
  // ----------------------
  // Trees
  std::vector<ttk::ftm::MergeTree<double>> baryMTree, inputMTrees;
  std::vector<vtkUnstructuredGrid *> baryTreeNodes, inputTreesNodes;
  std::vector<vtkUnstructuredGrid *> baryTreeArcs, inputTreesArcs;
  std::vector<vtkDataSet *> baryTreeSegmentation, inputTreesSegmentation;
  // Output
  std::vector<ttk::ftm::MergeTree<double>> reconstructedTrees,
    reconstructedTrees2;
  std::vector<double> reconstructionErrors, reconstructionErrors2;
  std::vector<std::vector<ttk::ftm::MergeTree<double>>> geodesicsTrees;
  std::vector<ttk::ftm::MergeTree<double>> geodesicsEllipses,
    geodesicsRectangle, geodesicsSurface, geodesicsSurface2;
  vtkFieldData *inputFieldData;
  // TODO remove this and copy all row data of correlation table instead
  std::vector<int> baryPointScalarFieldName;
  // Verify input changed
  vtkMTimeType tableCoefficientsMTime, tableVectorsMTime, tableCorrelationMTime,
    blockInputTreesMTime;

  void setDataVisualization(int numTrees) {
    // Trees
    baryMTree = std::vector<ttk::ftm::MergeTree<double>>(1);
    baryTreeNodes = std::vector<vtkUnstructuredGrid *>(1);
    baryTreeArcs = std::vector<vtkUnstructuredGrid *>(1);
    baryTreeSegmentation = std::vector<vtkDataSet *>(1);

    inputMTrees = std::vector<ttk::ftm::MergeTree<double>>(numTrees);
    inputTreesNodes = std::vector<vtkUnstructuredGrid *>(numTrees);
    inputTreesArcs = std::vector<vtkUnstructuredGrid *>(numTrees);
    inputTreesSegmentation = std::vector<vtkDataSet *>(numTrees);
  }

  bool isDataVisualizationFilled() {
    return (ReconstructInputTrees == (reconstructedTrees.size() != 0))
           and (ConstructGeodesicsTrees == (geodesicsTrees.size() != 0))
           and (ConstructEllipses == (geodesicsEllipses.size() != 0))
           and (ConstructRectangle == (geodesicsRectangle.size() != 0))
           and (ConstructSurface == (geodesicsSurface.size() != 0))
           and (not OutputInputTrees) and geodesicsDistances_.size() != 0;
  }

  void resetDataVisualization() {
    setDataVisualization(0);
    reconstructedTrees.clear();
    geodesicsTrees.clear();
    geodesicsEllipses.clear();
    geodesicsRectangle.clear();
    geodesicsSurface.clear();
    geodesicsDistances_.clear();
  }

public:
  /**
   * Automatically generate getters and setters of filter
   * parameters via vtkMacros.
   */
  // Input options
  void SetNumberOfGeodesicsIntervals(int k) {
    k_ = k;
    Modified();
    resetDataVisualization();
  }

  int GetNumberOfGeodesicsIntervals() {
    return k_;
  }

  // Output options
  vtkSetMacro(OutputInputTrees, bool);
  vtkGetMacro(OutputInputTrees, bool);

  vtkSetMacro(OutputBarycenter, bool);
  vtkGetMacro(OutputBarycenter, bool);

  vtkSetMacro(ReconstructInputTrees, bool);
  vtkGetMacro(ReconstructInputTrees, bool);

  vtkSetMacro(computeReconstructionError_, bool);
  vtkGetMacro(computeReconstructionError_, bool);

  vtkSetMacro(ConstructGeodesicsTrees, bool);
  vtkGetMacro(ConstructGeodesicsTrees, bool);

  vtkSetMacro(ConstructEllipses, bool);
  vtkGetMacro(ConstructEllipses, bool);

  vtkSetMacro(ConstructRectangle, bool);
  vtkGetMacro(ConstructRectangle, bool);

  void SetRectangleMultiplier(unsigned int mult) {
    RectangleMultiplier = mult;
    Modified();
    geodesicsRectangle.clear();
  }
  vtkGetMacro(RectangleMultiplier, unsigned int);

  vtkSetMacro(ConstructSurface, bool);
  vtkGetMacro(ConstructSurface, bool);

  vtkSetMacro(ProcessSecondInput, bool);
  vtkGetMacro(ProcessSecondInput, bool);

  /**
   * This static method and the macro below are VTK conventions on how to
   * instantiate VTK objects. You don't have to modify this.
   */
  static ttkMergeTreePrincipalGeodesicsDecoding *New();
  vtkTypeMacro(ttkMergeTreePrincipalGeodesicsDecoding, ttkAlgorithm);

protected:
  /**
   * Implement the filter constructor and destructor
   * (see cpp file)
   */
  ttkMergeTreePrincipalGeodesicsDecoding();
  ~ttkMergeTreePrincipalGeodesicsDecoding() override;

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
          std::vector<vtkMultiBlockDataSet *> &inputBary,
          std::vector<vtkMultiBlockDataSet *> &inputTrees);

  template <class dataType>
  int runCompute(vtkInformationVector *outputVector,
                 std::vector<vtkMultiBlockDataSet *> &inputBary,
                 std::vector<vtkMultiBlockDataSet *> &inputTrees);

  template <class dataType>
  int runOutput(vtkInformationVector *outputVector,
                std::vector<vtkMultiBlockDataSet *> &inputBary,
                std::vector<vtkMultiBlockDataSet *> &inputTrees);
};
