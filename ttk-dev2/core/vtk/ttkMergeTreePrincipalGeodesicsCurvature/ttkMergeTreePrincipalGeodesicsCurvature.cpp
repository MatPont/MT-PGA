#include <ttkFTMTreeUtils.h>
#include <ttkMergeTreePrincipalGeodesicsCurvature.h>

#include <vtkInformation.h>

#include <vtkCellType.h>
#include <vtkDataArray.h>
#include <vtkDataSet.h>
#include <vtkObjectFactory.h>
#include <vtkPointData.h>
#include <vtkSmartPointer.h>
#include <vtkStringArray.h>

#include <ttkMacros.h>
#include <ttkUtils.h>

// A VTK macro that enables the instantiation of this class via ::New()
// You do not have to modify this
vtkStandardNewMacro(ttkMergeTreePrincipalGeodesicsCurvature);

/**
 * Implement the filter constructor and destructor in the cpp file.
 *
 * The constructor has to specify the number of input and output ports
 * with the functions SetNumberOfInputPorts and SetNumberOfOutputPorts,
 * respectively. It should also set default values for all filter
 * parameters.
 *
 * The destructor is usually empty unless you want to manage memory
 * explicitly, by for example allocating memory on the heap that needs
 * to be freed when the filter is destroyed.
 */
ttkMergeTreePrincipalGeodesicsCurvature::
  ttkMergeTreePrincipalGeodesicsCurvature() {
  this->SetNumberOfInputPorts(2);
  this->SetNumberOfOutputPorts(1);
}

ttkMergeTreePrincipalGeodesicsCurvature::
  ~ttkMergeTreePrincipalGeodesicsCurvature() {
}

/**
 * Specify the required input data type of each input port
 *
 * This method specifies the required input object data types of the
 * filter by adding the vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE() key to
 * the port information.
 */
int ttkMergeTreePrincipalGeodesicsCurvature::FillInputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkMultiBlockDataSet");
  } else if(port == 1) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
  } else
    return 0;
  return 1;
}

/**
 * Specify the data object type of each output port
 *
 * This method specifies in the port information object the data type of the
 * corresponding output objects. It is possible to either explicitly
 * specify a type by adding a vtkDataObject::DATA_TYPE_NAME() key:
 *
 *      info->Set( vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid" );
 *
 * or to pass a type of an input port to an output port by adding the
 * ttkAlgorithm::SAME_DATA_TYPE_AS_INPUT_PORT() key (see below).
 *
 * Note: prior to the execution of the RequestData method the pipeline will
 * initialize empty output data objects based on this information.
 */
int ttkMergeTreePrincipalGeodesicsCurvature::FillOutputPortInformation(
  int port, vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
  } else
    return 0;
  return 1;
}

/**
 * Pass VTK data to the base code and convert base code output to VTK
 *
 * This method is called during the pipeline execution to update the
 * already initialized output data objects based on the given input
 * data objects and filter parameters.
 *
 * Note:
 *     1) The passed input data objects are validated based on the information
 *        provided by the FillInputPortInformation method.
 *     2) The output objects are already initialized based on the information
 *        provided by the FillOutputPortInformation method.
 */
int ttkMergeTreePrincipalGeodesicsCurvature::RequestData(
  vtkInformation *ttkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector) {
  // --------------------------------------------------------------------------
  // --- Get input object from input vector
  // --------------------------------------------------------------------------
  blockTrees = vtkMultiBlockDataSet::GetData(inputVector[0], 0);
  geodSurface = vtkUnstructuredGrid::GetData(inputVector[1], 0);

  // --------------------------------------------------------------------------
  // --- Load blocks
  // --------------------------------------------------------------------------
  std::vector<vtkMultiBlockDataSet *> inputTrees;
  ttk::ftm::loadBlocks(inputTrees, blockTrees);
  auto arrayToGet
    = vtkUnstructuredGrid::SafeDownCast(inputTrees[0]->GetBlock(0))
        ->GetPointData()
        ->GetArray("Scalar");
  if(arrayToGet == nullptr)
    arrayToGet = vtkUnstructuredGrid::SafeDownCast(inputTrees[0]->GetBlock(0))
                   ->GetPointData()
                   ->GetArray("Birth");
  int dataTypeInt = arrayToGet->GetDataType();

  // If we have already computed once but the input has changed
  if((inputTreesNodes.size() != 0
      and inputTrees[0]->GetBlock(0) != inputTreesNodes[0])
     or inputTrees.size() != inputTreesNodes.size()
     or (geodSurfaceMTime != geodSurface->GetMeshMTime()))
    resetDataVisualization();
  geodSurfaceMTime = geodSurface->GetMeshMTime();

  // Parameters
  printMsg("Load parameters from field data.");
  std::vector<std::string> paramNames;
  getParamNames(paramNames);
  for(auto paramName : paramNames) {
    auto array = blockTrees->GetFieldData()->GetArray(paramName.c_str());
    if(array) {
      double value = array->GetTuple1(0);
      setParamValueFromName(paramName, value);
      printMsg(" - " + paramName + " = " + std::to_string(value));
    } else
      printMsg(" - " + paramName + " was not found in the field data.");
  }
  if(normalizedWasserstein_)
    printMsg("Computation with normalized Wasserstein.");
  else
    printMsg("Computation without normalized Wasserstein.");
  if(ComputeSurfaceArea)
    ComputeGeodesicDistance = false;

  // Get number of surface trees
  inputIsSurface = std::vector<bool>(inputTrees.size(), false);
  noInputIsSurface = 0;
  for(unsigned int i = 0; i < inputTrees.size(); ++i) {
    auto fd = inputTrees[i]->GetBlock(0)->GetFieldData();
    inputIsSurface[i]
      = fd->GetAbstractArray("isSurface")->GetVariantValue(0).ToDouble();
    noInputIsSurface += inputIsSurface[i];
  }

  // Load input ts
  unsigned int numberOfGeodesics = 2;
  unsigned int noInputs = inputTrees.size() - noInputIsSurface;
  inputTs = std::vector<std::vector<double>>(
    noInputs, std::vector<double>(numberOfGeodesics));
  inputCorr = std::vector<int>(noInputs, -1);
  int inputCpt = 0;
  for(unsigned int i = 0; i < inputTrees.size(); ++i) {
    if(inputIsSurface[i])
      continue;
    auto fd = inputTrees[i]->GetBlock(0)->GetFieldData();
    for(unsigned int j = 0; j < numberOfGeodesics; ++j) {
      std::string name = getTableCoefficientName(numberOfGeodesics, j);
      inputTs[inputCpt][j]
        = fd->GetAbstractArray(name.c_str())->GetVariantValue(0).ToDouble();
    }
    inputCorr[inputCpt] = i;
    ++inputCpt;
  }

  // --------------------------------------------------------------------------
  // --- Load Surface
  // --------------------------------------------------------------------------
  auto noPoints = geodSurface->GetNumberOfPoints();
  surfacePoints
    = std::vector<std::vector<double>>(noPoints, std::vector<double>(3));
  surfaceTs = std::vector<std::vector<double>>(
    noPoints, std::vector<double>(numberOfGeodesics));
  for(unsigned int i = 0; i < noPoints; ++i) {
    double point[3];
    geodSurface->GetPoints()->GetPoint(i, point);
    for(unsigned int j = 0; j < 3; ++j)
      surfacePoints[i][j] = point[j];
    for(unsigned int j = 0; j < numberOfGeodesics; ++j) {
      std::string name = getTableCoefficientName(numberOfGeodesics, j);
      surfaceTs[i][j] = geodSurface->GetPointData()
                          ->GetAbstractArray(name.c_str())
                          ->GetVariantValue(i)
                          .ToDouble();
    }
  }

  std::vector<std::tuple<int, std::vector<double>>> pointsAndIds(noPoints);
  std::vector<std::tuple<int, std::vector<double>>> tsAndIds(noPoints);
  for(unsigned int i = 0; i < noPoints; ++i) {
    auto treeIdArray = geodSurface->GetPointData()->GetAbstractArray("TreeID");
    auto treeId = treeIdArray->GetVariantValue(i).ToInt();
    pointsAndIds[i] = std::make_tuple(treeId, surfacePoints[i]);
    tsAndIds[i] = std::make_tuple(treeId, surfaceTs[i]);
  }
  std::sort(pointsAndIds.begin(), pointsAndIds.end());
  std::sort(tsAndIds.begin(), tsAndIds.end());
  for(unsigned int i = 0; i < noPoints; ++i) {
    surfacePoints[i] = std::get<1>(pointsAndIds[i]);
    surfaceTs[i] = std::get<1>(tsAndIds[i]);
  }

  // --------------------------------------------------------------------------

  int res = 0;
  switch(dataTypeInt) {
    vtkTemplateMacro(res = run<VTK_TT>(outputVector, inputTrees););
  }
  return res;
}

template <class dataType>
int ttkMergeTreePrincipalGeodesicsCurvature::run(
  vtkInformationVector *outputVector,
  std::vector<vtkMultiBlockDataSet *> &inputTrees) {
  if(not isDataVisualizationFilled())
    runCompute<dataType>(outputVector, inputTrees);
  runOutput<dataType>(outputVector, inputTrees);
  return 1;
}

template <class dataType>
int ttkMergeTreePrincipalGeodesicsCurvature::runCompute(
  vtkInformationVector *ttkNotUsed(outputVector),
  std::vector<vtkMultiBlockDataSet *> &inputTrees) {
  // --------------------------------------------------------------------------
  // --- Construct trees
  // --------------------------------------------------------------------------
  const int numTrees = inputTrees.size();

  setDataVisualization(numTrees);

  std::vector<ttk::ftm::MergeTree<dataType>> inputDTrees(numTrees);

  isPersistenceDiagram_ = ttk::ftm::constructTrees<dataType>(
    inputTrees, inputDTrees, inputTreesNodes, inputTreesArcs,
    inputTreesSegmentation);

  //---------------------------------------------------------------------------
  // --- Call base
  //---------------------------------------------------------------------------
  if(ComputeGeodesicDistance or ComputeSurfaceArea) {
    // Get surface trees̉
    std::vector<ttk::ftm::MergeTree<dataType>> surfaceTrees(noInputIsSurface);
    int cpt = 0;
    for(unsigned int i = 0; i < inputIsSurface.size(); ++i)
      if(inputIsSurface[i]) {
        surfaceTrees[cpt] = inputDTrees[i];
        ++cpt;
      }

    // Compute
    if(ComputeGeodesicDistance and not ComputeSurfaceArea) {
      computeSurfaceDistances<dataType>(surfaceTrees);
    } else if(not ComputeGeodesicDistance and not ComputeSurfaceArea)
      distances_.clear();
    if(ComputeSurfaceArea) {
      computeSurfaceArea(
        surfaceTrees, surfacePoints, pointsArea, wassersteinArea, ratioArea);
    } else {
      pointsArea.clear();
      wassersteinArea.clear();
      ratioArea.clear();
      distances_.clear();
    }
  }
  computeInputPoints(inputTs, surfacePoints, surfaceTs, inputPoints);

  return 1;
}

template <class dataType>
int ttkMergeTreePrincipalGeodesicsCurvature::runOutput(
  vtkInformationVector *outputVector,
  std::vector<vtkMultiBlockDataSet *> &inputTrees) {
  auto output = vtkUnstructuredGrid::GetData(outputVector, 0);
  output->DeepCopy(geodSurface);

  std::set<std::string> toRemove;
  for(unsigned int i = 0; i < inputPoints.size(); ++i) {
    double point[3];
    for(unsigned int j = 0; j < 3; ++j) {
      point[j] = inputPoints[i][j];
    }
    output->GetPoints()->InsertNextPoint(point);

    // Merge data of decoded trees and geodesic surface
    for(int j = 0; j < output->GetPointData()->GetNumberOfArrays(); ++j) {
      auto outputArray = output->GetPointData()->GetAbstractArray(j);
      std::string name = outputArray->GetName();
      auto array = inputTrees[inputCorr[i]]
                     ->GetBlock(0)
                     ->GetFieldData()
                     ->GetAbstractArray(name.c_str());
      /*if(not array)
        array = inputTrees[i]->GetFieldData()->GetAbstractArray(
          name.c_str());*/
      if(array) {
        outputArray->InsertNextTuple(0, array);
      } else {
        auto dataArray = vtkDataArray::SafeDownCast(outputArray);
        if(dataArray) {
          const double val = std::nan("");
          dataArray->InsertNextTuple(&val);
        } else {
          auto stringArray = vtkStringArray::SafeDownCast(outputArray);
          if(stringArray)
            stringArray->InsertNextValue("");
          else
            toRemove.insert(name);
        }
      }
    }
  }

  // Add empty data for field data in decoded trees not in geodesic surface
  std::set<std::string> toGet;
  for(int j = 0;
      j < inputTrees[0]->GetBlock(0)->GetFieldData()->GetNumberOfArrays();
      ++j) {
    auto inputArray
      = inputTrees[0]->GetBlock(0)->GetFieldData()->GetAbstractArray(j);
    auto dataArray = vtkDataArray::SafeDownCast(inputArray);
    auto stringArray = vtkStringArray::SafeDownCast(inputArray);
    std::string name = inputArray->GetName();
    auto array = output->GetPointData()->GetAbstractArray(name.c_str());
    if(not array and (dataArray or stringArray))
      toGet.insert(name);
  }
  vtkNew<vtkFieldData> fd{};
  if(toGet.size() != 0)
    fd->DeepCopy(inputTrees[0]->GetBlock(0)->GetFieldData());
  for(auto &name : toGet) {
    auto inputArray = fd->GetAbstractArray(name.c_str());
    auto dataArray = vtkDataArray::SafeDownCast(inputArray);
    auto stringArray = vtkStringArray::SafeDownCast(inputArray);
    inputArray->SetNumberOfTuples(output->GetNumberOfPoints());
    int noPointsOri = output->GetNumberOfPoints() - inputPoints.size();
    for(int i = 0; i < output->GetNumberOfPoints(); ++i) {
      int inputCorrIndex = i - noPointsOri;
      // if is surface or tree without some fields (like polyLine trees)
      if(i < noPointsOri
         or (inputCorrIndex >= 0
             and not inputTrees[inputCorr[i - noPointsOri]]
                       ->GetBlock(0)
                       ->GetFieldData()
                       ->GetAbstractArray(name.c_str()))) {
        if(dataArray) {
          double val = std::nan("");
          dataArray->SetTuple1(i, val);
        } else if(stringArray) {
          stringArray->SetValue(i, "");
        } else
          printWrn("Can not convert " + name + " to dataArray or stringArray.");
      } else {
        int index = inputCorr[inputCorrIndex];
        inputArray->SetTuple(
          i, 0,
          inputTrees[index]->GetBlock(0)->GetFieldData()->GetAbstractArray(
            name.c_str()));
      }
    }
    output->GetPointData()->AddArray(inputArray);
  }

  for(auto &name : toRemove)
    output->GetPointData()->RemoveArray(name.c_str());

  // Cells
  auto noCells = output->GetNumberOfCells();
  vtkNew<vtkIntArray> isLineArray{};
  isLineArray->SetName("isLine");
  isLineArray->SetNumberOfTuples(noCells);
  vtkNew<vtkIntArray> isQuadArray{};
  isQuadArray->SetName("isQuad");
  isQuadArray->SetNumberOfTuples(noCells);

  vtkNew<vtkDoubleArray> surfaceAreaArray{};
  surfaceAreaArray->SetName("SurfaceArea");
  surfaceAreaArray->SetNumberOfTuples(noCells);
  vtkNew<vtkDoubleArray> wassersteinAreaArray{};
  wassersteinAreaArray->SetName("WassersteinArea");
  wassersteinAreaArray->SetNumberOfTuples(noCells);
  vtkNew<vtkDoubleArray> ratioAreaArray{};
  ratioAreaArray->SetName("RatioArea");
  ratioAreaArray->SetNumberOfTuples(noCells);

  vtkNew<vtkDoubleArray> wassersteinDistanceArray{};
  wassersteinDistanceArray->SetName("WassersteinDistance");
  wassersteinDistanceArray->SetNumberOfTuples(noCells);

  for(int i = 0; i < noCells; ++i) {
    auto type = output->GetCell(i)->GetCellType();
    bool isLine = (type == VTK_LINE);
    isLineArray->SetTuple1(i, isLine);
    bool isQuad = (type == VTK_QUAD);
    isQuadArray->SetTuple1(i, isQuad);

    std::vector<int> pointIds;
    auto pointIdsVTK = output->GetCell(i)->GetPointIds();
    for(unsigned int j = 0; j < pointIdsVTK->GetNumberOfIds(); ++j)
      pointIds.push_back(pointIdsVTK->GetId(j));
    std::sort(pointIds.begin(), pointIds.end());

    unsigned int dim = std::sqrt(noInputIsSurface); // k
    int i1 = pointIds[0] / dim;
    int i2 = pointIds[0] % dim;

    // Areas
    double surfaceAreaValue = std::nan("");
    double wassersteinAreaValue = std::nan("");
    double ratioAreaValue = std::nan("");
    if(isQuad and pointsArea.size() != 0 and wassersteinArea.size() != 0
       and ratioArea.size() != 0) {
      surfaceAreaValue = pointsArea[i1][i2];
      wassersteinAreaValue = wassersteinArea[i1][i2];
      ratioAreaValue = ratioArea[i1][i2];
    }
    surfaceAreaArray->SetTuple1(i, surfaceAreaValue);
    wassersteinAreaArray->SetTuple1(i, wassersteinAreaValue);
    ratioAreaArray->SetTuple1(i, ratioAreaValue);

    // Distances
    double wassersteinDistance = std::nan("");
    if(isLine and distances_.size() != 0) {
      int index = ((pointIds[0] + 1) == pointIds[1] ? 0 : 1);
      wassersteinDistance = distances_[i1][i2][index];
    }
    wassersteinDistanceArray->SetTuple1(i, wassersteinDistance);
  }
  output->GetCellData()->AddArray(isLineArray);
  output->GetCellData()->AddArray(isQuadArray);
  if(ComputeSurfaceArea) {
    output->GetCellData()->AddArray(surfaceAreaArray);
    output->GetCellData()->AddArray(wassersteinAreaArray);
    output->GetCellData()->AddArray(ratioAreaArray);
  }
  if(ComputeGeodesicDistance or ComputeSurfaceArea) {
    output->GetCellData()->AddArray(wassersteinDistanceArray);
  }

  /*
  std::vector<std::vector<double>> pointsArea, tsArea, ratioArea;
  std::vector<std::vector<std::vector<double>>> distances;
  */

  return 1;
}
