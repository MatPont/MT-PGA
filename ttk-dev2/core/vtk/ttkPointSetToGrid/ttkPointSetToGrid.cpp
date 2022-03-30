#include <ttkPointSetToGrid.h>

#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkPointSet.h>
#include <vtkUnstructuredGrid.h>

#include <vtkInformation.h>
#include <vtkInformationVector.h>

#include <ttkUtils.h>

#include <array>
#include <map>
#include <set>

// TODO use two parameters

vtkStandardNewMacro(ttkPointSetToGrid);

ttkPointSetToGrid::ttkPointSetToGrid() {
  this->setDebugMsgPrefix("PointSetToGrid");

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);
}

int ttkPointSetToGrid::FillInputPortInformation(int port,
                                                vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPointSet");
    return 1;
  }
  return 0;
}

int ttkPointSetToGrid::FillOutputPortInformation(int port,
                                                 vtkInformation *info) {
  if(port == 0) {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid");
    return 1;
  }
  return 0;
}

template <typename VTK_TT>
void ttkPointSetToGrid::dispatch(
  std::vector<std::pair<vtkIdType, double>> &storage,
  const VTK_TT *const values,
  const size_t nvalues) {

  for(size_t i = 0; i < nvalues; ++i) {
    storage.emplace_back(i, static_cast<double>(values[i]));
  }
}

int ttkPointSetToGrid::RequestData(vtkInformation *ttkNotUsed(request),
                                   vtkInformationVector **inputVector,
                                   vtkInformationVector *outputVector) {
  const auto input = vtkPointSet::GetData(inputVector[0]);
  auto output = vtkUnstructuredGrid::GetData(outputVector);

  if(input == nullptr || output == nullptr) {
    this->printErr("Null input data, aborting");
    return 0;
  }

  // ordering array
  const auto oa = this->GetInputArrayToProcess(0, inputVector);

  if(oa == nullptr) {
    this->printErr("Cannot find the required data array");
    return 0;
  }

  const auto nvalues = oa->GetNumberOfTuples();

  // store point index <-> ordering value in vector
  std::vector<std::pair<vtkIdType, double>> orderedValues{};

  switch(oa->GetDataType()) {
    vtkTemplateMacro(
      dispatch(orderedValues,
               static_cast<VTK_TT *>(ttkUtils::GetVoidPointer(oa)), nvalues));
  }

  // compare two pairs of index/value according to their values
  const auto cmp
    = [](const std::pair<vtkIdType, double> &a,
         const std::pair<vtkIdType, double> &b) { return a.second < b.second; };

  // sort the vector of indices/values in ascending order
  std::sort(orderedValues.begin(), orderedValues.end(), cmp);

  // deep-copy input into output
  output->DeepCopy(input);

  unsigned int k = std::sqrt(nvalues);

  for(unsigned int i = 0; i < k; ++i) {
    for(unsigned int j = 0; j < k; ++j) {
      int index = i * k + j;
      int index2 = (i - 1) * k + j;
      if(j != 0) {
        std::array<vtkIdType, 2> linePoints{
          orderedValues[index - 1].first, orderedValues[index].first};
        output->InsertNextCell(VTK_LINE, 2, linePoints.data());
      }

      if(i != 0) {
        std::array<vtkIdType, 2> linePoints{
          orderedValues[index2].first, orderedValues[index].first};
        output->InsertNextCell(VTK_LINE, 2, linePoints.data());
      }

      if(i != 0 and j != 0) {
        std::array<vtkIdType, 4> cellPoints{
          orderedValues[index2 - 1].first, orderedValues[index2].first,
          orderedValues[index].first, orderedValues[index - 1].first};
        output->InsertNextCell(VTK_QUAD, 4, cellPoints.data());
      }
    }
  }

  auto noCells = output->GetNumberOfCells();
  vtkNew<vtkIntArray> cellTypeArray{};
  cellTypeArray->SetName("CellType");
  cellTypeArray->SetNumberOfTuples(noCells);
  for(int i = 0; i < noCells; ++i) {
    cellTypeArray->SetTuple1(i, output->GetCellType(i));
  }
  output->GetCellData()->AddArray(cellTypeArray);

  return 1;
}
