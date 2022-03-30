/// \class ttkPointSetToGrid
/// \ingroup vtk
/// \author Pierre Guillou <pierre.guillou@lip6.fr>
/// \date March 2020
///
/// \brief TTK VTK-filter that reads a Cinema Spec D Database.
///
/// This filter can be used as any other VTK filter (for instance, by using the
/// sequence of calls SetInputData(), Update(), GetOutput()).
///
/// \param Output content of the data.csv file of the database in form of a
/// vtkTable

#pragma once

// Module include
#include <ttkPointSetToGridModule.h>

// VTK includes
#include <ttkAlgorithm.h>

class TTKPOINTSETTOGRID_EXPORT ttkPointSetToGrid : public ttkAlgorithm {

public:
  static ttkPointSetToGrid *New();
  vtkTypeMacro(ttkPointSetToGrid, ttkAlgorithm);

protected:
  ttkPointSetToGrid();
  ~ttkPointSetToGrid() = default;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  int RequestData(vtkInformation *request,
                  vtkInformationVector **inputVector,
                  vtkInformationVector *outputVector) override;
  template <typename VTK_TT>
  void dispatch(std::vector<std::pair<vtkIdType, double>> &storage,
                const VTK_TT *const values,
                const size_t nvalues);

private:
};
