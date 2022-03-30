/// \ingroup vtk
/// \class ttk::ttkFTMTreeUtils
/// \author Mathieu Pont (mathieu.pont@lip6.fr)
/// \date 2021.
///
/// Utils function for manipulating FTMTree class

#pragma once

#include <FTMTree.h>
#include <FTMTreeUtils.h>

#include <ttkUtils.h>

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>

namespace ttk {
  namespace ftm {

    template <class dataType>
    MergeTree<dataType> makeTree(vtkUnstructuredGrid *treeNodes,
                                 vtkUnstructuredGrid *treeArcs) {
      auto treeNodeIdArray = treeNodes->GetPointData()->GetArray("TreeNodeId");

      // Init Scalars
      Scalars scalars;
      vtkSmartPointer<vtkDataArray> nodesScalar
        = treeNodes->GetPointData()->GetArray("Scalar"); // 1: Scalar
      scalars.size = nodesScalar->GetNumberOfTuples();
      std::vector<dataType> scalarsValues(nodesScalar->GetNumberOfTuples());
      for(int i = 0; i < nodesScalar->GetNumberOfTuples(); ++i) {
        int index = (treeNodeIdArray ? treeNodeIdArray->GetTuple1(i) : i);
        scalarsValues[index] = nodesScalar->GetTuple1(i);
      }
      scalars.values = (void *)(scalarsValues.data());

      // Init Tree
      Params params;
      params.treeType = Join_Split;
      MergeTree<dataType> mergeTree(scalars, scalarsValues, params);

      // Add Nodes
      vtkSmartPointer<vtkDataArray> nodesId
        = treeNodes->GetPointData()->GetArray("NodeId"); // 0: NodeId
      vtkIdType nodesNumTuples = nodesId->GetNumberOfTuples();
      for(vtkIdType i = 0; i < nodesNumTuples; ++i) {
        mergeTree.tree.makeNode(i);
      }

      // Add Arcs
      vtkSmartPointer<vtkDataArray> arcsUp
        = treeArcs->GetCellData()->GetArray("upNodeId"); // 1: upNodeId
      vtkSmartPointer<vtkDataArray> arcsDown
        = treeArcs->GetCellData()->GetArray("downNodeId"); // 2: downNodeId
      vtkIdType arcsNumTuples = arcsUp->GetNumberOfTuples();
      vtkSmartPointer<vtkDataArray> dummyArcArray
        = treeArcs->GetCellData()->GetArray("isDummyArc");
      std::set<std::tuple<double, double>> added_arcs; // Avoid duplicates
      for(vtkIdType i = 0; i < arcsNumTuples; ++i) {
        if(dummyArcArray != nullptr and dummyArcArray->GetTuple1(i) == 1)
          continue;
        double downId = arcsDown->GetTuple1(i);
        double upId = arcsUp->GetTuple1(i);
        if(treeNodeIdArray) {
          downId = treeNodeIdArray->GetTuple1(downId);
          upId = treeNodeIdArray->GetTuple1(upId);
        }
        auto it = added_arcs.find(std::make_tuple(downId, upId));
        if(it == added_arcs.end()) { // arc not added yet
          mergeTree.tree.makeSuperArc(downId, upId); // (down, Up)
          added_arcs.insert(std::make_tuple(downId, upId));
        }
      }

      // Manage inconsistent arcs
      manageInconsistentArcsMultiParent(&(mergeTree.tree));

      // Remove self link
      removeSelfLink(&(mergeTree.tree));

      return mergeTree;
    }

    // Take a branch decomposition tree in input
    template <class dataType>
    void makePDGridFromBDTree(ftm::MergeTree<dataType> &barycenter,
                              vtkUnstructuredGrid *vtkOutputPD) {
      vtkNew<vtkPoints> points{};

      vtkNew<vtkDoubleArray> birthArray{};
      birthArray->SetName("Birth");
      vtkNew<vtkDoubleArray> deathArray{};
      deathArray->SetName("Death");
      vtkNew<vtkIntArray> isMinMaxPairArray{};
      isMinMaxPairArray->SetName("isMinMaxPair");
      vtkNew<vtkIntArray> treeNodeIdArray{};
      treeNodeIdArray->SetName("TreeNodeId");
      vtkNew<vtkIntArray> treeNodeId2Array{};
      treeNodeId2Array->SetName("TreeNodeIdOrigin");

      ftm::FTMTree_MT *tree = &(barycenter.tree);
      std::queue<ftm::idNode> queue;
      queue.emplace(tree->getRoot());
      while(!queue.empty()) {
        ftm::idNode node = queue.front();
        queue.pop();

        auto birthDeath = tree->getBirthDeath<dataType>(node);
        double birth = std::get<0>(birthDeath);
        double death = std::get<1>(birthDeath);
        birthArray->InsertNextTuple1(birth);
        deathArray->InsertNextTuple1(death);

        bool isMinMaxPair = tree->isRoot(node);
        isMinMaxPairArray->InsertNextTuple1(isMinMaxPair);

        /*auto birthDeathNode = tree->getBirthDeathNode<dataType>(node);
        auto birthNode = std::get<0>(birthDeathNode);
        auto deathNode = std::get<1>(birthDeathNode);*/
        treeNodeIdArray->InsertNextTuple1(node);
        treeNodeId2Array->InsertNextTuple1(tree->getNode(node)->getOrigin());

        double point[3] = {birth, death, 0};
        points->InsertNextPoint(point);

        std::vector<ftm::idNode> children;
        tree->getChildren(node, children);
        for(auto child : children)
          queue.emplace(child);
      }

      vtkOutputPD->SetPoints(points);
      vtkOutputPD->GetPointData()->AddArray(birthArray);
      vtkOutputPD->GetPointData()->AddArray(deathArray);
      vtkOutputPD->GetPointData()->AddArray(isMinMaxPairArray);
      vtkOutputPD->GetPointData()->AddArray(treeNodeIdArray);
      vtkOutputPD->GetPointData()->AddArray(treeNodeId2Array);
    }

    // Returns a branch decomposition tree
    template <class dataType>
    MergeTree<dataType>
      makeBDTreeFromPDGrid(vtkUnstructuredGrid *persistenceDiagram) {
      auto birthArray = persistenceDiagram->GetPointData()->GetArray("Birth");
      auto deathArray = persistenceDiagram->GetPointData()->GetArray("Death");
      auto isMinMaxPairArray
        = persistenceDiagram->GetPointData()->GetArray("isMinMaxPair");

      auto treeNodeIdArray
        = persistenceDiagram->GetPointData()->GetArray("TreeNodeId");
      auto treeNodeId2Array
        = persistenceDiagram->GetPointData()->GetArray("TreeNodeIdOrigin");
      bool gotNodeArrays = (treeNodeIdArray and treeNodeId2Array);

      auto noPairs = birthArray->GetNumberOfTuples();
      int noNodes = noPairs * 2;
      if(gotNodeArrays) {
        for(vtkIdType i = 0; i < treeNodeIdArray->GetNumberOfTuples(); ++i) {
          int val = std::max(treeNodeIdArray->GetTuple1(i),
                             treeNodeId2Array->GetTuple1(i))
                    + 1;
          noNodes = std::max(noNodes, val);
        }
      }
      std::vector<dataType> scalarsVector(noNodes);

      // Init Tree
      MergeTree<dataType> mergeTree
        = ttk::ftm::createEmptyMergeTree<dataType>(scalarsVector.size());

      // Init Tree Structure
      int minMaxPairIndex = -1;
      for(vtkIdType i = 0; i < noPairs; ++i)
        if(isMinMaxPairArray->GetTuple1(i) == 1)
          minMaxPairIndex = i;
      for(vtkIdType i = 0; i < noNodes; ++i)
        mergeTree.tree.makeNode(i);
      for(vtkIdType i = 0; i < noPairs; ++i) {
        int index1 = (gotNodeArrays ? treeNodeIdArray->GetTuple1(i) : i * 2);
        int index2
          = (gotNodeArrays ? treeNodeId2Array->GetTuple1(i) : i * 2 + 1);
        mergeTree.tree.getNode(index1)->setOrigin(index2);
        mergeTree.tree.getNode(index2)->setOrigin(index1);
        scalarsVector[index1] = birthArray->GetTuple1(i);
        scalarsVector[index2] = deathArray->GetTuple1(i);

        if(i != minMaxPairIndex) {
          auto up
            = (gotNodeArrays ? treeNodeId2Array->GetTuple1(minMaxPairIndex)
                             : minMaxPairIndex * 2 + 1);
          // mergeTree.tree.makeSuperArc(index1, up);
          mergeTree.tree.makeSuperArc(index2, up);
        }
      }

      // Init scalars
      ttk::ftm::setTreeScalars<dataType>(mergeTree, scalarsVector);

      return mergeTree;
    }

    inline void loadBlocks(std::vector<vtkMultiBlockDataSet *> &inputTrees,
                           vtkMultiBlockDataSet *blocks) {
      if(blocks != nullptr) {
        inputTrees.resize(blocks->GetNumberOfBlocks());
        for(size_t i = 0; i < inputTrees.size(); ++i) {
          inputTrees[i]
            = vtkMultiBlockDataSet::SafeDownCast(blocks->GetBlock(i));
        }
      }
    }

    template <class dataType>
    bool constructTrees(std::vector<vtkMultiBlockDataSet *> &inputTrees,
                        std::vector<MergeTree<dataType>> &intermediateTrees,
                        std::vector<vtkUnstructuredGrid *> &treesNodes,
                        std::vector<vtkUnstructuredGrid *> &treesArcs,
                        std::vector<vtkDataSet *> &treesSegmentation) {
      bool isPersistenceDiagram = false;
      const int numInputs = inputTrees.size();
      intermediateTrees = std::vector<MergeTree<dataType>>(numInputs);
      treesNodes = std::vector<vtkUnstructuredGrid *>(numInputs);
      treesArcs = std::vector<vtkUnstructuredGrid *>(numInputs);
      treesSegmentation = std::vector<vtkDataSet *>(numInputs);
      for(int i = 0; i < numInputs; i++) {
        if(inputTrees[i]->GetNumberOfBlocks() >= 2) {
          treesNodes[i]
            = vtkUnstructuredGrid::SafeDownCast(inputTrees[i]->GetBlock(0));
          treesArcs[i]
            = vtkUnstructuredGrid::SafeDownCast(inputTrees[i]->GetBlock(1));
          if(inputTrees[i]->GetNumberOfBlocks() > 2)
            treesSegmentation[i]
              = vtkDataSet::SafeDownCast(inputTrees[i]->GetBlock(2));
          intermediateTrees[i]
            = makeTree<dataType>(treesNodes[i], treesArcs[i]);
        } else {
          vtkUnstructuredGrid *persistenceDiagram
            = vtkUnstructuredGrid::SafeDownCast(inputTrees[i]->GetBlock(0));
          intermediateTrees[i]
            = makeBDTreeFromPDGrid<dataType>(persistenceDiagram);
          isPersistenceDiagram = true;
        }
      }
      return isPersistenceDiagram;
    }

    template <class dataType>
    bool constructTrees(std::vector<vtkMultiBlockDataSet *> &inputTrees,
                        std::vector<MergeTree<dataType>> &intermediateTrees) {
      std::vector<vtkUnstructuredGrid *> treesNodes;
      std::vector<vtkUnstructuredGrid *> treesArcs;
      std::vector<vtkDataSet *> treesSegmentation;
      return constructTrees(inputTrees, intermediateTrees, treesNodes,
                            treesArcs, treesSegmentation);
    }
  } // namespace ftm
} // namespace ttk
