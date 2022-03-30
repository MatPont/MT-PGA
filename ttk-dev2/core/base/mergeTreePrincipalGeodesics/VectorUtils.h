/// \ingroup base
/// \class ttk::VectorUtils
/// \author XXX
/// \date 2021
///
/// This module defines the TODO
///

#pragma once
// other includes
#include <math.h>
#include <vector>

// ----------------------------------------------------------------------------
// Vector Utils
// ----------------------------------------------------------------------------

double distanceL2(std::vector<double> &v1, std::vector<double> &v2);

double distanceL2Flatten(std::vector<std::vector<double>> &v1,
                         std::vector<std::vector<double>> &v2);

double scalarProduct(std::vector<double> &v1, std::vector<double> &v2);

double scalarProductFlatten(std::vector<std::vector<double>> &v1,
                            std::vector<std::vector<double>> &v2);

double norm(std::vector<double> &v);

double normFlatten(std::vector<std::vector<double>> &v);

// Project v2 on v1
void vectorProjection(std::vector<double> &v1,
                      std::vector<double> &v2,
                      std::vector<double> &projec);

void sumProjection(std::vector<double> &v1,
                   std::vector<double> &v2,
                   std::vector<double> &v1Out,
                   std::vector<double> &v2Out);

void gramSchmidt(std::vector<std::vector<double>> &vS,
                 std::vector<double> &v,
                 std::vector<double> &newV);

void sumVector(std::vector<double> &v1,
               std::vector<double> &v2,
               std::vector<double> &sumV);

void subVector(std::vector<double> &v1,
               std::vector<double> &v2,
               std::vector<double> &subV);

void multVectorByScalar(std::vector<double> &v,
                        double scalar,
                        std::vector<double> &multV);

void multVectorByScalarFlatten(std::vector<std::vector<double>> &v,
                               double scalar,
                               std::vector<std::vector<double>> &multV);

double sumVectorElements(std::vector<double> &v);

double sumVectorElementsFlatten(std::vector<std::vector<double>> &v);

void multiSumVector(std::vector<std::vector<double>> &v1,
                    std::vector<std::vector<double>> &v2,
                    std::vector<std::vector<double>> &sumV);

void multiSumVectorFlatten(std::vector<std::vector<std::vector<double>>> &v1,
                           std::vector<std::vector<std::vector<double>>> &v2,
                           std::vector<std::vector<double>> &sumV);

void flatten(std::vector<std::vector<double>> &v, std::vector<double> &newV);

void multiFlatten(std::vector<std::vector<std::vector<double>>> &v,
                  std::vector<std::vector<double>> &newV);

void unflatten(std::vector<double> &v, std::vector<std::vector<double>> &newV);

bool isVectorUniform(std::vector<double> &v);

bool isVectorNull(std::vector<double> &v);

bool isVectorNullFlatten(std::vector<std::vector<double>> &v);

// ----------------------------------------------------------------------------
// Matrix Utils
// ----------------------------------------------------------------------------

void matrixDot(std::vector<std::vector<double>> &m1,
               std::vector<std::vector<double>> &m2,
               std::vector<std::vector<double>> &newM);

void subMatrix(std::vector<std::vector<double>> &m1,
               std::vector<std::vector<double>> &m2,
               std::vector<std::vector<double>> &newM);

void sumMatrix(std::vector<std::vector<double>> &m1,
               std::vector<std::vector<double>> &m2,
               std::vector<std::vector<double>> &newM);

void multMatrix(std::vector<std::vector<double>> &m1,
                double mult,
                std::vector<std::vector<double>> &newM);

void transpose(std::vector<std::vector<double>> &m,
               std::vector<std::vector<double>> &newM);

// ----------------------------------------------------------------------------
// Statistics Utils
// ----------------------------------------------------------------------------

double mean(std::vector<double> &v);

double var(std::vector<double> &v);

double cov(std::vector<double> &v1, std::vector<double> &v2);

double corr(std::vector<double> &v1, std::vector<double> &v2);
