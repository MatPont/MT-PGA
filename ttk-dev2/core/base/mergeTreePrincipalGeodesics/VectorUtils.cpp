#include <VectorUtils.h>

double distanceL2(std::vector<double> &v1, std::vector<double> &v2) {
  double distance = 0.0;
  for(unsigned int i = 0; i < v1.size(); ++i)
    distance += (v1[i] - v2[i]) * (v1[i] - v2[i]);
  return sqrt(distance);
}

double distanceL2Flatten(std::vector<std::vector<double>> &v1,
                         std::vector<std::vector<double>> &v2) {
  std::vector<double> v1_flatten, v2_flatten;
  flatten(v1, v1_flatten);
  flatten(v2, v2_flatten);
  return distanceL2(v1_flatten, v2_flatten);
}

double scalarProduct(std::vector<double> &v1, std::vector<double> &v2) {
  double res = 0.0;
  for(unsigned int i = 0; i < v1.size(); ++i)
    res += v1[i] * v2[i];
  return res;
}

double scalarProductFlatten(std::vector<std::vector<double>> &v1,
                            std::vector<std::vector<double>> &v2) {
  std::vector<double> v1_flatten, v2_flatten;
  flatten(v1, v1_flatten);
  flatten(v2, v2_flatten);
  return scalarProduct(v1_flatten, v2_flatten);
}

double norm(std::vector<double> &v) {
  return sqrt(scalarProduct(v, v));
}

double normFlatten(std::vector<std::vector<double>> &v) {
  std::vector<double> v_flatten;
  flatten(v, v_flatten);
  return norm(v_flatten);
}

void vectorProjection(std::vector<double> &v1,
                      std::vector<double> &v2,
                      std::vector<double> &projec) {
  projec = std::vector<double>(v1.size(), 0.0);
  double normV = norm(v1);
  double squared_norm = normV * normV;
  if(squared_norm > 1e-12) {
    double projecV2 = scalarProduct(v1, v2);
    projecV2 /= squared_norm;
    for(unsigned int i = 0; i < v1.size(); ++i)
      projec[i] = v1[i] * projecV2;
  }
}

void sumProjection(std::vector<double> &v1,
                   std::vector<double> &v2,
                   std::vector<double> &v1Out,
                   std::vector<double> &v2Out) {
  std::vector<double> sumV;
  sumVector(v1, v2, sumV);
  vectorProjection(sumV, v1, v1Out);
  vectorProjection(sumV, v2, v2Out);
}

void gramSchmidt(std::vector<std::vector<double>> &vS,
                 std::vector<double> &v,
                 std::vector<double> &newV) {
  std::vector<std::vector<double>> allVs = vS, uS;
  allVs.push_back(v);
  uS = allVs;
  for(unsigned int i = 1; i < allVs.size(); ++i) {
    std::vector<double> projecSum;
    vectorProjection(uS[0], allVs[i], projecSum);
    for(unsigned int j = 1; j < i; ++j) {
      std::vector<double> projecTemp, projecSumTemp;
      vectorProjection(uS[j], allVs[i], projecTemp);
      sumVector(projecSum, projecTemp, projecSumTemp);
      projecSum = projecSumTemp;
    }
    subVector(allVs[i], projecSum, uS[i]);
  }
  newV = uS[uS.size() - 1];
}

void sumVector(std::vector<double> &v1,
               std::vector<double> &v2,
               std::vector<double> &sumV) {
  if(sumV.size() != v1.size())
    sumV = std::vector<double>(v1.size());
  for(unsigned int i = 0; i < v1.size(); ++i)
    sumV[i] = v1[i] + v2[i];
}

void subVector(std::vector<double> &v1,
               std::vector<double> &v2,
               std::vector<double> &subV) {
  subV = std::vector<double>(v1.size());
  for(unsigned int i = 0; i < v1.size(); ++i)
    subV[i] = v1[i] - v2[i];
}

void multVectorByScalar(std::vector<double> &v,
                        double scalar,
                        std::vector<double> &multV) {
  multV = v;
  for(unsigned int i = 0; i < multV.size(); ++i)
    multV[i] *= scalar;
}

void multVectorByScalarFlatten(std::vector<std::vector<double>> &v,
                               double scalar,
                               std::vector<std::vector<double>> &multV) {
  std::vector<double> v_flat, multV_flat;
  flatten(v, v_flat);
  multVectorByScalar(v_flat, scalar, multV_flat);
  unflatten(multV_flat, multV);
}

double sumVectorElements(std::vector<double> &v) {
  double sum = 0;
  for(auto e : v)
    sum += e;
  return sum;
}

double sumVectorElementsFlatten(std::vector<std::vector<double>> &v) {
  std::vector<double> v_flatten;
  flatten(v, v_flatten);
  return sumVectorElements(v_flatten);
}

void multiSumVector(std::vector<std::vector<double>> &v1,
                    std::vector<std::vector<double>> &v2,
                    std::vector<std::vector<double>> &sumV) {
  sumV = std::vector<std::vector<double>>(v1.size());
  for(unsigned int i = 0; i < v1.size(); ++i)
    sumVector(v1[i], v2[i], sumV[i]);
}

void multiSumVectorFlatten(std::vector<std::vector<std::vector<double>>> &v1,
                           std::vector<std::vector<std::vector<double>>> &v2,
                           std::vector<std::vector<double>> &sumV) {
  std::vector<std::vector<double>> v1_flatten, v2_flatten;
  multiFlatten(v1, v1_flatten);
  multiFlatten(v2, v2_flatten);
  multiSumVector(v1_flatten, v2_flatten, sumV);
}

void flatten(std::vector<std::vector<double>> &v, std::vector<double> &newV) {
  newV = std::vector<double>(v.size() * v[0].size());
  for(unsigned int i = 0; i < v.size(); ++i)
    for(unsigned int j = 0; j < v[0].size(); ++j)
      newV[i * v[0].size() + j] = v[i][j];
}

void multiFlatten(std::vector<std::vector<std::vector<double>>> &v,
                  std::vector<std::vector<double>> &newV) {
  newV = std::vector<std::vector<double>>(v.size());
  for(unsigned int i = 0; i < v.size(); ++i)
    flatten(v[i], newV[i]);
}

void unflatten(std::vector<double> &v, std::vector<std::vector<double>> &newV) {
  newV = std::vector<std::vector<double>>(v.size() / 2);
  for(unsigned int i = 0; i < v.size(); i += 2)
    newV[i / 2] = std::vector<double>{v[i], v[i + 1]};
}

bool isVectorUniform(std::vector<double> &v) {
  for(unsigned int i = 0; i < v.size() - 1; ++i)
    if(not(std::abs(v[i] - v[i + 1]) < 1e-6))
      return false;
  return true;
}

bool isVectorNull(std::vector<double> &v) {
  for(unsigned int i = 0; i < v.size(); ++i)
    if(not(std::abs(v[i]) < 1e-6))
      return false;
  return true;
}

bool isVectorNullFlatten(std::vector<std::vector<double>> &v) {
  std::vector<double> v_flat;
  flatten(v, v_flat);
  return isVectorNull(v_flat);
}

// Matrix utils

void matrixDot(std::vector<std::vector<double>> &m1,
               std::vector<std::vector<double>> &m2,
               std::vector<std::vector<double>> &newM) {
  newM = std::vector<std::vector<double>>(
    m1.size(), std::vector<double>(m2[0].size(), 0.0));
  for(unsigned int i = 0; i < newM.size(); ++i)
    for(unsigned int j = 0; j < newM[i].size(); ++j)
      for(unsigned int k = 0; k < m1[i].size(); ++k)
        newM[i][j] += m1[i][k] * m2[k][j];
}

void subMatrix(std::vector<std::vector<double>> &m1,
               std::vector<std::vector<double>> &m2,
               std::vector<std::vector<double>> &newM) {
  newM = std::vector<std::vector<double>>(
    m1.size(), std::vector<double>(m1[0].size()));
  for(unsigned int i = 0; i < m1.size(); ++i)
    for(unsigned int j = 0; j < m1[0].size(); ++j)
      newM[i][j] = m1[i][j] - m2[i][j];
}

void sumMatrix(std::vector<std::vector<double>> &m1,
               std::vector<std::vector<double>> &m2,
               std::vector<std::vector<double>> &newM) {
  newM = std::vector<std::vector<double>>(
    m1.size(), std::vector<double>(m1[0].size()));
  for(unsigned int i = 0; i < m1.size(); ++i)
    for(unsigned int j = 0; j < m1[0].size(); ++j)
      newM[i][j] = m1[i][j] + m2[i][j];
}

void multMatrix(std::vector<std::vector<double>> &m1,
                double mult,
                std::vector<std::vector<double>> &newM) {
  newM = m1;
  for(unsigned int i = 0; i < newM.size(); ++i)
    for(unsigned int j = 0; j < newM[i].size(); ++j)
      newM[i][j] *= mult;
}

void transpose(std::vector<std::vector<double>> &m,
               std::vector<std::vector<double>> &newM) {
  newM = std::vector<std::vector<double>>(
    m[0].size(), std::vector<double>(m.size()));
  for(unsigned int i = 0; i < m.size(); ++i)
    for(unsigned int j = 0; j < m[0].size(); ++j)
      newM[j][i] = m[i][j];
}

// Statistics utils

double mean(std::vector<double> &v) {
  double mean = 0.0;
  for(auto e : v)
    mean += e;
  return mean / v.size();
}

double var(std::vector<double> &v) {
  return cov(v, v);
}

double cov(std::vector<double> &v1, std::vector<double> &v2) {
  double cov = 0.0;
  double meanV1 = mean(v1);
  double meanV2 = mean(v2);
  for(unsigned int i = 0; i < v1.size(); ++i)
    cov += (v1[i] - meanV1) * (v2[i] - meanV2);
  return cov / (v1.size() - 1);
}

double corr(std::vector<double> &v1, std::vector<double> &v2) {
  return cov(v1, v2) / (std::sqrt(var(v1)) * std::sqrt(var(v2)));
}
