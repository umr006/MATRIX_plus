#include <math.h>
#include <cstring>
#include <iostream>

class S21Matrix {
 public:
  S21Matrix();
  S21Matrix(int rows, int cols);
  S21Matrix(const S21Matrix& other);
  S21Matrix(S21Matrix &&other) noexcept;
  ~S21Matrix();

  bool EqMatrix(const S21Matrix& o) const noexcept;
  void SumMatrix(const S21Matrix& other);
  void SubMatrix(const S21Matrix& other);
  void MulMatrix(const S21Matrix& other);
  void MulNumber(const double num) noexcept;
  S21Matrix Transpose() noexcept;
  S21Matrix FindMinor(int i, int j);
  double Determinant();
  S21Matrix CalcComplements();
  S21Matrix InverseMatrix();

  void FreeMem() noexcept;
  void AllocMem(int rows, int cols) noexcept;
  void MatrixCP(const S21Matrix& other) noexcept;
  void CheckMatrix(const S21Matrix& other);

  int GetRows() const;
  int GetCols() const;
  double** GetMatrix() const;
  double GetMatrixElement(int i, int j);
  void SetRows(int rows);
  void SetCols(int cols);
  void SetMatrixElement(int i, int j, double value);

  S21Matrix& operator=(const S21Matrix& other) noexcept;
  S21Matrix& operator=(S21Matrix&& other) noexcept;
  S21Matrix operator+(const S21Matrix& other);
  S21Matrix operator-(const S21Matrix& other);
  S21Matrix operator*(const S21Matrix& other);
  friend S21Matrix operator*(const S21Matrix& other, double num);
  bool operator==(const S21Matrix& other) const noexcept;
  S21Matrix& operator+=(const S21Matrix& other);
  S21Matrix& operator-=(const S21Matrix& other);
  S21Matrix& operator*=(const S21Matrix& other);
  S21Matrix& operator*=(const double num);
  double operator()(int i, int j) const;
  double& operator()(int i, int j);

 private:
  int rows_, cols_;
  double** matrix_;
};