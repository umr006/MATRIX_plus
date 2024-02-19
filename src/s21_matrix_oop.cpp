#include "s21_matrix_cpp.h"

S21Matrix::S21Matrix() {
  rows_ = cols_ = 0;
  matrix_ = nullptr;
}

S21Matrix::S21Matrix(int rows, int cols) : rows_(rows), cols_(cols) {
  if (rows < 0 || cols < 0) {
    rows_ = cols_ = 0;
    matrix_ = nullptr;
  } else {
    AllocMem(rows, cols);
  }
}

S21Matrix::S21Matrix(const S21Matrix &other)
    : rows_(other.rows_), cols_(other.cols_) {
      AllocMem(rows_, cols_);
      MatrixCP(other);
  }

S21Matrix::S21Matrix(S21Matrix &&other) noexcept : rows_(0), cols_(0), matrix_(nullptr) {
  std::swap(other.rows_, rows_);
  std::swap(other.cols_, cols_);
  std::swap(other.matrix_, matrix_);
}

S21Matrix::~S21Matrix() {
  FreeMem();
  rows_ = cols_ = 0;
}

bool S21Matrix::EqMatrix(const S21Matrix &o) const noexcept {
  bool return_code = true;
  if (rows_ != o.rows_ || cols_ != o.cols_ || rows_ <= 0 || o.rows_ <= 0 ||
      cols_ <= 0 || o.cols_ <= 0) {
    return_code = false;
  } else {
    for (int i = 0; i < rows_; i++) {
      for (int j = 0; j < cols_; j++) {
        if (fabs(matrix_[i][j]) != fabs(o.matrix_[i][j])) {
          return_code = false;
          break;
        }
        if (!return_code) break;
      }
    }
  }
  return return_code;
}

void S21Matrix::SumMatrix(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::length_error("различная размерность матриц");
  }
  CheckMatrix(*this);
  CheckMatrix(other);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] += other.matrix_[i][j];
    }
  }
}

void S21Matrix::SubMatrix(const S21Matrix &other) {
  if (rows_ != other.rows_ || cols_ != other.cols_) {
    throw std::length_error("различная размерность матриц");
  }
  CheckMatrix(*this);
  CheckMatrix(other);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] -= other.matrix_[i][j];
    }
  }
}

void S21Matrix::MulMatrix(const S21Matrix &other) {
  if (cols_ != other.rows_)
    throw std::length_error(
        "число столбцов первой матрицы не равно числу строк второй матрицы");
  // CheckMatrix(*this);
  // CheckMatrix(other);
  S21Matrix res(rows_, cols_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < other.cols_; j++) {
      for (int k = 0; k < cols_; k++) {
        res.matrix_[i][j] += matrix_[i][k] * other.matrix_[k][j];
      }
    }
  }
  *this = res;
}

void S21Matrix::MulNumber(const double num) noexcept {
  CheckMatrix(*this);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] *= num;
    }
  }
}

S21Matrix S21Matrix::Transpose() noexcept {
  CheckMatrix(*this);

  S21Matrix result_matrix(cols_, rows_);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      result_matrix.matrix_[j][i] = matrix_[i][j];
    }
  }
  return result_matrix;
}

S21Matrix S21Matrix::FindMinor(int i, int j) {
  S21Matrix res(rows_ - 1, cols_ - 1);
  int cnt_rows = 0;
  int cnt_cols = 0;
  for (int k = 0; k < rows_; k++) {
    if (k == i) continue;
    for (int l = 0; l < cols_; l++) {
      if (j == l) continue;
      res.matrix_[cnt_rows][cnt_cols] = matrix_[k][l];
      cnt_cols++;
    }
    cnt_rows++;
    cnt_cols = 0;
  }
  return res;
}

double S21Matrix::Determinant() {
  if (rows_ != cols_) {
    throw std::length_error("матрица не является квадратной");
  }
  CheckMatrix(*this);
  double determinant = 0;
  if (rows_ == 1) {
    determinant = matrix_[0][0];
  } else if (rows_ == 2) {
    determinant = matrix_[0][0] * matrix_[1][1] - matrix_[0][1] * matrix_[1][0];
  } else {
    for (int i = 0; i < cols_; i++) {
      S21Matrix minor = FindMinor(0, i);
      double tmp_determinant = minor.Determinant();
      double tmp_res = pow(-1, 2 + i) * matrix_[0][i] * tmp_determinant;
      determinant += tmp_res;
    }
  }
  return determinant;
}

S21Matrix S21Matrix::CalcComplements() {
  S21Matrix res(rows_, cols_);
  if (rows_ != cols_) {
    throw std::length_error("матрица не является квадратной");
  }
  CheckMatrix(*this);
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      S21Matrix minor = FindMinor(i, j);
      double tmp_determinant = minor.Determinant();
      tmp_determinant *= pow(-1, 2 + i + j);
      res.matrix_[i][j] = tmp_determinant;
    }
  }
  return res;
}

S21Matrix S21Matrix::InverseMatrix() {
  double det = Determinant();
  if (fabs(det) < 1e-10)
    throw std::length_error("определитель матрицы равен 0");
  CheckMatrix(*this);

  S21Matrix res(rows_, cols_);
  S21Matrix tmp(rows_, cols_);
  tmp = CalcComplements();
  res = tmp.Transpose();
  res.MulNumber(1. / fabs(det));
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < rows_; j++) {
      res.matrix_[i][j] *= -1;
    }
  }
  return res;
}

int S21Matrix::GetRows() const { return rows_; }

int S21Matrix::GetCols() const { return cols_; }


void S21Matrix::SetRows(int rows) {
  if (rows < 0)
    throw std::length_error("индекс матрицы не может быть меньше 0");
  if (rows_ != rows) {
    S21Matrix temp(rows, cols_);
    for (int i = 0; i < temp.rows_; i++) {
      for (int j = 0; j < temp.cols_; j++) {
        if (i < rows_) temp.matrix_[i][j] = matrix_[i][j];
      }
    }
    FreeMem();
    AllocMem(temp.rows_, temp.cols_);
    MatrixCP(temp);
  }
}

void S21Matrix::SetCols(int cols) {
  if (cols < 0)
    throw std::length_error("индекс матрицы не может быть меньше 0");
  if (cols_ != cols) {
    S21Matrix temp(rows_, cols);
    for (int i = 0; i < temp.rows_; i++) {
      for (int j = 0; j < temp.cols_; j++) {
        if (j < cols_) temp.matrix_[i][j] = matrix_[i][j];
      }
    }
    FreeMem();
    AllocMem(temp.rows_, temp.cols_);
    MatrixCP(temp);
  }
}

void S21Matrix::SetMatrixElement(int i, int j, double value) {
  if (i < 0 || j < 0)
    throw std::length_error("индекс матрицы не может быть меньше 0");
  matrix_[i][j] = value;
}

double S21Matrix::GetMatrixElement(int i, int j) { return matrix_[i][j]; }

S21Matrix &S21Matrix::operator=(const S21Matrix &other) noexcept {
  if (this != &other) {
    FreeMem();
    AllocMem(other.rows_, other.cols_);
    MatrixCP(other);
  }
  return *this;
}

S21Matrix &S21Matrix::operator=(S21Matrix &&other) noexcept {
  if (this != &other) {
    std::swap(cols_, other.cols_);
    std::swap(rows_, other.rows_);
    std::swap(matrix_, other.matrix_);
  }
  return *this;
}

S21Matrix S21Matrix::operator+(const S21Matrix &other) {
  S21Matrix res = *this;
  res.SumMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator-(const S21Matrix &other) {
  S21Matrix res = *this;
  res.SubMatrix(other);
  return res;
}

S21Matrix S21Matrix::operator*(const S21Matrix &other) {
  S21Matrix res = *this;
  res.MulMatrix(other);
  return res;
}

S21Matrix operator*(const S21Matrix &other, double num) {
  S21Matrix res = other;
  res.MulNumber(num);
  return res;
}

bool S21Matrix::operator==(const S21Matrix &other) const noexcept { return EqMatrix(other); }

S21Matrix &S21Matrix::operator+=(const S21Matrix &other) {
  SumMatrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator-=(const S21Matrix &other) {
  SubMatrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator*=(const S21Matrix &other) {
  MulMatrix(other);
  return *this;
}

S21Matrix &S21Matrix::operator*=(const double num) {
  MulNumber(num);
  return *this;
}

double S21Matrix::operator()(int i, int j) const {
  if (i < 0 || i >= rows_ || j < 0 || j >= cols_)
    throw std::length_error("индекс за пределами матрицы");

  return matrix_[i][j];
}

double &S21Matrix::operator()(int i, int j) {
  if (i < 0 || i >= rows_ || j < 0 || j >= cols_)
    throw std::length_error("индекс за пределами матрицы");

  return matrix_[i][j];
}

void S21Matrix::FreeMem() noexcept {
  if (matrix_ != nullptr) {
    for (int i = 0; i < rows_; i++) {
      delete[] matrix_[i];
    }
    delete[] matrix_;
  }
}

void S21Matrix::AllocMem(int rows, int cols) noexcept {
  rows_ = rows;
  cols_ = cols;
  matrix_ = new double *[rows_];
  for (int i = 0; i < rows_; i++) {
    matrix_[i] = new double[cols_];
  }
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = 0;
    }
  }
}

void S21Matrix::MatrixCP(const S21Matrix &other) noexcept {
  cols_ = other.cols_;
  rows_ = other.rows_;
  for (int i = 0; i < rows_; i++) {
    for (int j = 0; j < cols_; j++) {
      matrix_[i][j] = other(i, j);
    }
  }
}

void S21Matrix::CheckMatrix(const S21Matrix &other) {
  if (other.rows_ < 0 || other.cols_ < 0 || matrix_ == nullptr) {
    throw std::length_error("матрица некорректна");
  }
}