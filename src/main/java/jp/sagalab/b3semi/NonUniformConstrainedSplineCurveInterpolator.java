package jp.sagalab.b3semi;

import jp.sagalab.jftk.LeastSquares;
import jp.sagalab.jftk.Matrix;
import jp.sagalab.jftk.Point;
import jp.sagalab.jftk.curve.Range;
import jp.sagalab.jftk.curve.SplineCurve;
import jp.sagalab.jftk.fuzzy.FuzzyUtils;

import java.util.Arrays;

/**
 * 制約付きスプライン曲線補間を行うためのクラスです。
 *
 * @author Takumi Inagaki
 */

public final class NonUniformConstrainedSplineCurveInterpolator {

  /**
   * 指定された点列に対して制約付きスプライン曲線補間を行います。
   *
   * @param _points       点列
   * @param _degree       次数
//   * @param _knotInterval 節点間隔
   * @return スプライン曲線
   * @throws IllegalArgumentException 次数が0以下の場合
   * @throws IllegalArgumentException 節点間隔が0以下の場合
   * @throws IllegalArgumentException 点列がnullの場合
   * @throws IllegalArgumentException 点列にnullが含まれる場合
   * @throws IllegalArgumentException 点列の要素数が1以下の場合
   * @throws IllegalArgumentException 点列中の時刻がNaN、もしくは無限大の場合
   * @throws IllegalArgumentException 点列中に時間的に逆行している箇所があった場合
   */
  public static SplineCurve interpolate(Point[] _points, int _degree, double[] _knots) {
    // 次数のチェック
    if (_degree < 1) {
      throw new IllegalArgumentException(" degree is must be greater than 0 ");
    }
    // 節点間隔チェック
//    if (_knotInterval <= 0.0) {
//      throw new IllegalArgumentException(" knot's interval is must be greater than 0 ");
//    }
    if (_points == null) {
      throw new IllegalArgumentException("_points is null.");
    }
    // 入力点列にnullが混入していないかチェック
    if (Arrays.asList(_points).contains(null)) {
      throw new IllegalArgumentException(" points include null ");
    }
    // 点列の要素数チェック
    if (_points.length < 2) {
      throw new IllegalArgumentException(" points's length must be greater than 1 ");
    }

    // 時系列チェック
    double preTime = Double.NEGATIVE_INFINITY;
    boolean isFuzzy = false;
    for (Point p : _points) {
      double t = p.time();
      if (Double.isNaN(t) || Double.isInfinite(t)) {
        throw new IllegalArgumentException("point's time include NaN or infinite");
      }
      if (!(preTime <= t)) {
        throw new IllegalArgumentException("time series is not a positive order");
      }
      if (!isFuzzy) {
        isFuzzy = (p.fuzziness() > 0.0);
      }
      preTime = t;
    }

    Range range = Range.create(_points[0].time(), _points[_points.length - 1].time());

    // 節点系列の生成
//    double[] knots = createKnots(range, _degree, _knotInterval);
    double[] knots = _knots;

    // PointsGraphの生成
    Main.createPointsGraph(_points, knots);

    // 重み行列の生成
    Matrix wmat = createWeightMatrix(_points, _degree, knots);

    // 制御点列の導出
    Point[] controlPoints = calculateControlPoints(wmat, _points, _degree, knots);

    //点列中にファジィ点が含まれていた場合はファジィスプライン曲線補間を行う
    if (isFuzzy) {
      double[] observations = new double[_points.length];
      for (int i = 0; i < observations.length; ++i) {
        observations[i] = _points[i].fuzziness();
      }
      double[] fuzzinessElements = FuzzyUtils.nnls(wmat, observations);
      for (int i = 0; i < controlPoints.length; ++i) {
        controlPoints[i] = Point.createXYZTF(controlPoints[i].x(), controlPoints[i].y(), controlPoints[i].z(),
                controlPoints[i].time(), fuzzinessElements[i]);
      }
    }

    // スプライン曲線構築
    return SplineCurve.create(_degree, controlPoints, knots, range);
  }

  /**
   * 節点系列を生成します。
   *
   * @param _range        存在範囲
   * @param _degree       次数
   * @param _knotInterval 節点間隔
   * @return 節点系列
   */
  private static double[] createKnots(Range _range, int _degree, double _knotInterval) {
    // 節点系列の生成
    double start = _range.start();
    double end = _range.end();
    // 有効定義域の節点区間数
    int knotIntervalNum = (int) Math.round((end - start) / _knotInterval);
    double[] knots = new double[knotIntervalNum + 2 * _degree - 1];

    for (int i = 0; i < knots.length; ++i) {
      double w = (i - _degree + 1) / (double) knotIntervalNum;
      knots[i] = (1.0 - w) * start + w * end;
    }

    return knots;
  }

  /**
   * スプライン曲線の重み行列を生成します。<br>
   * 生成する行列は行数：入力点数、列数：制御点数となります。
   *
   * @param _points 入力点列
   * @param _degree 次数
   * @param _knots  節点系列
   * @return 重み行列
   */
  public static Matrix createWeightMatrix(Point[] _points, int _degree, double[] _knots) {
    // 生成する行列は行数：入力点数、列数：制御点数
    final int pointsNum = _points.length;
    double[][] elements = new double[pointsNum][];

    // 各入力点の時刻での重み列を導出し、重み行列として構成する
    for (int i = 0; i < pointsNum; ++i) {
      // ある時刻における重み列（各制御点に対応する重みの列）の導出
      elements[i] = calculateWeights(_knots, _degree, _points[i].time());
    }

    return Matrix.create(elements);
  }

  /**
   * ある時刻における重み列を導出します。
   *
   * @param _knots  節点系列
   * @param _degree 次数
   * @param _time   時刻
   * @return 重み列
   */
  private static double[] calculateWeights(double[] _knots, int _degree, double _time) {
    // 時刻に対応する節点番号( _knots[ num ] <= _time <= _knots[ num + 1 ] )の取得
    int num = _degree;
    int end = _knots.length - _degree;
    while (num < end && _time > _knots[num]) {
      ++num;
    }

    double[] part = new double[]{1.0};

    for (int i = 1; i <= _degree; ++i) {
      double[] now = new double[i + 1];
      for (int j = 0; j <= i; ++j) {
        double tmp = 0;
        int base = num + j - 1;
        if (j != 0) {
          final double d = _knots[base - i];
          tmp += (_time - d) * part[j - 1] / (_knots[base] - d);
        }
        if (j != i) {
          final double d = _knots[base + 1];
          tmp += (d - _time) * part[j] / (d - _knots[base + 1 - i]);
        }
        now[j] = tmp;
      }
      part = now;
    }

    double[] weights = new double[_knots.length - _degree + 1];
    System.arraycopy(part, 0, weights, num - _degree, _degree + 1);
    return weights;
  }

  /**
   * 制御点列を導出します。
   *
   * @param _mat          重み行列
   * @param _points       通過点列
   * @param _degree       次数
   * @param _knots        節点系列
   //* @param _knotInterval 節点間隔
   * @return 制御点列
   */
  private static Point[] calculateControlPoints(Matrix _mat, Point[] _points, int _degree, double[] _knots) {
    double[][] elements = new double[_points.length][];
    for (int i = 0; i < _points.length; ++i) {
      Point p = _points[i];
      elements[i] = new double[]{p.x(), p.y(), p.z()};
    }
    Matrix elementsMatrix = Matrix.create(elements);

    Matrix result;
    //↓本来_knotIntervalが存在する
    Matrix cMat = createCMat(_points, _degree, _knots);

    if (cMat != null) {
      Matrix qMat = createQMat(cMat.rowSize(), elementsMatrix.columnSize());
      result = LeastSquares.solveConstrained(_mat, elementsMatrix,
              cMat, qMat);
    } else {
      result = LeastSquares.solve(_mat, elementsMatrix);
    }

    // 制御点列の構成
    Point[] controlPoints = new Point[_knots.length - _degree + 1];
    for (int i = 0; i < controlPoints.length; ++i) {
      controlPoints[i] = Point.createXYZT(result.get(i, 0), result.get(i, 1), result.get(i, 2), cpTime(_degree, _knots, i));
    }

    return controlPoints;
  }

  /**
   * グレビル横座標を導出します。
   *
   * @param _degree 次数
   * @param _knots  節点列
   * @param _i      制御点のインデックス
   * @return 制御点のグレビル横座標
   */
  private static double cpTime(int _degree, double[] _knots, int _i) {
    double sum = 0.0;
    for (int j = _i; j <= _i + _degree - 1; j++) {
      sum += _knots[j];
    }

    return sum / _degree;
  }

  /**
   * 目的関数 Cx = d の左辺行列 C を生成します。
   *
   * @param _points       点列
   * @param _degree       次数
   * @param _knots        節点列
   //* @param _knotInterval 節点間隔
   * @return 目的関数の左辺行列 C
   */
  private static Matrix createCMat(Point[] _points, int _degree, double[] _knots) {
    int cpSize = _knots.length - _degree + 1;
    Matrix cMat = null;
    for (int i = 0; i < cpSize; i++) {

      // 定義域外の処理
      if (i < Math.floor(_degree / 2.0) || cpSize - Math.floor(_degree / 2.0) <= i) {
        double[][] constraint = new double[1][cpSize];
        // 始点側の処理
        if (i < Math.floor(_degree / 2.0)) {
          constraint[0][i] = 1;
          constraint[0][i + 1] = -1;
          // 終点側の処理
        } else {
          constraint[0][i - 1] = 1;
          constraint[0][i] = -1;
        }
        if (cMat != null) {
          cMat = Matrix.concatVertical(cMat, Matrix.create(constraint));
        } else {
          cMat = Matrix.create(constraint);
        }

      } else {
        int containPointsNum = 0;
        // i番目の制御点のグレビル横座標を取得
        double cpTime = cpTime(_degree, _knots, i);
        for (Point p : _points) {
          // 有効範囲内にある点列数を数える.
//          if (cpTime - (_knotInterval / 2.0) <= p.time() && p.time() < cpTime + (_knotInterval / 2.0)) {
//            containPointsNum++;
//          }
          for (int x=0; x <= _knots.length-2; x++){
            if (_knots[x] <= p.time() && p.time() < _knots[x+1]) {
              containPointsNum++;
            }
          }
        }
        // 有効区間内の点列数が0の場合(1, -2, 1)制約を付ける.
        if (containPointsNum == 0) {
          double[][] constraint = new double[1][cpSize];
          constraint[0][i - 1] = 1;
          constraint[0][i] = - 2;
          constraint[0][i + 1] = 1;
          if (cMat != null) {
            cMat = Matrix.concatVertical(cMat, Matrix.create(constraint));
          } else {
            cMat = Matrix.create(constraint);
          }
        }
      }
    }

    return cMat;
  }

  /**
   * 目的関数 Cx = d の右辺行列 d を生成します。
   *
   * @param _rowSize    行列の行数
   * @param _columnSize 行列の列数
   * @return 目的関数の右辺行列 d
   */
  private static Matrix createQMat(int _rowSize, int _columnSize) {
    double[][] constraintAns = new double[_rowSize][_columnSize];
    return Matrix.create(constraintAns);
  }

  private NonUniformConstrainedSplineCurveInterpolator() {
    throw new UnsupportedOperationException("can not create instance.");
  }
}
