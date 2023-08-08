package jp.sagalab.b3semi;

import jp.sagalab.b3semi.graph.PointsGraph;

import java.util.Arrays;

/**
 * スプライン曲線補間を行うためのクラスです。
 *
 * @author Akira Nishikawa
 */
public final class SplineCurveInterpolator {

  /**
   * 指定された点列に対してスプライン曲線補間を行います。
   *
   * @param _points       点列
   * @param _degree       次数
   * @param _knotInterval 節点間隔
   * @return スプライン曲線
   * @throws IllegalArgumentException 次数が0以下の場合
   * @throws IllegalArgumentException 節点間隔が0以下の場合
   * @throws IllegalArgumentException 点列がnullの場合
   * @throws IllegalArgumentException 点列にnullが含まれる場合
   * @throws IllegalArgumentException 点列の要素数が1以下の場合
   * @throws IllegalArgumentException 点列中の時刻がNaN、もしくは無限大の場合
   * @throws IllegalArgumentException 点列中に時間的に逆行している箇所があった場合
   */
  public static SplineCurve interpolate(Point[] _points, int _degree, double _knotInterval) {
    // 次数のチェック
    if (_degree < 1) {
      throw new IllegalArgumentException(" degree is must be greater than 0 ");
    }
    // 節点間隔チェック
    if (_knotInterval <= 0.0) {
      throw new IllegalArgumentException(" knot's interval is must be greater than 0 ");
    }
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
    for (Point p : _points) {
      double t = p.time();
      if (Double.isNaN(t) || Double.isInfinite(t)) {
        throw new IllegalArgumentException("point's time include NaN or infinite");
      }
      if (t <= preTime) {
        throw new IllegalArgumentException("time series is not a positive order");
      }
    }

    Range range = Range.create(_points[0].time(), _points[_points.length - 1].time());

    // 節点系列の生成
    double[] knots = createKnots(range, _degree, _knotInterval);

    // PointsGraphの生成
    Main.createPointsGraph(_points, knots);

    // 重み行列の生成
    Matrix wmat = createWeightMatrix(_points, _degree, knots);
    System.out.println(wmat);

    // 制御点列の導出
    Point[] controlPoints = calculateControlPoints(wmat, _points);

    // スプライン曲線構築
    return SplineCurve.create(_degree, controlPoints, knots, range);
  }

  /**
   * 指定された点列に対してスプライン曲線補間を行います。
   *
   * @param _points       点列
   * @param _degree       次数
   * @param _knots        節点
   * @return スプライン曲線
   * @throws IllegalArgumentException 次数が0以下の場合
   * @throws IllegalArgumentException 節点列がnullの場合
   * @throws IllegalArgumentException 節点列の要素数が2*degree-1以下の場合
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
    // 節点列のチェック
    if (_knots == null) {
      throw new IllegalArgumentException(" _knots is null");
    }
    // 節点列の要素数チェック
    if (_knots.length <= 2 * _degree - 1) {
      throw new IllegalArgumentException(" knots's length must be greater than 2 * _degree - 1 ");
    }
    //入力点列のチェック
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
    for (Point p : _points) {
      double t = p.time();
      if (Double.isNaN(t) || Double.isInfinite(t)) {
        throw new IllegalArgumentException("point's time include NaN or infinite");
      }
      if (t <= preTime) {
        throw new IllegalArgumentException("time series is not a positive order");
      }
    }

    Range range = Range.create(_points[0].time(), _points[_points.length - 1].time());

    // PointsGraphの生成
    Main.createPointsGraph(_points, _knots);

    // 重み行列の生成
    Matrix wmat = createWeightMatrix(_points, _degree, _knots);
    System.out.println(wmat);

    // 制御点列の導出
    Point[] controlPoints = calculateControlPoints(wmat, _points);

    // スプライン曲線構築
    return SplineCurve.create(_degree, controlPoints, _knots, range);
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
    int knotIntervalNum = (int) Math.ceil((end - start) / _knotInterval);
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
   * @param _knots 節点系列
   * @param _time  時刻
   * @return 重み列
   */
  private static double[] calculateWeights(double[] _knots, int _degree, double _time) {
    int size = _knots.length - _degree + 1;

    double[] weights = new double[size];
    for (int i = 0; i < size; ++i) {
      weights[i] = calculateWeight(_knots, _degree, i, _time);
    }

    return weights;
  }

  /**
   * ある時刻における重みを導出します。
   *
   * @param _knots  節点系列
   * @param _degree 次数
   * @param _i      制御点のインデックス
   * @param _time   時刻
   * @return 重み
   */
  private static double calculateWeight(double[] _knots, int _degree, int _i, double _time)
  {
    { // 特別な場合の処理
      int knotsSize = _knots.length;

      // 左端ではブレンドの左側を考慮しない（p.48 図3.14、p.50 図3.16、図3.17）
      if (_i == 0) {
        double coeff = (_knots[_i + _degree] - _time) / (_knots[_i + _degree] - _knots[_i]);
        return coeff * calculateWeight(_knots, _degree- 1, _i + 1, _time);
      }

      // 右端ではブレンドの右側を考慮しない（p.48 図3.14、p.50 図3.16、図3.17）
      if (_i == knotsSize - _degree) {
        double coeff = (_time - _knots[_i - 1]) / (_knots[_i + _degree- 1] - _knots[_i - 1]);
        return coeff * calculateWeight(_knots, _degree- 1, _i, _time);
      }

      if (_degree== 0) {
        return (_knots[_i - 1] <= _time && _time < _knots[_i]) ? 1.0 : 0.0;
      }
    }

    // 通常処理
    // 分母を先に計算して0になったら（0除算が発生しそうなら）その項の係数は0とする
    double denom1 = _knots[_i + _degree- 1] - _knots[_i - 1];
    double denom2 = _knots[_i + _degree] - _knots[_i];
    double coeff1 = (denom1 != 0.0) ? (_time - _knots[_i - 1]) / denom1 : 0.0;
    double coeff2 = (denom2 != 0.0) ? (_knots[_i + _degree] - _time) / denom2 : 0.0;

    return coeff1 * calculateWeight(_knots, _degree- 1, _i, _time)
            + coeff2 * calculateWeight(_knots, _degree- 1, _i + 1, _time);
  }

  /**
   * 制御点列を導出します。
   *
   * @param _mat    重み行列
   * @param _points 通過点列
   * @return 制御点列
   */
  private static Point[] calculateControlPoints(Matrix _mat, Point[] _points) {
    double[][] elements = new double[_points.length][];
    for (int i = 0; i < _points.length; ++i) {
      Point p = _points[i];
      elements[i] = new double[]{p.x(), p.y()};
    }

    Matrix result = LeastSquares.solve(_mat, Matrix.create(elements));

    // 制御点列の構成
    Point[] controlPoints = new Point[result.rowSize()];
    for (int i = 0; i < controlPoints.length; ++i) {
      controlPoints[i] = Point.createXY(result.get(i, 0), result.get(i, 1));
    }

    return controlPoints;
  }

  private SplineCurveInterpolator() {
    throw new UnsupportedOperationException("can not create instance.");
  }
}