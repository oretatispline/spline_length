//package jp.sagalab.b3semi;
//
///**
// * 平面上の点を表す. 平面上の点のx座標とy座標,時刻を保持する.
// */
//public class Point {
//
//  /**
//   * 指定した座標にある点の生成を行う.
//   * @param _x x座標
//   * @param _y y座標
//   * @return 点
//   */
//  public static Point createXY(double _x, double _y) {
//    return Point.createXYT(_x, _y, 0.0);
//  }
//
//  /**
//   * 指定した座標と時刻を持った点の生成を行う.
//   * @param _x x座標
//   * @param _y y座標
//   * @param _t 時刻
//   * @return 点
//   */
//  public static Point createXYT(double _x, double _y, double _t) {
//    if (Double.isNaN(_x)) {
//      throw new IllegalArgumentException("_x is NaN.");
//    }
//    if (Double.isInfinite(_x)) {
//      throw new IllegalArgumentException("_x is Inf.");
//    }if (Double.isNaN(_y)) {
//      throw new IllegalArgumentException("_y is NaN.");
//    }
//    if (Double.isInfinite(_y)) {
//      throw new IllegalArgumentException("_y is Inf.");
//    }if (Double.isNaN(_t)) {
//      throw new IllegalArgumentException("_t is NaN.");
//    }
//    if (Double.isInfinite(_t)) {
//      throw new IllegalArgumentException("_t is Inf.");
//    }
//
//    return new Point(_x, _y, _t);
//  }
//
//  /**
//   * 指定された点と内分計算を行います。
//   * 自身と他方の点を _tA : _tB に内分する点を求めます.
//   * また、内部で自動的に内分比の合計が１になるように正規化します.
//   * ex) 「1 : 4」 は 「0.2 : 0.8」 となります.
//   * また、内分比の合計が0になる場合は自身の点を返します。
//   *
//   * @param _other 点
//   * @param _tA    内分比A
//   * @param _tB    内分比B
//   * @return 内分点
//   * @throws IllegalArgumentException 内分比にNaNまたはInfiniteが指定された場合
//   */
//  public Point internalDivision(Point _other, double _tA, double _tB) {
//    // 内分比のチェック
//    // 非数のチェック
//    if (Double.isNaN(_tA) || Double.isNaN(_tB)) {
//      throw new IllegalArgumentException("Internal ratio is NaN.");
//    }
//    // 無限大のチェックと補正
//    if (Double.isInfinite(_tA) || Double.isInfinite(_tB)) {
//      throw new IllegalArgumentException("Internal ratio is Infinite.");
//    }
//
//    // 内分比の合計
//    double sumOfRatio = _tA + _tB;
//
//    // 内分計算
//    double x = (_tB * m_x + _tA * _other.m_x) / sumOfRatio;
//    double y = (_tB * m_y + _tA * _other.m_y) / sumOfRatio;
//    double time = (_tB * m_t + _tA * _other.m_t) / sumOfRatio;
//
//    if (Double.isInfinite(x) || Double.isNaN(x)
//            || Double.isInfinite(y) || Double.isNaN(y)) {
//      return this;
//    }
//
//    return createXYT(x, y, time);
//  }
//
//  /**
//   * 点のx座標を取得する.
//   * @return x座標
//   */
//  public double x() {
//    return m_x;
//  }
//
//  /**
//   * 点のy座標を取得する.
//   * @return y座標
//   */
//  public double y() {
//    return m_y;
//  }
//
//  /**
//   * 点の時刻を取得する.
//   * @return 時刻
//   */
//  public double time() {
//    return m_t;
//  }
//
//  /**
//   * この Point の文字列表現を返します。
//   *
//   * @return x座標、y座標、時刻を表す String
//   */
//  @Override
//  public String toString() {
//    return String.format("x:%.3f y:%.3f t:%.3f", m_x, m_y, m_t);
//  }
//
//
//  /**
//   * コンストラクタ
//   *
//   * @param _x x座標
//   * @param _y y座標
//   * @param _t 時刻パラメータ
//   */
//  public Point(double _x, double _y, double _t) {
//    m_x = _x;
//    m_y = _y;
//    m_t = _t;
//  }
//
//  /** x座標 */
//  private final double m_x;
//  /** y座標 */
//  private final double m_y;
//  /** 時刻 */
//  private final double m_t;
//}