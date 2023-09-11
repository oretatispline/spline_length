package jp.sagalab.b3semi;

import java.util.Arrays;
import jp.sagalab.jftk.Point;

/**
 * スプライン曲線を表すクラスです。
 *
 * @author Akira Nishikawa
 */
public class SplineCurve {

  /**
   * スプライン曲線を生成します。
   *
   * @param _degree        次数
   * @param _controlPoints 制御点列
   * @param _knots         節点系列
   * @param _range         存在範囲
   * @return スプライン曲線
   * @throws IllegalArgumentException スプライン曲線の次数が1未満の場合
   * @throws IllegalArgumentException 存在範囲の始点が節点系列の(次数 - 1)番目よりも小さい場合、
   *                                  または、存在範囲の終点が節点系列の(節点系列の要素数 - 次数)番目よりも大きい場合
   */
  public static SplineCurve create(int _degree, Point[] _controlPoints, double[] _knots, Range _range) {
    // 次数のチェック
    if (_degree < 1) {
      throw new IllegalArgumentException("_degree < 1");
    }
    // 存在範囲と節点系列の整合性チェック
    if (_range.start() < _knots[_degree - 1] || _knots[_knots.length - _degree] < _range.end()) {
      throw new IllegalArgumentException("There is no consistency of _range and _knots.");
    }
    // 節点系列と制御点列の整合性チェック
    if (_knots.length != _controlPoints.length + _degree - 1) {
      throw new IllegalArgumentException("_knots.length NOT equals (_controlPoints.length + _degree - 1).");
    }
    // 節点系列のチェック
    double pre = _knots[0];
    for (double d : _knots) {
      if (Double.isInfinite(d) || Double.isNaN(d)) {
        throw new IllegalArgumentException("_knots is included in infinity or NaN.");
      }
      if (d < pre) {
        throw new IllegalArgumentException("There are counter flowed _knots.");
      }
      pre = d;
    }
    return new SplineCurve(_degree, _controlPoints, _knots, _range);
  }

  /**
   * 指定されたパラメータでの点を評価します。
   *
   * @param _t パラメータ
   * @return 評価点
   */
  public Point evaluate(double _t) {
    // 対象となる節点番号を求める
    int knotNum = searchKnotNum(_t, m_degree - 1, m_knots.length - m_degree);
    // 部分制御点列の抽出
    Point[] part = new Point[m_degree + 1];
    System.arraycopy(m_cp, knotNum - m_degree, part, 0, part.length);
    // de Boor による評価
    for (int i = 0; i < m_degree; ++i) {
      for (int j = 0; j < m_degree - i; ++j) {
        int k = knotNum - j - 1;
        double w = (_t - m_knots[k]) / (m_knots[k + m_degree - i] - m_knots[k]);
        part[m_degree - j] = part[m_degree - j].internalDivision(part[m_degree - j - 1], 1 - w, w);
      }
    }
    Point p = part[m_degree];

    return Point.createXYT(p.x(), p.y(), _t);
  }

  /**
   * 節点番号の探索を行います。
   *
   * @param _t        パラメータ
   * @param _minIndex 探索範囲の最小節点番号
   * @param _maxIndex 探索範囲の最大節点番号
   * @return 節点番号
   * @throws ArrayIndexOutOfBoundsException 指定された節点番号が不正な値の場合
   * @throws IllegalArgumentException       _minIndex が _maxIndex より大きい場合
   */
  public int searchKnotNum(double _t, int _minIndex, int _maxIndex) {
    if (_minIndex < 0 || _maxIndex >= m_knots.length) {
      throw new ArrayIndexOutOfBoundsException("_minIndex < 0 || _maxIndex >= m_knots.length");
    }
    if (_minIndex > _maxIndex) {
      throw new IllegalArgumentException("_minIndex > _maxIndex");
    }

    if (m_knots[_maxIndex] <= _t) {
      return _maxIndex;
    }
    if (_minIndex < _maxIndex) { // 探索区間がある場合
      // 二分探索法
      do {
        int i = (_minIndex + _maxIndex) / 2;
        if (m_knots[i] <= _t && _t < m_knots[i + 1]) {
          return i + 1;
        } else if (_t < m_knots[i]) {
          _maxIndex = i - 1;
        } else /* if ( _knots[i + 1] <= _t ) */ {
          _minIndex = i + 1;
        }
      } while (_minIndex <= _maxIndex);
    }

    return _minIndex + 1;
  }

  /**
   * 次数を返します。
   *
   * @return 次数
   */
  public int degree() {
    return m_degree;
  }

  /**
   * 制御点列を返します。
   *
   * @return 制御点列
   */
  public Point[] controlPoints() {
    return m_cp.clone();
  }

  /**
   * 節点系列を返します。
   *
   * @return 節点系列
   */
  public double[] knots() {
    return m_knots.clone();
  }

  /**
   * パラメータの範囲を返します。
   *
   * @return パラメータの範囲
   */
  public Range range() {
    return Range.create(m_range.start(), m_range.end());
  }

  /**
   * この SplineCurve と指定された Object が等しいかどうかを比較します。
   *
   * @param obj この SplineCurve と比較される Object
   * @return 指定された Object が、このオブジェクトと
   * 次数、制御点列、節点系列、パラメータ範囲がまったく同じ SplineCurve である限りtrue
   */
  @Override
  public boolean equals(Object obj) {
    if (obj == null) {
      return false;
    }
    if (getClass() != obj.getClass()) {
      return false;
    }
    final SplineCurve other = (SplineCurve) obj;
    if (this.m_degree != other.m_degree) {
      return false;
    }
    if (!Arrays.deepEquals(this.m_cp, other.m_cp)) {
      return false;
    }
    if (!Arrays.equals(this.m_knots, other.m_knots)) {
      return false;
    }
    return this.m_range != null && this.m_range.equals(other.m_range);
  }

  /**
   * この SplineCurve の文字列表現を返します。
   *
   * @return 次数、制御点列、節点系列、パラメータ範囲を表す String
   */
  @Override
  public String toString() {
    return String.format(
            "cp:%s knots:%s degree:%d range:%s", Arrays.deepToString(m_cp),
            Arrays.toString(m_knots), m_degree, m_range.toString());
  }

  /**
   * スプライン曲線の複製を取得します.
   *
   * @return 複製したスプライン曲線
   */
  public SplineCurve copy() {
    return SplineCurve.create(m_degree, m_cp.clone(), m_knots.clone(), m_range.copy());
  }


  /**
   * コンストラクタ
   *
   * @param _degree        次数
   * @param _controlPoints 制御点列
   * @param _knots         節点列
   * @param _range         定義域
   */
  public SplineCurve(int _degree, Point[] _controlPoints, double[] _knots, Range _range) {
    m_degree = _degree;
    m_cp = _controlPoints;
    m_knots = _knots;
    m_range = _range;
  }

  /** 次数 */
  private final int m_degree;
  /** 制御点列 */
  private final Point[] m_cp;
  /** 節点系列 */
  private final double[] m_knots;
  /** パラメータ範囲 */
  private final Range m_range;
}