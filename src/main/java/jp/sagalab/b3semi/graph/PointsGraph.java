package jp.sagalab.b3semi.graph;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.Line2D;

import jp.sagalab.jftk.Point;

/**
 * 点列と節点列を表すグラフを作成するクラスです.
 *
 * @author inagaki
 */
/*
 適宜改良して使ってください.
 不具合があれば教えてください.
 */
public class PointsGraph extends JPanel {
  /**
   * 点列と節点列を表すグラフを生成します.
   *
   * @param _points 点列
   * @param _knots  節点列
   * @return 点列と節点列を表すグラフ
   * @throws IllegalArgumentException 点列がnullの場合
   * @throws IllegalArgumentException WIDTHが0以下の場合
   * @throws IllegalArgumentException HEIGHTが0以下の場合
   */
  public static PointsGraph create(Point[] _points, double[] _knots) {
    if (_points == null) {
      throw new IllegalArgumentException("_points is null");
    }
    if (WIDTH <= 0) {
      throw new IllegalArgumentException("WIDTH is must be greater than 0");
    }
    if (HEIGHT <= 0) {
      throw new IllegalArgumentException("HEIGHT is must be greater than 0");
    }
    return new PointsGraph(_points, _knots);
  }

  /**
   * グラフの描画を行います.
   *
   * @param g  the <code>Graphics</code> context in which to paint
   */
  @Override
  public void paint(Graphics g) {
    super.paint(g);
    Graphics2D g2 = (Graphics2D) g;
    if (getWidth() > 100) {
      //横軸の描画
      int width = getWidth();
      int height = getHeight();
      int x_min = 50;
      int x_max = width - 50;
      double y = height * 0.6;
      g2.setColor(Color.BLACK);
      g2.setStroke(new BasicStroke(1.5f));
      g2.draw(new Line2D.Double(x_min, y, x_max, y));
      for (double knot : m_knots) {
        //目盛りの描画
        double x_knot = x_min + (knot - m_knots[0]) * (x_max - x_min) / (m_knots[m_knots.length - 1] - m_knots[0]);
        double length;
        if (height < 100) {
          length = height * 0.2;
        } else {
          length = 20;
        }
        g2.draw(new Line2D.Double(x_knot, y, x_knot, y - length));
        //フォーマットは各自いじって使う
        String s = String.format("%2.2f", knot);
        g.drawString(s, (int) x_knot, (int) y + 20);
      }

      //点列の描画
      g2.setColor(Color.RED);
      g2.setStroke(new BasicStroke(1.0f));
      for (Point p : m_points) {
        double time = p.time();
        //点列の描画
        double x_point = x_min + (time - m_knots[0]) * (x_max - x_min) / (m_knots[m_knots.length - 1] - m_knots[0]);
        double length;
        if (height < 100) {
          length = height * 0.1;
        } else {
          length = 10;
        }
        g2.draw(new Line2D.Double(x_point, y, x_point, y - length));
      }
    }
  }

  /**
   * 点列のコピーを取得します.
   *
   * @return 点列のコピー
   */
  public Point[] getPoints() {
    return m_points.clone();
  }

  /**
   * 節点列のコピーを取得します.
   *
   * @return 節点列のコピー
   */
  public double[] getKnots() {
    return m_knots.clone();
  }

  /**
   * コンストラクタ
   *
   * @param _points 点列
   * @param _knots  節点列
   */
  private PointsGraph(Point[] _points, double[] _knots) {
    m_points = _points;
    m_knots = _knots;
    //サイズは使いやすい大きさに調整してください.小さすぎると何も表示しないです.
    setPreferredSize(new Dimension(WIDTH, HEIGHT));
  }

  /** 点列 */
  private final Point[] m_points;
  /** 節点列 */
  private final double[] m_knots;
  /** 横方向のサイズ */
  private static final int WIDTH = 800;
  /** 縦方向のサイズ */
  private static final int HEIGHT = 100;
}
