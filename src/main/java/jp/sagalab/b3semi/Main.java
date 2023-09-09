package jp.sagalab.b3semi;

import jp.sagalab.b3semi.graph.PointsGraph;
import jp.sagalab.b3semi.io.csv.ReadCSV;
import jp.sagalab.b3semi.io.csv.WriteToCSV;

import jp.sagalab.jftk.Point;
import jp.sagalab.jftk.curve.interporation.ConstrainedSplineCurveInterpolator;
import jp.sagalab.jftk.curve.Range;
import jp.sagalab.jftk.curve.SplineCurve;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionAdapter;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Line2D;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 * @author yako, takashima, inagaki
 */
public class Main extends JFrame {

  /**
   * @param _args the command line arguments
   */
  public static void main(String[] _args) {
    new Main();
  }

  /**
   * コンストラクタ
   */
  public Main() {
    setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    m_canvas.setSize(800, 600);
    m_canvas.setBackground(Color.WHITE);
    setTitle("NonUniformConstrainedSplineCurveInterpolator");

    // 入力点列を保存する"save points"ボタン
    JButton savePoints = new JButton("save points");
    // ボタンを押したときの処理
    savePoints.addActionListener(e -> WriteToCSV.writePoints(m_points));

    // 点列を読み込むための"load points"ボタン
    JButton loadPoints = new JButton("load points");
    // ボタンを押したときの処理
    loadPoints.addActionListener(e -> {
      // 読み込みたい点列のCSVファイルのパスをfilePathに入れる.
      String filePath = ".\\Points\\samplePoints.csv";
      File pointsFile = new File(filePath);
      // readPointsメソッドを用いてPointのリストに変換
      List<Point> points = ReadCSV.readPoints(pointsFile);
      // m_Pointsを上書き
      setPoints(points);
      //Canvasを初期化
      clear();
      // 点列を描画
      for (Point point : points) {
        drawPoint(point.x(), point.y(), 1.5, Color.BLACK);
      }
      // スプライン曲線を描画
      drawSplineCurve();
    });

    add(m_canvas, BorderLayout.CENTER);
    JPanel buttonPanel = new JPanel();
    buttonPanel.add(savePoints);
    buttonPanel.add(loadPoints);
    add(buttonPanel, BorderLayout.SOUTH);
    pack();
    setVisible(true);

    // ドラッグする場合マウスのボタン入力は,
    // mousePressed → mouseDragged → mouseReleased の順で処理が行われる.
    // またドラッグ中はマウスカーソルが動くたびにmouseDraggedが呼び出される.
    m_canvas.addMouseListener(
            new MouseAdapter() {
              @Override
              public void mousePressed(MouseEvent e) {
                clear();
                m_points.clear();
              }

              @Override
              public void mouseReleased(MouseEvent e) {
                drawSplineCurve();
              }
            }
    );

    m_canvas.addMouseMotionListener(
            new MouseMotionAdapter() {
              @Override
              public void mouseDragged(MouseEvent e) {
                Point point = Point.createXYT(e.getX(), e.getY(), System.currentTimeMillis() * 0.001);
                m_points.add(point);
                drawPoint(point.x(), point.y(), 1.5, Color.BLACK);
              }
            }
    );
  }

  public static void createPointsGraph(Point[] _points, double[] _knots) {
    PointsGraph pointsGraph = PointsGraph.create(_points, _knots);
    POINTS_GRAPH_FRAME.getContentPane().removeAll();
    POINTS_GRAPH_FRAME.getContentPane().add(pointsGraph);
    POINTS_GRAPH_FRAME.pack();
    POINTS_GRAPH_FRAME.setVisible(true);
  }

  /**
   * スプライン曲線を求め、描画を行う.
   */
  public void drawSplineCurve() {
    double L = 0;
    for (int i = 0; i<m_points.size()-1; i++){
      L += distance(m_points.get(i).x(),m_points.get(i+1).x(),m_points.get(i).y(),m_points.get(i+1).y());
    }

    //距離を合計したリスト
    double dis = 0.0;
    List<Double> disList = new ArrayList<>();
    for (int i = 0; i < m_points.size();i++){
      if (i == 0){
        disList.add(0.0);
      }
      else {
        dis += distance(m_points.get(i-1).x(),m_points.get(i).x(),m_points.get(i-1).y(),m_points.get(i).y());
        disList.add(dis);
      }
      System.out.println(disList.get(i));
    }

    // ---------- ↓knotを指定しない場合↓ (節点間隔に合わせて節点列を自動で生成) ----------
    // 分かりやすいように時刻パラメータを0から始まるようにシフトしておく.
    Range LRange = Range.create(0.0, L);
//    // 点列の時系列を正規化する.
    List<Point> normalizedPoints = normalizePoints(LRange);
    List<Point> shiftedPoints = shiftPointsTimeZero();
    //System.out.println(shiftedPoints.get(3).time());
    for (int i = 0; i< shiftedPoints.size(); i++) {
      System.out.println(shiftedPoints.get(i));
    }

    List<Double> kotyolist = new ArrayList<>();
    for (int i = 0; i<shiftedPoints.size()-1;i++){
      for(int n = 0; n <=100; n++){
        if((Math.floor(shiftedPoints.get(i+1).time()*10)/10)-(Math.floor(shiftedPoints.get(i).time()*10)/10) < 0.5){
          if(shiftedPoints.get(i).time() <= 0.1 * n && 0.1 * n < shiftedPoints.get(i+1).time()) {
            double a = (shiftedPoints.get(i + 1).time() - 0.1 * n) / (shiftedPoints.get(i + 1).time() - shiftedPoints.get(i).time());
            double b = (0.1 * n - shiftedPoints.get(i).time()) / (shiftedPoints.get(i + 1).time() - shiftedPoints.get(i).time());
            double v = a * disList.get(i) + b * disList.get(i + 1);
            kotyolist.add(v);
            System.out.println(v);
          }
        }
      }
    }

    //kotyolistをdouble配列に代入
    double[] knot2 = new double[kotyolist.size()+5];
    knot2[0] = 0;
    knot2[1] = 0;
    for(int i=0;i<kotyolist.size();i++){
      knot2[i+2]= kotyolist.get(i);
    }
    knot2[kotyolist.size()+2] = L;
    knot2[kotyolist.size()+3] = L;
    knot2[kotyolist.size()+4] = L;

    for (int i = 0; i < m_points.size() - 1; i++){
      normalizedPoints.set(i,Point.createXYT(normalizedPoints.get(i).x(),normalizedPoints.get(i).y(),disList.get(i)));
    }

    // リストを配列に変換する.
    Point[] points = normalizedPoints.toArray(new Point[0]);

    // 次数
    int degree = 3;

    // 節点間隔
    double knotInterval = L/10;
    //double[] knot = new double[]{-0.2, -0.1, 0.0,L/2 ,(L/10)*6,(L/10)*7,(L/10)*8, (L/10)*9, L,L+10, L+20};

    // スプライン補間を行う
    // SplineCurveInterpolator.interpolateの引数は(点列(Point[]型), 次数(int型), 節点間隔(double型))にする.
    SplineCurve splineCurve = NonUniformConstrainedSplineCurveInterpolator.interpolate(points, degree, knot2);

    //稲垣さんのスプライン補間
    //jp.sagalab.jftk.curve.SplineCurve splineCurve2 = ConstrainedSplineCurveInterpolator.interpolate(points,degree,0.1);
    // ---------- ↑knotを指定しない場合↑ (節点間隔に合わせて節点列を自動で生成) ----------


////    // ++++++++++ ↓knotを指定する場合↓ ++++++++++
////    // 時刻パラメータを正規化しておくと節点を自分で定義しやすい.
//    Range timeRange = Range.create(0.0, L);
////    // 点列の時系列を正規化する.
//    List<Point> normalizedPoints = normalizePoints(timeRange);
////    // リストを配列に変換する.
//    Point[] points = normalizedPoints.toArray(new Point[0]);
////
////    // 次数
//    int degree = 3;
////
////    // 節点を定義する.
//    double[] knot = new double[]{-0.2, -0.1, 0.0,L/2 ,(L/10)*6,(L/10)*7,(L/10)*8, (L/10)*9, L,L+10, L+20};
////    // スプライン補間を行う
////    //SplineCurveInterpolator.interpolateの引数は(点列(Point[]型), 次数(int型), 節点列(double[]型))
//    SplineCurve splineCurve = SplineCurveInterpolator.interpolate(points, degree, knot);
////    // ++++++++++ ↑knotを指定する場合↑ ++++++++++


    // スプライン曲線の評価点を求める↓
    double start = splineCurve.range().start();
    double end = splineCurve.range().end();
    List<Point> evaluateList = new ArrayList<>();



//      for (int i = 0; i < m_points.size() - 1; i++){
//        evaluateList.add(splineCurve.evaluate(disList.get(i)));
//      }

     for (double t = start; t < end; t += 0.01) {
       evaluateList.add(splineCurve.evaluate(t));
    }


    // SplineCurveの描画
    for (int i = 1; i < evaluateList.size(); i++) {
      drawLine(evaluateList.get(i-1), evaluateList.get(i), Color.RED);
    }

    Point[] controlList = splineCurve.controlPoints();
    for (int i = 0; i< controlList.length ;i++){
      drawPoint(controlList[i].x(),controlList[i].y(),3,Color.blue);
    }
//    for (int i = 1; i <= controlList.length - 1; i++) {
//      drawLine(controlList[i-1],controlList[i] , Color.blue);
//    }
    //System.out.println(L);
  }

  public double distance(double _x1, double _x2, double _y1, double _y2){
    double X = _x2 - _x1;
    double Y = _y2 - _y1;
    double L = Math.sqrt(Math.pow(X, 2) + Math.pow(Y, 2));
    return L;
  }

  /**
   * 点を描画する.
   *
   * @param _x      x座標
   * @param _y      y座標
   * @param _radius 点の半径
   * @param _color  点の色
   */
  public void drawPoint(double _x, double _y, double _radius, Color _color) {
    Graphics2D g = (Graphics2D) m_canvas.getGraphics();
    g.setColor(_color);

    Ellipse2D.Double oval = new Ellipse2D.Double(_x - _radius, _y - _radius, _radius * 2, _radius * 2);
    g.draw(oval);
  }

  /**
   * 線を描画する.
   *
   * @param _p1    始点
   * @param _p2    終点
   * @param _color 線の色
   */
  public void drawLine(Point _p1, Point _p2, Color _color) {
    Graphics2D g = (Graphics2D)m_canvas.getGraphics();
    g.setColor(_color);

    Line2D.Double line = new Line2D.Double(_p1.x(), _p1.y(), _p2.x(), _p2.y());
    g.draw(line);
  }

  /**
   * 描画を初期化する.
   */
  public void clear() {
    Graphics g = m_canvas.getGraphics();

    g.clearRect(0, 0, m_canvas.getWidth(), m_canvas.getHeight());
  }

  /**
   * 点列の時刻パラメータが0始まりになるように全体をシフトします.
   */
  public List<Point> shiftPointsTimeZero() {
    return normalizePoints(Range.create(0, m_points.get(m_points.size() - 1).time() - m_points.get(0).time()));
  }

  /**
   * 点列の時刻パラメータの正規化をします.
   * m_points全体の時刻パラメータが_range区間に収まるように正規化します.
   *
   * @param _range 正規化後の時刻パラメータの範囲
   */
  public List<Point> normalizePoints(Range _range) {
    double startTime = m_points.get(0).time();
    double timeLength = m_points.get(m_points.size() - 1).time() - startTime;
    double rangeLength = _range.length();
    List<Point> points = new ArrayList<>();
    for (Point point : m_points) {
      points.add(Point.createXYT(point.x(), point.y()
              , _range.start() + (point.time() - startTime) * (rangeLength / timeLength)));
    }

    return points;
  }

  /**
   * m_pointsを引数の点列で上書きする.
   *
   * @param _points 点列
   */
  public void setPoints(List<Point> _points) {
    //入力点列のチェック
    if (_points == null) {
      throw new IllegalArgumentException("_points is null.");
    }
    // 入力点列にnullが混入していないかチェック
    if (_points.contains(null)) {
      throw new IllegalArgumentException(" points include null ");
    }
    // 点列の要素数チェック
    if (_points.size() < 2) {
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

    m_points = new ArrayList<>(_points);
  }

  /** PointsGraphを保持するためのJFrame */
  private static final JFrame POINTS_GRAPH_FRAME = new JFrame();

  /** キャンバスを表す変数 */
  private final Canvas m_canvas = new Canvas();

  /** ドラッグで打たれた点列を保持するリスト */
  private List<Point> m_points = new ArrayList<>();
}