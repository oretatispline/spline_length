package jp.sagalab.b3semi.io.csv;

//import jp.sagalab.b3semi.Point;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.List;

import jp.sagalab.jftk.Point;

public class WriteToCSV {
  public static File writePoints(List<Point> _points) {
    if (_points == null) {
      throw new IllegalArgumentException("_pointsList is null");
    }
    String fileName = "points_" + System.currentTimeMillis();
    PrintWriter pw = null;
    File pointsDir = new File(".\\Points");
    if (pointsDir.mkdir()) {
      System.out.println("create Points directory");
    }
    File pointsFile = new File(".\\Points\\" + fileName + ".csv");
    try {
      pw = new PrintWriter(pointsFile);
      pw.println("x,y,parameter");
      for (Point point: _points) {
        pw.println(point.x() + "," + point.y() + "," + point.time());
      }
      pw.flush();
    } catch (FileNotFoundException e) {
      throw new RuntimeException(e);
    } finally {
      if (pw != null) {
        pw.close();
      }
    }
    return pointsFile;
  }

  public static File writeKnots(double[] _knot) {
    String fileName = "knot_" + System.currentTimeMillis();
    PrintWriter pw = null;
    File knotDir = new File(".\\Knot");
    if (knotDir.mkdir()) {
      System.out.println("create knot directory");
    }
    File knotFile = new File(".\\Knot\\" + fileName + ".csv");
    try {
      pw = new PrintWriter(knotFile);
      pw.println("knot");
      for (int i = 0; i < _knot.length; i++) {
        if (i != _knot.length - 1) {
          pw.print(_knot[i]+",");
        } else {
          pw.println(_knot[i]);
        }
      }
      pw.flush();
    } catch (FileNotFoundException e) {
      throw new RuntimeException(e);
    } finally {
      if (pw != null) {
        pw.close();
      }
    }
    return knotFile;
  }
}
