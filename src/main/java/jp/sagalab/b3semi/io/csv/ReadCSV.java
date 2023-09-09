package jp.sagalab.b3semi.io.csv;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import jp.sagalab.jftk.Point;

/**
 * CSVファイルの読み取りを行うクラスです.
 *
 */
public class ReadCSV {
  /**
   * CSV形式の点列ファイルの読み込みを行う.
   *
   * @param _pointsFile 点列のCSVファイル
   * @return 点列のリスト
   */
  public static List<Point> readPoints(File _pointsFile) {
    if (_pointsFile == null || !(_pointsFile.exists())) {
      throw new IllegalArgumentException("File does not exist");
    }
    if (!(_pointsFile.getName().endsWith(".csv"))) {
      throw new IllegalArgumentException("Points file is not csv file");
    }
    BufferedReader br = null;
    List<Point> pointsList = new ArrayList<>();
    try {
      br = new BufferedReader(new FileReader(_pointsFile));
      String line;
      String[] data;
      int row = 0;
      while ((line = br.readLine()) != null) {
        data = line.split(",");
        if (row == 0) {
          if (data.length != 3) {
            throw new IllegalArgumentException("Column size is not 3");
          }
          if (!data[0].strip().equals("x")) {
            throw new IllegalArgumentException("First element is not x");
          }
          if (!data[1].strip().equals("y")) {
            throw new IllegalArgumentException("Second element is not y");
          }
          if (!data[2].strip().equals("parameter")) {
            throw new IllegalArgumentException("Third element is not parameter");
          }
        } else {
          double x = Double.parseDouble(data[0]);
          double y = Double.parseDouble(data[1]);
          double t = Double.parseDouble(data[2]);
          pointsList.add(Point.createXYT(x, y, t));
        }
        row++;
      }
    } catch (IOException e) {
      throw new RuntimeException(e);
    } finally {
      if (br != null) {
        try {
          br.close();
        } catch (IOException e) {
          throw new RuntimeException(e);
        }
      }
    }
    return pointsList;
  }

  /**
   * CSV形式のknot列の読み込みを行う.
   *
   * @param _knotFile knot列のCSVファイル
   * @return knot列
   */
  public static double[] readKnots(File _knotFile) {
    if (_knotFile == null || !(_knotFile.exists())) {
      throw new IllegalArgumentException("File does not exist");
    }
    if (!(_knotFile.getName().endsWith(".csv"))) {
      throw new IllegalArgumentException("Knot file is not csv file");
    }
    BufferedReader br = null;
    List<double[]> knotList = new ArrayList<>();
    try {
      br = new BufferedReader(new FileReader(_knotFile));
      String line;
      String[] data;
      int row = 0;
      while ((line = br.readLine()) != null) {
        data = line.split(",");
        if (row == 0) {
          if (data.length != 1) {
            throw new IllegalArgumentException("column size is not 1");
          }
          if (!data[0].strip().equals("knot")) {
            throw new IllegalArgumentException("first row must be 'knot'");
          }
        } else {
          double[] knot = new double[data.length];
          for (int i = 0; i < data.length; i++) {
            knot[i] = Double.parseDouble(data[i]);
          }
          knotList.add(knot);
        }
        row++;
      }
    } catch (IOException e) {
      throw new RuntimeException(e);
    } finally {
      if (br != null) {
        try {
          br.close();
        } catch (IOException e) {
          throw new RuntimeException(e);
        }
      }
    }
    if (knotList.isEmpty()) {
      throw new IllegalArgumentException("knotList is empty");
    }
    if (knotList.size() != 1) {
      throw new IllegalArgumentException("knotList size is not 1");
    }
    return knotList.get(0);
  }
}
