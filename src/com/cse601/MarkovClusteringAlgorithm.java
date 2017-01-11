package com.cse601;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;
import java.util.TreeMap;

public class MarkovClusteringAlgorithm {

	static double[][] attweb = new double[180][180];

	static double[][] yeast = new double[359][359];

	static double[][] physicsData = new double[142][142];

	static int[] listAtt = new int[180];

	static int[] listYeast = new int[359];

	static int[] listPhysicsData = new int[142];
	static TreeMap<Integer, Integer> yeastMap = new TreeMap<>();
	static TreeMap<String, Integer> PhysicsMap = new TreeMap<>();

	public static void main(String[] args) throws IOException,
			NumberFormatException {

		Scanner scanner = new Scanner(System.in);
		System.out.println("Enter the value of parameter e (Integer value) > 1: ");
		int e = scanner.nextInt();
		System.out.println("Enter the value of parameter r (Double value): ");
		double r = scanner.nextDouble();
		System.out.println();
		System.out.println("<<<AT&T Network Data>>>");
		generateMatrixForNumbers("input/attweb_net.txt", attweb, " ");
		ArrayList<ArrayList<Integer>> opAtt = executeAttWeb(attweb, r, e);

		listAtt = fillOutput(opAtt, listAtt);

		writeToFile(listAtt, "att_output", "input/attweb_net.net");
		System.out.println();
		System.out.println("<<<Yeast Metabolic Data>>>");
		
		generateMatrixForYeast("input/yeast_undirected_metabolic.txt", yeast,
				"\t");
		ArrayList<ArrayList<Integer>> opYeast = executeAttWeb(yeast, r, e);
		listYeast = fillOutput(opYeast, listYeast);
		writeToFileYeast(listYeast, "yeast_output",
				"input/yeast_undirected_metabolic.net");
		
		System.out.println();
		System.out.println("<<<Physics Collaboration Data>>>");
		generateMatrixForPhysics("input/physics_collaboration_net.txt",
				physicsData, " ");
		ArrayList<ArrayList<Integer>> opPhysics = executeAttWeb(physicsData,
				r, e);
		listPhysicsData = fillOutput(opPhysics, listPhysicsData);
		writeToFilePhy(listPhysicsData, "Phy_output",
				"input/physics_collaboration_net.net");
	}

	private static int[] fillOutput(ArrayList<ArrayList<Integer>> clusters,
			int[] list) {

		for (int i = 0; i < clusters.size(); i++) {
			for (int j = 0; j < clusters.get(i).size(); j++) {
				list[clusters.get(i).get(j)] = i + 1;
			}

		}

		return list;
	}

	private static ArrayList<ArrayList<Integer>> executeAttWeb(
			double[][] matrix, double r, int expansionFactor) {
		matrix = normalizeMatrix(matrix);

		double[][] output;
		double[][] prev = matrix;
		int count = 0;
		for (int i = 0; i < 100; i++) {
			output = expandMatrix(matrix, expansionFactor);
			matrix = inflateMatrix(output, r);
			if (count > 1 && checkIfConverged(prev, matrix)) {
				break;
			}
			count++;
			prev = matrix;
		}
		System.out.println("Loop terminated after: "+count + " iterations");

		return getClusters(matrix);
	}

	private static boolean checkIfConverged(double[][] matrix1,
			double[][] matrix2) {
		int rows = matrix1.length;
		int cols = matrix1[0].length;
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < cols; j++) {
				if (Math.abs(matrix1[i][j] - matrix2[i][j]) > 0) {
					return false;
				}
			}
		}
		return true;
	}

	private static ArrayList<ArrayList<Integer>> getClusters(double[][] matrix) {
		int rows = matrix.length;
		int cols = matrix[0].length;
		ArrayList<ArrayList<Integer>> clusters = new ArrayList<>();
		for (int i = 0; i < rows; i++) {
			ArrayList<Integer> list = new ArrayList<Integer>();
			for (int j = 0; j < cols; j++) {
				if (matrix[i][j] > 0) {
					list.add(j);
				}
			}
			if (list.size() > 0) {
				clusters.add(list);
			}
		}
		System.out.println("Number of clusters: "+clusters.size());
		for (int i = 0; i < clusters.size(); i++) {
			System.out.println(clusters.get(i));
		}

		return clusters;

	}

	private static double[][] expandMatrix(double[][] matrix, int factor) {
		int rows = matrix.length;
		int columns = matrix[0].length;
		double[][] expandedMatrix = new double[rows][columns];
		double[][] prev;
		prev = matrix;
		for (int l = 0; l < factor - 1; l++) {
			for (int i = 0; i < rows; i++) {
				for (int j = 0; j < columns; j++) {
					double sum = 0;
					for (int k = 0; k < columns; k++) {
						sum += prev[i][k] * matrix[k][j];
					}
					expandedMatrix[i][j] = sum;
				}
			}
			prev = expandedMatrix;
		}
		return expandedMatrix;
	}

	private static double[][] inflateMatrix(double[][] matrix, double factor) {
		int rows = matrix.length;
		int columns = matrix[0].length;
		double[][] inflatedMatrix = new double[rows][columns];
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {
				matrix[i][j] = Math.pow(matrix[i][j], factor);
			}
		}

		inflatedMatrix = normalizeMatrix(matrix);

		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {
				if (inflatedMatrix[i][j] < 0.001) {
					inflatedMatrix[i][j] = 0;
				}
			}
		}

		return inflatedMatrix;
	}

	private static double[][] normalizeMatrix(double[][] matrix) {
		int rows = matrix.length;
		int columns = matrix[0].length;
		double[][] inflatedMatrix = new double[rows][columns];
		double[] columnSums = new double[matrix[0].length];
		for (int i = 0; i < columns; i++) {
			double sum = 0;
			for (int j = 0; j < rows; j++) {
				sum += matrix[j][i];
			}
			columnSums[i] = sum;
		}

		for (int i = 0; i < columns; i++) {
			for (int j = 0; j < rows; j++) {
				inflatedMatrix[j][i] = matrix[j][i] / columnSums[i];
			}
		}

		return inflatedMatrix;
	}

	public static void printMatrix(double[][] matrix) {
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix[0].length; j++) {
				System.out.print(matrix[i][j] + " ");
			}
			// System.out.println();
		}
	}

	public static double[][] generateMatrixForNumbers(String fileName,
			double[][] matrix, String delimiter) throws IOException {
		String line = null;

		FileReader fileReader = new FileReader(fileName);
		@SuppressWarnings("resource")
		BufferedReader bufferedReader = new BufferedReader(fileReader);
		while ((line = bufferedReader.readLine()) != null) {
			String sA[] = line.split(delimiter);
			int num1 = Integer.parseInt(sA[0]);
			int num2 = Integer.parseInt(sA[1]);
			matrix[num1][num2] = matrix[num2][num1] = 1;

		}

		for (int i = 0; i < matrix.length; i++) {
			matrix[i][i] = 1;
		}

		return matrix;
	}

	public static TreeMap<Integer, Integer> fillMapForYeast(String fileName,
			double[][] matrix, String delimiter, TreeMap<Integer, Integer> map)
			throws NumberFormatException, IOException {
		int counter = 0;
		String line = null;

		FileReader fileReader = new FileReader(fileName);
		@SuppressWarnings("resource")
		BufferedReader bufferedReader = new BufferedReader(fileReader);
		while ((line = bufferedReader.readLine()) != null) {
			String sA[] = line.split(delimiter);
			int num1 = Integer.parseInt(sA[0]);
			int num2 = Integer.parseInt(sA[1]);

			if (!map.containsKey(num1)) {
				map.put(num1, counter++);
			}
			if (!map.containsKey(num2)) {
				map.put(num2, counter++);
			}
		}

		return map;
	}

	public static double[][] generateMatrixForYeast(String fileName,
			double[][] matrix, String delimiter) throws NumberFormatException,
			IOException {

		yeastMap = fillMapForYeast(fileName, matrix, delimiter, yeastMap);
		// System.out.println(yeastMap);
		String line = null;

		FileReader fileReader = new FileReader(fileName);
		@SuppressWarnings("resource")
		BufferedReader bufferedReader = new BufferedReader(fileReader);
		while ((line = bufferedReader.readLine()) != null) {
			String sA[] = line.split(delimiter);
			int num1 = yeastMap.get(Integer.parseInt(sA[0]));
			int num2 = yeastMap.get(Integer.parseInt(sA[1]));

			matrix[num1][num2] = matrix[num2][num1] = 1;

		}

		for (int i = 0; i < matrix.length; i++) {
			matrix[i][i] = 1;
		}

		return matrix;
	}

	// 0 1 138 24
	public static TreeMap<String, Integer> fillMapForPhysics(String fileName,
			double[][] matrix, String delimiter, TreeMap<String, Integer> map)
			throws NumberFormatException, IOException {
		int counter = 0;
		String line = null;

		FileReader fileReader = new FileReader(fileName);
		@SuppressWarnings("resource")
		BufferedReader bufferedReader = new BufferedReader(fileReader);
		while ((line = bufferedReader.readLine()) != null) {
			String sA[] = line.split(delimiter);

			if (!map.containsKey(sA[0])) {
				map.put(sA[0], counter++);
			}
			if (!map.containsKey(sA[1])) {
				map.put(sA[1], counter++);
			}
		}

		return map;
	}

	public static double[][] generateMatrixForPhysics(String fileName,
			double[][] matrix, String delimiter) throws NumberFormatException,
			IOException {

		PhysicsMap = fillMapForPhysics(fileName, matrix, delimiter, PhysicsMap);
		// System.out.println(PhysicsMap);
		String line = null;

		FileReader fileReader = new FileReader(fileName);
		@SuppressWarnings("resource")
		BufferedReader bufferedReader = new BufferedReader(fileReader);
		while ((line = bufferedReader.readLine()) != null) {
			String sA[] = line.split(delimiter);
			int num1 = PhysicsMap.get(sA[0]);
			int num2 = PhysicsMap.get(sA[1]);

			matrix[num1][num2] = matrix[num2][num1] = 1;

		}

		for (int i = 0; i < matrix.length; i++) {
			matrix[i][i] = 1;
		}

		return matrix;
	}

	private static void writeToFile(int[] list, String fileName, String netFile)
			throws NumberFormatException, IOException {
		try {
			PrintWriter writer = new PrintWriter("output/" + fileName + ".clu",
					"UTF-8");
			ArrayList<Integer> vList = readfile(netFile);
			writer.println("*Vertices 180");
			for (int i : vList) {
				writer.println(list[i]);
//				System.out.print(list[i] + " ");
			}

			writer.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private static void writeToFileYeast(int[] list, String fileName,
			String netFile) throws NumberFormatException, IOException {
		try {
			PrintWriter writer = new PrintWriter("output/" + fileName + ".clu",
					"UTF-8");
			ArrayList<Integer> vList = readfile(netFile);
			writer.println("*Vertices 359");
			for (int i : vList) {
				writer.println(list[yeastMap.get(i)]);
//				System.out.print(list[yeastMap.get(i)] + " ");
			}

			writer.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private static void writeToFilePhy(int[] list, String fileName,
			String netFile) throws NumberFormatException, IOException {
		try {
			PrintWriter writer = new PrintWriter("output/" + fileName + ".clu",
					"UTF-8");
			ArrayList<String> vList = readfilePhy(netFile);
			writer.println("*Vertices 142");
			for (String i : vList) {
				writer.println(list[PhysicsMap.get(i)]);
//				System.out.print(list[PhysicsMap.get(i)] + " ");
			}

			writer.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private static ArrayList<Integer> readfile(String fileName)
			throws NumberFormatException, IOException {
		FileReader fileReader = new FileReader(fileName);
		ArrayList<Integer> vertexArray = new ArrayList<>();
		String line = "";
		@SuppressWarnings("resource")
		BufferedReader bufferedReader = new BufferedReader(fileReader);
		while ((line = bufferedReader.readLine()) != null) {
			String[] lineArray = line.split("\"");
			if (lineArray.length > 1) {
				int vertex = Integer.parseInt(lineArray[1]);
				// System.out.println(vertex);
				vertexArray.add(vertex);

			}
		}
		return vertexArray;

	}

	private static ArrayList<String> readfilePhy(String fileName)
			throws NumberFormatException, IOException {
		FileReader fileReader = new FileReader(fileName);
		ArrayList<String> vertexArray = new ArrayList<>();
		String line = "";
		@SuppressWarnings("resource")
		BufferedReader bufferedReader = new BufferedReader(fileReader);
		while ((line = bufferedReader.readLine()) != null) {
			String[] lineArray = line.split("\"");
			if (lineArray.length > 1) {
				String vertex = lineArray[1];
				// System.out.println(vertex);
				vertexArray.add(vertex);

			}
		}
		return vertexArray;

	}

}
