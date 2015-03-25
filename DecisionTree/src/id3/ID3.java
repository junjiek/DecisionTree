package id3;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Pattern;
import java.util.regex.Matcher;

import javax.jws.Oneway;

public class ID3 {
	
	private enum AttributeType {
		CONTINUOUS, DISCRETE
	}
	private static final int UNKNOWN = Integer.MAX_VALUE;
	private enum CompareType {
		EQ, LT, GE
	}
	
	class TreeNode {
		
		public TreeNode parent = null;			 //父节点
		public int decompositionAttribute = -1;  //当前节点分类属性
		public double pDecompositionValue = -1;  //父节点分类属性值
		public CompareType type = CompareType.EQ;
		public ArrayList<TreeNode> children = new ArrayList<TreeNode>();  //子节点列表
		public String classLabel = "";      //若当前节点为叶节点，则该节点表示的类别
	}

	private TreeNode treeRoot = null;

	// private ArrayList<String> classTypes = new ArrayList<String>();  //存储分类属性种类
	private ArrayList<String> attributes = new ArrayList<String>();  //存储属性名
	private ArrayList<AttributeType> attributeTypes = new ArrayList<AttributeType>();  //存储属性类型 (DISCRETE, CONTINUOUS)
	private ArrayList<ArrayList<String>> attributeValues = new ArrayList<ArrayList<String>>();  //存储每个属性的取值
	private ArrayList<String[]> data = new ArrayList<String[]>();  //存储String格式数据
	private int classAttributeIdx = -1;    //分类属性在data列表中的索引
	private double[] newSplitPoint;
	//读取ARFF格式数据文件
	public void readARFF(String filename) {
		try {
			BufferedReader reader = new BufferedReader(new FileReader(filename));
			Pattern pattern = Pattern.compile("@attribute(.*)[{](.*?)[}]");
			String line = "";
			while ((line = reader.readLine()) != null) {
				Matcher matcher = pattern.matcher(line);
				if (matcher.find()) {
					attributes.add(matcher.group(1).trim());
					String[] values = matcher.group(2).split(",");;
					ArrayList<String> list = new ArrayList<String>(values.length);
					for (String value : values) {
						list.add(value.trim());
					}
					attributeValues.add(list);
				} else if (line.startsWith("@data")) {
					while ((line = reader.readLine()) != null) {
						if (line == "") {
							continue;
						}
						String[] row = line.split(",");
						data.add(row);
					}
				} else {
					continue;
				}
			}
			reader.close();
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}

	public void readC45(String namesFile, String dataFile) {
		// read *.names file
		try {
			BufferedReader reader = new BufferedReader(new FileReader(namesFile));
			String line = "";
			boolean firstEntry = true;
			ArrayList<String> classValues = new ArrayList<String>();
			while ((line = reader.readLine()) != null) {
				line = line.replaceAll("\\s*", "");
				int comment = line.indexOf('|');
				if (comment >= 0) {
					line = line.substring(0, comment);
				}
				if (line.length() == 0) continue;
				if (line.endsWith(".")) line = line.substring(0, line.length()-1);
				if (firstEntry) {
					String[] values = line.split(",");
					for (String value : values) {
						classValues.add(value.trim());
					}
					firstEntry = false;
					continue;
				}
				int nameIndex = line.indexOf(':');
				attributes.add(line.substring(0, nameIndex));
				line = line.substring(nameIndex+1, line.length());
				// System.out.println(line);
				if(line.toUpperCase().compareTo("CONTINUOUS") == 0) {
					attributeTypes.add(AttributeType.CONTINUOUS);
					ArrayList<String> list = new ArrayList<String>();
					attributeValues.add(list);
				} else {
					attributeTypes.add(AttributeType.DISCRETE);
					String[] values = line.split(",");
					ArrayList<String> list = new ArrayList<String>();
					for (String value : values) {
						list.add(value.trim());
					}
					attributeValues.add(list);
				}
			}
			attributes.add("classLabel");
			attributeTypes.add(AttributeType.DISCRETE);
			attributeValues.add(classValues);
			setClassAttribute("classLabel");
			newSplitPoint = new double[attributes.size()];
			reader.close();
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}

		// read *.data file
		try {
			BufferedReader reader = new BufferedReader(new FileReader(dataFile));
			String line = "";
			
			while ((line = reader.readLine()) != null) {		
				if (line.endsWith(".")) line = line.substring(0, line.length()-1);
				String[] row = line.trim().split(",");
				data.add(row);
			}
			reader.close();
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}
	}

	//打印数据，以确定数据读取是否正确
	public void printData() {
		System.out.println("@attributes");
		for (String attribute : attributes) {
			System.out.println(attribute);
		}
		System.out.println();

		System.out.println("@attributeValues");
		for (ArrayList list : attributeValues) {
			for (Object value : list) {
				System.out.print(value.toString() + " ");
			}
			System.out.println();
		}
		System.out.println();

		System.out.println("@data");
		for (String[] array : data) {
			for (String d : array) {
				System.out.print(d + " ");
			}
			System.out.println();
		}
		System.out.println();
		System.out.println("@data size: " + data.size());
		System.out.println("@attribute size: " + attributes.size());
		System.out.println("@attributeTypes size: " + attributeTypes.size());
		System.out.println("@attributeValues size: " + attributeValues.size());
		System.out.println("@classAttributeIdx: " + classAttributeIdx);
	}

	//设置分类属性在attributes中的索引
	public void setClassAttributeIdx(int idx) {
		if (idx < 0 || idx >= attributes.size()) {
			System.err.println("class attribute set error!");
			System.exit(1);
		}
		classAttributeIdx = idx;
	}

	public void setClassAttribute(String classAttribute) {
		int n = attributes.indexOf(classAttribute);
		setClassAttributeIdx(n);
	}

	//判断subset中数据是否属于同一类
	public boolean classPure(ArrayList<Integer> subset) {
		String classLabel = data.get(subset.get(0))[classAttributeIdx];
		for (int i = 1; i < subset.size(); i++) {
			String nextClassLabel = data.get(subset.get(i))[classAttributeIdx];
			if (!nextClassLabel.equals(classLabel)) {
				return false;
			}
		}
		return true;
	}

	//统计不同类别计数
	public int[] classCount(ArrayList<Integer> subset) {
		ArrayList<String> attrval = attributeValues.get(classAttributeIdx);
		int[] count = new int[attrval.size()];
		for (int i = 0; i < subset.size(); i++) {
			String classLabel = data.get(subset.get(i))[classAttributeIdx];
			count[attrval.indexOf(classLabel)]++;
		}
		return count;
	}

	//多数表决判定结点类别
	public String MajorityVoting(ArrayList<Integer> subset) {
		int[] count = classCount(subset);

		int maxIdx = 0;
		for (int i = 1; i < count.length; i++) {
			if (count[i] > count[maxIdx]) {
				maxIdx = i;
			}
		}

		return (String)attributeValues.get(classAttributeIdx).get(maxIdx);
	}

	//计算熵
	public double calEntropy(int[] count) {
		double entropy = 0.0;
		int total = 0;
		for (int num : count) {
			total += num;
		}

		for (int num : count) {
			if (num == 0 || total == 0) {
				return 0;
			}
			entropy += -(num * 1.0 / total) * Math.log(num * 1.0 / total) / Math.log(2); 
		}

		return entropy;
	}

	public int compareTo(Object attributeValue, int attributeIndex, String[] oneRecord) {
		if (oneRecord[attributeIndex].compareTo("?") == 0) return UNKNOWN;
		if (attributeTypes.get(attributeIndex) == AttributeType.DISCRETE) {
			return oneRecord[attributeIndex].compareTo((String)attributeValue);
		}
		Double recordValue = Double.parseDouble(oneRecord[attributeIndex]);
		return recordValue.compareTo((Double)attributeValue);
	}

	private class Pair implements Comparable{
		Double value;
		Integer index;
		public Pair(Double v, Integer i) {
			value = v;
			index = i;
		}
		public int compareTo(Object obj) {
			Pair tmp = (Pair)obj;
			if (tmp.value.compareTo(value) != 0){
				return value.compareTo(tmp.value);	
			}
				return index.compareTo(tmp.index);
		}
		public String toString() {
			return "(" + value.toString() + ", " + index.toString() + ")";
		}
	}

	public double calSubsetInfo(TreeSet<Pair> subset, double splitPoint) {
		ArrayList<String> classattrval = attributeValues.get(classAttributeIdx);
		// < splitPoint
		SortedSet<Pair> lowerSet = subset.headSet(new Pair(splitPoint, -1));
		// >= splitPoint
		SortedSet<Pair> upperSet = subset.tailSet(new Pair(splitPoint, -1));
		int info[][] = new int[2][classattrval.size()];
		Iterator<Pair> iter = lowerSet.iterator();
		while(iter.hasNext()) {
			int instance = iter.next().index;
			String classAttr = data.get(instance)[classAttributeIdx];
			info[0][classattrval.indexOf(classAttr)]++;	
		}
		iter = upperSet.iterator();
		while(iter.hasNext()) {
			int instance = iter.next().index;
			String classAttr = data.get(instance)[classAttributeIdx];
			info[1][classattrval.indexOf(classAttr)]++;	
		}
		return (calEntropy(info[0])*lowerSet.size()+
				calEntropy(info[1])*upperSet.size())/subset.size();
	}

	public double calContinuousInfoGain(ArrayList<Integer> subset, int index) {
		//整体的熵
		double infoD = calEntropy(classCount(subset));
		
		TreeSet<Pair> values = new TreeSet<Pair>();
		TreeSet<Double> diffValues = new TreeSet<Double>(); 
		for (int i : subset) {
			String dataStr = data.get(i)[index];
			if (dataStr.compareTo("?") == 0) continue;
			Double tmp = Double.parseDouble(dataStr);
			values.add(new Pair(tmp, i));
			diffValues.add(tmp);
		}
		if(diffValues.size() < 2) return 0;
		// System.out.println("------\n"  + diffValues);
		double min_info = Double.MAX_VALUE;
		double splitPoint = Double.MAX_VALUE;
		Iterator iter = diffValues.iterator();
		iter.next();  //跳过第一个分割点
		while(iter.hasNext()) {
			Double point = (Double)iter.next();
			double entropy = calSubsetInfo(values, point);
			// System.out.println(entropy + " " + splitPoint);
			if (entropy < min_info){
				splitPoint = point;
				min_info = entropy;
			}
		}
		newSplitPoint[index] = splitPoint;
		System.out.println(attributes.get(index) + ": " + splitPoint);
		//由属性index划分后的熵
		return infoD - min_info;
	}

	//计算信息增益
	public double calDiscreteInfoGain(ArrayList<Integer> subset, int index) {
		//整体的熵
		double infoD = calEntropy(classCount(subset));
		//由属性index划分后的熵
		ArrayList<String> classattrval = attributeValues.get(classAttributeIdx);
		ArrayList<String> attrval = attributeValues.get(index);
		int[][] info = new int[attrval.size()][classattrval.size()];
		int[] count = new int[attrval.size()];
		for (int i = 0; i < subset.size(); i++) {
			int n = subset.get(i);
			int attrvalIndex;
			if (data.get(n)[index].compareTo("?") == 0) continue;
			attrvalIndex = attrval.indexOf(data.get(n)[index]);
			int classattrvalIndex = classattrval.indexOf(data.get(n)[classAttributeIdx]);
			info[attrvalIndex][classattrvalIndex]++;
			count[attrvalIndex]++;
		}
		int sum = subset.size();
		double infoDA = 0.0;
		for (int i = 0; i < attrval.size(); i++) {
			infoDA += calEntropy(info[i]) * count[i] / sum;
		}
		return infoD - infoDA;
	}

	// 构建分类决策树
	public void buildDecisionTree() {
		HashSet<Integer> selattr = new HashSet<Integer>();  //初始可用分类属性集为全集
		for (int i = 0; i < attributes.size(); i++) {
			if (i != classAttributeIdx) {
				selattr.add(i);
			}
		}
		ArrayList<Integer> subset = new ArrayList<Integer>(data.size());  //初始数据集
		for (int i = 0; i < data.size(); i++) {
			subset.add(i);
		}
		treeRoot = buildSubDecisionTree(selattr, subset);
		treeRoot.pDecompositionValue = -1.0;
	}

	//构建子集的分类决策树
	public TreeNode buildSubDecisionTree(HashSet<Integer> selattr, ArrayList<Integer> subset) {
		TreeNode node = new TreeNode();
		//如果subset中所有数据都属于同一类
		if (classPure(subset)) {
			node.classLabel = data.get(subset.get(0))[classAttributeIdx];
			return node;
		} 

		//如果selattr候选分类属性集为空
		if (selattr.size() == 0) {
			node.classLabel = MajorityVoting(subset);//多数表决
			return node;
		}

		//计算各属性的信息增益，并从中选择信息增益最大的属性作为分类属性
		System.out.println("Calculating max infoGain...");
		int maxIndex = -1;
		double maxInfoGain = -1.0;
		for (int i : selattr) {
			double infoGain;
			System.out.println("----- " + attributes.get(i) + " ------");
			if (attributeTypes.get(i) == AttributeType.CONTINUOUS)
				infoGain = calContinuousInfoGain(subset, i);
			else
				infoGain = calDiscreteInfoGain(subset, i);
			System.out.println("infoGain: " + infoGain);
			if (infoGain > maxInfoGain) {
				maxIndex = i;
				maxInfoGain = infoGain;
			}
		}

		//划分
		System.out.println("Choose \"" + attributes.get(maxIndex)
						  + "\" to decompose");
		node.decompositionAttribute = maxIndex;
		selattr.remove(new Integer(maxIndex));
		if (attributeTypes.get(maxIndex) == AttributeType.CONTINUOUS) {
			double splitPoint = newSplitPoint[maxIndex];
			ArrayList<Integer> lowerSubset = new ArrayList<Integer>();
			ArrayList<Integer> upperSubset = new ArrayList<Integer>();
			for (int j : subset) {
				if (compareTo(splitPoint, maxIndex, data.get(j)) >= 0)
					upperSubset.add(j);
				else
					lowerSubset.add(j);
			}

			TreeNode lowerChild;
			if (lowerSubset.size() != 0)
				lowerChild = buildSubDecisionTree(new HashSet<Integer>(selattr), lowerSubset);
			else {
				lowerChild = new TreeNode();
				lowerChild.classLabel = MajorityVoting(subset);
			}
			lowerChild.pDecompositionValue = splitPoint;
			lowerChild.type = CompareType.LT;
			lowerChild.parent = node;
			node.children.add(lowerChild);

			TreeNode upperChild;
			if (upperSubset.size() != 0)
				upperChild = buildSubDecisionTree(new HashSet<Integer>(selattr), upperSubset);
			else {
				upperChild = new TreeNode();
				upperChild.classLabel = MajorityVoting(subset);
			}
			upperChild.pDecompositionValue = splitPoint;
			upperChild.type = CompareType.GE;
			upperChild.parent = node;
			node.children.add(upperChild);
		} else {
			ArrayList<String> attrval = attributeValues.get(maxIndex);
			for (int i = 0; i < attrval.size(); i++) {
				ArrayList<Integer> subsubset = new ArrayList<Integer>();
				for (int j : subset) {
					if (compareTo(attrval.get(i), maxIndex, data.get(j)) == 0) {
						subsubset.add(j);
					}
				}
				System.out.println(attrval.get(i) + ": " + subsubset.size());
				TreeNode child;
				if (subsubset.size() != 0)
					child = buildSubDecisionTree(new HashSet<Integer>(selattr), subsubset);
				else {
					child = new TreeNode();
					child.classLabel = MajorityVoting(subset);
				}
				child.pDecompositionValue = i;
				child.type = CompareType.EQ;
				child.parent = node;
				node.children.add(child);
			}
		}
		

		return node;
	}


	public void printTree(TreeNode node, String tab, BufferedWriter out) throws IOException{
		if (node.children.size() == 0) {
			out.write(tab + "\t" + attributes.get(classAttributeIdx)
					+ " = \"" + node.classLabel+ "\";");
			out.newLine();
			return;
		}
		int childsize = node.children.size();
		for (int i = 0; i < childsize; i++) {
			TreeNode child = node.children.get(i);
			String classifier = "";
			switch (child.type) {
				case EQ: classifier = " == "; break;
				case GE: classifier = " >= "; break;
				case LT: classifier = " < "; break;
			}
			String pDecompositionValue;
			if (attributeTypes.get(node.decompositionAttribute) == AttributeType.CONTINUOUS)
				pDecompositionValue = child.pDecompositionValue + "";
			else {
				ArrayList<String> value = attributeValues.get(node.decompositionAttribute);
				pDecompositionValue = value.get((int)child.pDecompositionValue);
			}
			out.write(tab + "if( "
					+ attributes.get(node.decompositionAttribute) + classifier
					+ "\"" + pDecompositionValue + "\") {");
			out.newLine();
			printTree(child, tab + "\t", out);
			if (i != childsize - 1) {
				out.write(tab + "} else ");
				out.newLine();
			}
			else {
				out.write(tab + "}");
				out.newLine();
			}
		}

	}

	public static void main(String[] args) {
		long startTime = System.currentTimeMillis();
		ID3 id3 = new ID3();
		// 读取ARFF格式数据文件
		// id3.readARFF("./data/weather.nominal.arff");
		// id3.setClassAttribute("play");
		// 读取C4.5格式数据文件
		id3.readC45("./data/adult.names", "./data/adult.data");
		// id3.printData();
		// 构建分类决策树
		id3.buildDecisionTree();

		try {
			BufferedWriter myout = new BufferedWriter(new FileWriter(new File("./tree.txt")));       
			id3.printTree(id3.treeRoot, "", myout); 
		} catch(Exception e) {
			System.out.println(e);
		}
		
		long endTime = System.currentTimeMillis();
		System.out.println("========= Tree built, running time: "
						    + (endTime - startTime)/1000.0 + "s =========");

		// id3.readC45("")
	}
}
