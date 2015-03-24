package id3;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.regex.Pattern;
import java.util.regex.Matcher;

import javax.jws.Oneway;

public class ID3 {

	class TreeNode {
		public TreeNode parent = null;			 //父节点
		public int decompositionAttribute = -1;  //当前节点分类属性
		public int pDecompositionValue = -1;     //父节点分类属性值
		public ArrayList<TreeNode> children = new ArrayList<TreeNode>();  //子节点列表
		public String classLabel = "";      //若当前节点为叶节点，则该节点表示的类别
	}

	private TreeNode treeRoot = null;
	private enum attributeType {
		CONTINUOUS, DESCRETE
	}

	// private ArrayList<String> classTypes = new ArrayList<String>();  //存储分类属性种类
	private ArrayList<String> attributes = new ArrayList<String>();  //存储属性名
	private ArrayList<attributeType> attributeTypes = new ArrayList<attributeType>();  //存储属性类型 (DESCRETE, CONTINUOUS)
	private ArrayList<ArrayList> attributeValues = new ArrayList<ArrayList>();  //存储每个属性的取值
	private ArrayList<String[]> data = new ArrayList<String[]>();  //存储String格式数据
	private int classAttributeIdx = -1;    //分类属性在data列表中的索引

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
					attributeTypes.add(attributeType.CONTINUOUS);
					ArrayList<Double> list = new ArrayList<Double>();
					attributeValues.add(list);
				} else {
					attributeTypes.add(attributeType.DESCRETE);
					String[] values = line.split(",");
					ArrayList<String> list = new ArrayList<String>();
					for (String value : values) {
						list.add(value.trim());
					}
					attributeValues.add(list);
				}
			}
			attributes.add("classLabel");
			attributeTypes.add(attributeType.DESCRETE);
			attributeValues.add(classValues);
			setClassAttribute("classLabel");
			reader.close();
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}

		// read *.data file
		try {
			BufferedReader reader = new BufferedReader(new FileReader(dataFile));
			String line = "";
			ArrayList<Integer> contiAttributesIndex = new ArrayList<Integer>();
			for (int i = 0; i < attributeTypes.size(); i++) {
				if(attributeTypes.get(i) == attributeType.CONTINUOUS)
					contiAttributesIndex.add(i);
			}
			while ((line = reader.readLine()) != null) {		
				if (line.endsWith(".")) line = line.substring(0, line.length()-1);
				String[] row = line.trim().split(",");
				for (Integer index : contiAttributesIndex) {
					attributeValues.get(index).add(Double.parseDouble(row[index]));
				}
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

	//计算信息增益
	public double calInformationGain(ArrayList<Integer> subset, int index) {
		//整体的熵
		double infoD = calEntropy(classCount(subset));

		//由属性index划分后的熵
		ArrayList<String> classattrval = attributeValues.get(classAttributeIdx);
		ArrayList<String> attrval = attributeValues.get(index);
		int[][] info = new int[attrval.size()][classattrval.size()];
		int[] count = new int[attrval.size()];
		for (int i = 0; i < subset.size(); i++) {
			int n = subset.get(i);
			int attrvalIndex = attrval.indexOf(data.get(n)[index]);
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
		LinkedList<Integer> selattr = new LinkedList<Integer>();  //初始可用分类属性集为全集
		for (int i = 0; i < attributes.size(); i++) {
			if (i != classAttributeIdx) {
				selattr.add(i);
			}
		}
		ArrayList<Integer> subset = new ArrayList<Integer>(data.size());  //初始数据集
		for (int i = 0; i < data.size(); i++) {
			subset.add(i);
		}
		treeRoot = buildSubDecisionTree(selattr, subset, -1);
	}

	//构建子集的分类决策树
	public TreeNode buildSubDecisionTree(LinkedList<Integer> selattr, ArrayList<Integer> subset, int pDecompositionValue) {
		TreeNode node = new TreeNode();
		node.pDecompositionValue = pDecompositionValue;

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
		int maxIndex = -1;
		double maxInfoGain = Double.MIN_VALUE;
		for (int i = 0; i < selattr.size(); i++) {
			double infoGain = calInformationGain(subset, selattr.get(i));
			if (infoGain > maxInfoGain) {
				maxIndex = selattr.get(i);
				maxInfoGain = infoGain;
			}
		}

		//划分
		node.decompositionAttribute = maxIndex;
		selattr.remove(new Integer(maxIndex));
		ArrayList<String> attrval = attributeValues.get(maxIndex);
		for (int i = 0; i < attrval.size(); i++) {
			ArrayList<Integer> subsubset = new ArrayList<Integer>();
			for (int j = 0; j < subset.size(); j++) {
				if (data.get(subset.get(j))[maxIndex].equals(attrval.get(i))) {
					subsubset.add(subset.get(j));
				}
			}
			if (subsubset.size() != 0) {
				TreeNode child = buildSubDecisionTree(selattr, subsubset, i);
				child.parent = node;
				node.children.add(child);
			} else {
				TreeNode child = new TreeNode();
				child.parent = node;
				child.classLabel = MajorityVoting(subset);
				node.children.add(child);
			}
		}

		return node;
	}

	//决策树的先根遍历
	public void printTree(TreeNode node) {
		String pDepositionValue = "", decompositionAttribute = "";
		if (node.pDecompositionValue != -1)
			pDepositionValue = attributes.get(node.pDecompositionValue);
		if (node.decompositionAttribute != -1)
			decompositionAttribute = attributes.get(node.decompositionAttribute);

		System.out.println(pDepositionValue + "\t" + decompositionAttribute);
		if (node.children.size() != 0) {
			for (TreeNode tnode : node.children) {
				printTree(tnode);
			}
		}
	}

	public void printTree(TreeNode node, String tab) {
		if (node.pDecompositionValue != -1) {
			TreeNode parent = node.parent;
			System.out.println(attributeValues.get(parent.decompositionAttribute).get(node.pDecompositionValue)
					+ "\") {");
		}
		if (node.children.size() == 0) {
			System.out.println(tab + "\t" + attributes.get(classAttributeIdx)
					+ " = \"" + node.classLabel+ "\";");
			return;
		}
		int childsize = node.children.size();
		for (int i = 0; i < childsize; i++) {
			System.out.print(tab + "if( "
					+ attributes.get(node.decompositionAttribute) + " == \"");
			printTree(node.children.get(i), tab + "\t");
			if (i != childsize - 1)
				System.out.print(tab + "} else ");
			else
				System.out.println(tab + "}");
		}

	}

	//决策树的后根遍历
	public void nprintTree(TreeNode node) {
		if (node.children.size() != 0) {
			for (TreeNode tnode : node.children) {
				printTree(tnode);
			}
		}
		System.out.println(node.pDecompositionValue + "\t" + node.decompositionAttribute);
	}

	public static void main(String[] args) {

		ID3 id3 = new ID3();
		// 读取ARFF格式数据文件
		// id3.readARFF("./data/weather.nominal.arff");
		// id3.setClassAttribute("play");
		// 读取C4.5格式数据文件
		id3.readC45("./data/weather.names", "./data/weather.data");
		// id3.printData();

		// 构建分类决策树
		id3.buildDecisionTree();
		id3.printTree(id3.treeRoot, "");
		 // System.out.println();
		 // id3.nprintTree(id3.treeRoot);
	}
}
