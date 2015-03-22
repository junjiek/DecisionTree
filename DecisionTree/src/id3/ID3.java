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
		public TreeNode parent = null;			    //父节点
		public String decompositionAttribute = "";  //当前节点分类属性
		public String pDecompositionValue = "";     //父节点分类属性值
		public ArrayList<TreeNode> children = new ArrayList<TreeNode>();  //子节点列表
		public String classLabel = "";      //若当前节点为叶节点，则该节点表示的类别
	}

	private TreeNode treeRoot = null;
	private enum attributeType {
		CONTINUOUS, DESCRETE
	}
	private ArrayList<String> attributes = new ArrayList<String>(); //存储属性名
	private ArrayList<attributeType> attributeTypes = new ArrayList<attributeType>(); //存储属性类型 (DESCRETE, CONTINUOUS)
	private ArrayList<ArrayList> attributeValues = new ArrayList<ArrayList>();  //存储每个属性的取值
	private ArrayList<String[]> data = new ArrayList<String[]>();  //存储String格式数据

	private int classAttributeIdx = -1;    //分类属性在attributes列表中的索引

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

	//构建分类决策树
	public TreeNode buildDecisionTree(LinkedList<Integer> selattr, ArrayList<Integer> subset, String pDecompositionValue) {
		TreeNode node = new TreeNode();
		node.pDecompositionValue = pDecompositionValue;

		//如果subset中所有数据都属于同一类
		if (classPure(subset)) {
			node.classLabel = data.get(subset.get(0))[classAttributeIdx];
			//System.out.println(node.pDecompositionValue + "\t" + node.classLabel);
			return node;
		} 

		//如果selattr候选分类属性集为空
		if (selattr.size() == 0) {
			node.classLabel = MajorityVoting(subset);//多数表决
			return node;
		}

		//计算各属性的信息增益，并从中选择信息增益最大的属性作为分类属性
		int maxIndex = -1;
		double maxEntropy = Double.MIN_VALUE;
		for (int i = 0; i < selattr.size(); i++) {
			double entropy = calInformationGain(subset, selattr.get(i));
			if (entropy > maxEntropy) {
				maxIndex = selattr.get(i);
				maxEntropy = entropy;
			}
		}
		//划分
		node.decompositionAttribute = attributes.get(maxIndex);
		selattr.remove(new Integer(maxIndex));
		ArrayList<String> attrval = attributeValues.get(maxIndex);
		for (String val : attrval) {
			ArrayList<Integer> subsubset = new ArrayList<Integer>();
			for (int i = 0; i < subset.size(); i++) {
				if (data.get(subset.get(i))[maxIndex].equals(val)) {
					subsubset.add(subset.get(i));
				}
			}
			if (subsubset.size() != 0) {
				TreeNode child = buildDecisionTree(selattr, subsubset, val);
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
		System.out.println(node.pDecompositionValue + "\t" + node.decompositionAttribute);
		if (node.children.size() != 0) {
			for (TreeNode tnode : node.children) {
				printTree(tnode);
			}
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

	// 需要参数：ARFF格式数据文件的文件名
	public static void main(String[] args) {

		ID3 id3 = new ID3();
		//读取ARFF格式数据文件
		id3.readC45("./data/adult.names", "./data/small.data");
//		id3.readARFF("./data/weather.nominal.arff");
		id3.printData();
//		id3.setClassAttribute("play");
//
//		// 构建分类决策树
//		LinkedList<Integer> selattr = new LinkedList<Integer>();//当前节点可用分类属性集
//		for (int i = 0; i < id3.attributes.size(); i++) {
//			if (i != id3.classAttributeIdx) {
//				selattr.add(i);
//			}
//		}
//		ArrayList<Integer> subset = new ArrayList<Integer>(id3.data.size());//当前节点训练数据集
//		for (int i = 0; i < id3.data.size(); i++) {
//			subset.add(i);
//		}
//		id3.treeRoot = id3.buildDecisionTree(selattr, subset, "");
//		id3.printTree(id3.treeRoot);
//		System.out.println();
//		id3.nprintTree(id3.treeRoot);
	}
}
