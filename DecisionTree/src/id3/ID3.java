package id3;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Random;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Pattern;
import java.util.regex.Matcher;

import javax.jws.Oneway;
import javax.naming.ldap.Rdn;

public class ID3 {
	
	private static final int UNKNOWN = Integer.MAX_VALUE;
	private static final double zalpha2 = 1.150;  // p = 0.25 two-sided
	private enum AttributeType {
		CONTINUOUS, DISCRETE
	}
	private enum CompareType {
		EQ, LT, GE
	}

	
	class TreeNode {
		public TreeNode parent = null;       //父节点
		public int decomposeAttribute = -1;  //当前节点分类属性
		public String pDecomposeValue = "";  //父节点分类属性值
		public CompareType type = CompareType.EQ;
		public ArrayList<TreeNode> children = new ArrayList<TreeNode>();  //子节点列表
		public boolean leaf = false;
		public String classLabel = "";  //若当前节点为叶节点，则该节点表示的类别
		public int[] data = new int[2];  // 保存validation set的信息，方便起见，直接默认是两个分类
	}
	private int treeSize = 0;
	private TreeNode treeRoot = null;
	private HashSet<Integer> trainSet;
	private HashSet<Integer> validationSet;
	// private ArrayList<String> classTypes = new ArrayList<String>();  //存储分类属性种类
	private ArrayList<String> attributes = new ArrayList<String>();  //存储属性名
	private ArrayList<AttributeType> attributeTypes = new ArrayList<AttributeType>();  //存储属性类型 (DISCRETE, CONTINUOUS)
	private ArrayList<ArrayList<String>> attributeValues = new ArrayList<ArrayList<String>>();  //存储每个属性的取值
	private ArrayList<String[]> traindata = new ArrayList<String[]>();  //存储String格式数据
	private ArrayList<String[]> testdata = new ArrayList<String[]>();
	private int classAttributeIdx = -1;    //分类属性在data列表中的索引
	private double[] newSplitPoint;

	public void readC45(String namesFile, String dataFile, String testFile) {
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
				traindata.add(row);
			}
			reader.close();
		} catch (IOException ioe) {
			ioe.printStackTrace();
		}

		// read *.test file
		try {
			BufferedReader reader = new BufferedReader(new FileReader(testFile));
			String line = "";
			
			while ((line = reader.readLine()) != null) {		
				if (line.endsWith(".")) line = line.substring(0, line.length()-1);
				String[] row = line.trim().split(",");
				testdata.add(row);
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
		for (String[] array : traindata) {
			for (String d : array) {
				System.out.print(d + " ");
			}
			System.out.println();
		}
		System.out.println();
		System.out.println("@data size: " + traindata.size());
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
	public boolean classPure(HashSet<Integer> subset) {
		String classLabel = traindata.get(subset.iterator().next())[classAttributeIdx];
		for (int i : subset) {
			String nextClassLabel = traindata.get(i)[classAttributeIdx];
			if (!nextClassLabel.equals(classLabel)) {
				return false;
			}
		}
		return true;
	}

	//统计不同类别计数
	public int[] classCount(HashSet<Integer> subset) {
		ArrayList<String> attrval = attributeValues.get(classAttributeIdx);
		int[] count = new int[attrval.size()];
		for (int i : subset) {
			String classLabel = traindata.get(i)[classAttributeIdx];
			count[attrval.indexOf(classLabel)]++;
		}
		return count;
	}

	//多数表决判定结点类别
	public String MajorityVoting(HashSet<Integer> subset) {
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

	public int compareTo(String attributeValue, int attributeIndex, String[] oneRecord) {
		if (oneRecord[attributeIndex].compareTo("?") == 0) return UNKNOWN;
		if (attributeTypes.get(attributeIndex) == AttributeType.DISCRETE) {
			return oneRecord[attributeIndex].compareTo((String)attributeValue);
		}
		Double recordValue = Double.parseDouble(oneRecord[attributeIndex]);
		return recordValue.compareTo(Double.parseDouble(attributeValue));
	}

	private class Pair implements Comparable {
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
			String classAttr = traindata.get(instance)[classAttributeIdx];
			info[0][classattrval.indexOf(classAttr)]++;	
		}
		iter = upperSet.iterator();
		while(iter.hasNext()) {
			int instance = iter.next().index;
			String classAttr = traindata.get(instance)[classAttributeIdx];
			info[1][classattrval.indexOf(classAttr)]++;	
		}
		return (calEntropy(info[0])*lowerSet.size()+
				calEntropy(info[1])*upperSet.size())/subset.size();
	}

	public double calContinuousInfoGain(HashSet<Integer> subset, int index) {
		//整体的熵
		double infoD = calEntropy(classCount(subset));
		
		TreeSet<Pair> values = new TreeSet<Pair>();
		TreeSet<Double> diffValues = new TreeSet<Double>(); 
		for (int i : subset) {
			String dataStr = traindata.get(i)[index];
			if (dataStr.compareTo("?") == 0) continue;
			Double tmp = Double.parseDouble(dataStr);
			values.add(new Pair(tmp, i));
			diffValues.add(tmp);
		}
		if(diffValues.size() < 2) return 0;
		// 取中点为划分点
		ArrayList<Double> splitPoints = new ArrayList<Double>();
		Iterator<Double> iter = diffValues.iterator();
		Double p1 = iter.next();
		while (iter.hasNext()) {
			Double p2 = iter.next();
			splitPoints.add((p1 + p2) / 2);
			p1 = p2;
		}
		double min_info = Double.MAX_VALUE;
		double splitPoint = Double.MAX_VALUE;
		for (double point : splitPoints) {
			double entropy = calSubsetInfo(values, point);
			if (entropy < min_info){
				splitPoint = point;
				min_info = entropy;
			}
		}
		// 直接取右值为划分点
		// Iterator<Double> iter = diffValues.iterator();
		// iter.next();  //跳过第一个分割点
		// while(iter.hasNext()) {
		// 	Double point = (Double)iter.next();
		// 	double entropy = calSubsetInfo(values, point);
		// 	if (entropy < min_info){
		// 		splitPoint = point;
		// 		min_info = entropy;
		// 	}
		// }

		newSplitPoint[index] = splitPoint;
		// System.out.println("splitPoint: " + splitPoint);
		//由属性index划分后的熵
		return infoD - min_info;
	}

	double splitInfoDiscrete(HashSet<Integer> subset, int attr) {
		ArrayList<String> attrvals = attributeValues.get(attr);
		int counts[] = new int[attrvals.size()];
		int size = subset.size();
		int missDataSize = 0;
		for (int i : subset) {
			String val = traindata.get(i)[attr];
			if (val.compareTo("?") == 0) {
				missDataSize ++;
				continue;
			}
			int index = attrvals.indexOf(val);
			counts[index] ++;
		}
		// 将缺失数据取为最可能的属性
		int mostCommon = getMostCommonIdx(counts);
		counts[mostCommon] += missDataSize;

		double info = 0.0;
		for (int count : counts) {            
			double rate = count / (double) size;
			if (rate != 0 && rate != 1) {
				info += (-1 * rate) * (Math.log(rate) / Math.log(2.0));
			}
			//System.out.println(info);
		}
		if (info == 0.0) {
			return .000001;
		}
		return info;
	}

	double splitInfoContinuous(HashSet<Integer> subset, int attr, double split) {
		ArrayList<String> attrvals = attributeValues.get(attr);
		int counts[] = new int[2];
		int missDataSize = 0;
		for (int i : subset) {
			String valstr = traindata.get(i)[attr];
			if (valstr.compareTo("?") == 0) {
				missDataSize ++;
				continue;
			}
			Double val = Double.parseDouble(valstr);
			if (val < split)
				counts[0] ++;
			else
				counts[1] ++;
		}
		// 将缺失数据取为最可能的属性
		int mostCommon = getMostCommonIdx(counts);
		counts[mostCommon] += missDataSize;

		double info = 0.0;
		for (int count : counts) {            
			double rate = count / (double) subset.size();
			if (rate != 0 && rate != 1) {
				info += (-1 * rate) * (Math.log(rate) / Math.log(2.0));
			}
		}
		if (info == 0.0) {
			return .000001;
		}
		return info;
	}

	private int getMostCommonIdx(int[] count) {
		int max = count[0];
		int index = 0;
		for (int i = 1; i < count.length; i++) {
			if (max < count[i]) {
				index = i;
				max = count[i];
			}
		}
		return index;
	}

	//计算信息增益
	public double calDiscreteInfoGain(HashSet<Integer> subset, int index) {
		// 整体的熵
		double infoD = calEntropy(classCount(subset));
		// 统计划分后的取值情况
		ArrayList<String> classattrval = attributeValues.get(classAttributeIdx);
		ArrayList<String> attrval = attributeValues.get(index);
		int[][] info = new int[attrval.size()][classattrval.size()];
		int[] count = new int[attrval.size()];
		int[] missData = new int[attrval.size()];
		for (int i : subset) {
			int classattrvalIndex = classattrval.indexOf(traindata.get(i)[classAttributeIdx]);
			if (traindata.get(i)[index].compareTo("?") == 0) {
				missData[classattrvalIndex] ++;
				continue;
			}
			int attrvalIndex;
			attrvalIndex = attrval.indexOf(traindata.get(i)[index]);
			info[attrvalIndex][classattrvalIndex]++;
			count[attrvalIndex]++;
		}
		// 将缺失数据取值为最可能的属性值
		int mostCommon = getMostCommonIdx(count);
		for (int  i : missData)
			count[mostCommon] += i;
		for (int i = 0; i < classattrval.size(); i++) {
			info[mostCommon][i] += missData[i];
		}
		// 计算由属性index划分后的熵
		int sum = subset.size();
		double infoDA = 0.0;
		for (int i = 0; i < attrval.size(); i++) {
			infoDA += calEntropy(info[i]) * count[i] / sum;
		}
		return infoD - infoDA;
	}

	private void generateTrainTestSet(double trainRate, double validationRate) {
		// 选取train set
		int trainSetSize = (int)(trainRate * traindata.size());
		trainSet = new HashSet<Integer>(trainSetSize);
		Random rdm = new Random(System.currentTimeMillis());
		for (int i = 0; i < trainSetSize; i++) {
			Integer train = rdm.nextInt() % traindata.size();
			if (train < 0) train += traindata.size();
			while (trainSet.contains(train)) {
				train = rdm.nextInt() % traindata.size();
				if (train < 0) train += traindata.size();
			}
			trainSet.add(train);
		}
		// 在train set中选取validation set
		if (validationRate <= 0) return;
		int validationSetSize = (int)(validationRate * trainSetSize);
		validationSet = new HashSet<Integer>(validationSetSize);
		ArrayList<Integer> trainsetList = new ArrayList<Integer>(trainSet);
		for (int i = 0; i < validationSetSize; i++) {
			Integer rdmIdx = rdm.nextInt() % trainSetSize;
			if (rdmIdx < 0) rdmIdx += trainSetSize;
			while (validationSet.contains(trainsetList.get(rdmIdx))) {
				rdmIdx = rdm.nextInt() % trainSetSize;
				if (rdmIdx < 0) rdmIdx += trainSetSize;
			}
			Integer validation = trainsetList.get(rdmIdx);
			validationSet.add(validation);
			trainSet.remove(validation);
		}
		return;
	}
	// 构建分类决策树, 随机选择一部分数据作为训练集
	public void train() {
		HashSet<Integer> selattr = new HashSet<Integer>();  //初始可用分类属性集为全集
		for (int i = 0; i < attributes.size(); i++) {
			if (i != classAttributeIdx) {
				selattr.add(i);
			}
		}
		treeSize = 0;
		treeRoot = buildDecisionTree(selattr, trainSet);
		treeRoot.pDecomposeValue = "";
	}

	private String classifyOneRecord(String[] record) throws IOException{
		// String debug = "---------------------------\n";
		// for (String str : record)
		// 	debug += (str + ", ");
		// debug += "\n";
		TreeNode node = treeRoot;
		while (!node.leaf) {
			int attr = node.decomposeAttribute;
			// debug += attributes.get(attr);
			if (record[attr].compareTo("?") == 0) {
				// debug += ( "--( ? )-->");
				break;
			}
			TreeNode child;
			for (int i = 0; i < node.children.size(); i++) {
				child = node.children.get(i);
				int cmp = compareTo(child.pDecomposeValue, attr, record);
				boolean match1 = child.type == CompareType.EQ && cmp == 0;
				boolean match2 = child.type == CompareType.GE && cmp >= 0;
				boolean match3 = child.type == CompareType.LT && cmp < 0;
				// String sym = "";
				// if (match2) {
				// 	sym = record[attr] + " >= ";
				// }
				// if (match3) {
				// 	sym = record[attr] + " < ";
				// }
				if (match1 || match2 || match3) {
					// debug += ( "--(" + sym + child.pDecomposeValue + ")-->");
					node = child;
					break;
				}
			}
		}
		// debug += " " + node.classLabel;
		// out.write(debug);
		// out.newLine();
		return node.classLabel;
	}

	public double test(BufferedWriter out) throws IOException{
		int correct = 0;
		for (String[] test : testdata) {
			String predict = classifyOneRecord(test);
			String real = test[classAttributeIdx];
			out.write("p: " + predict + ", r: " + real);
			out.newLine();
			if (predict.compareTo(real) == 0) 
				correct ++;
		}
		return correct * 1.0 / testdata.size();
	}

	private void traceOneRecord(String[] record) {
		ArrayList<String> classattrval = attributeValues.get(classAttributeIdx);
		int classval = classattrval.indexOf(record[classAttributeIdx]);
		TreeNode node = treeRoot;
		while (!node.leaf) {
			node.data[classval] ++;
			int attr = node.decomposeAttribute;
			if (record[attr].compareTo("?") == 0) break;
			TreeNode child;
			for (int i = 0; i < node.children.size(); i++) {
				child = node.children.get(i);
				int cmp = compareTo(child.pDecomposeValue, attr, record);
				boolean match1 = child.type == CompareType.EQ && cmp == 0;
				boolean match2 = child.type == CompareType.GE && cmp >= 0;
				boolean match3 = child.type == CompareType.LT && cmp < 0;
				if (match1 || match2 || match3) {
					node = child;
					break;
				}
			}
		}
		node.data[classval] ++;
		return;
	}

	private void classifyValidationSet() {
		for(int i : validationSet) {
			traceOneRecord(traindata.get(i));
		}
	}

	public void pruning() {
		classifyValidationSet();
		// PEPrune(treeRoot);
		REPrune(treeRoot);
	}


	private int REPrune(TreeNode node) {
		ArrayList<String> classattrval = attributeValues.get(classAttributeIdx);
		int label = classattrval.indexOf(node.classLabel);
		int rightNum = node.data[label];
		if (node.leaf) {
			return rightNum;
		}
		int tot = node.data[0] + node.data[1];
		if (tot == 0) return rightNum;
		int childCorrect = 0;
		for (TreeNode child : node.children) {
			childCorrect += REPrune(child);
		}
		if (rightNum > childCorrect) {
			//prune
			treeSize -= node.children.size();
			node.leaf = true;
			// node.classLabel = classattrval.get(label);
			return rightNum;
		}
		// don't prune
		return childCorrect;
	}

	private double PEPrune(TreeNode node) {
		double prunedError = pessimisticError(node);
		if (node.leaf)
			return prunedError;
		double childError = 0.0;
		for(TreeNode child : node.children) {
			childError += PEPrune(child);
		}
		if (prunedError < childError) {
			// prune
			// System.out.println("pruned");
			treeSize -= node.children.size();
			node.leaf = true;
			return prunedError;
		}
		// don't prune
		return childError;
	}

	private double pessimisticError(TreeNode node) {
		int n = node.data[0] + node.data[1]; 
		double wrong = Math.min(node.data[0], node.data[1]);
		double p = (wrong + 1.0)/(n + 2.0);
		return n * (p + zalpha2 * Math.sqrt(p*(1.0-p) / (n+2.0)));
	}


	//构建子集的分类决策树
	public TreeNode buildDecisionTree(HashSet<Integer> selattr, HashSet<Integer> subset) {
		TreeNode node = new TreeNode();
		//如果subset中所有数据都属于同一类
		if (classPure(subset)) {
			node.classLabel = traindata.get(subset.iterator().next())[classAttributeIdx];
			node.leaf = true;
			return node;
		} 

		node.classLabel = MajorityVoting(subset);//多数表决
		//如果selattr候选分类属性集为空
		if (selattr.size() == 0) {
			node.leaf = true;
			return node;
		}

		//计算各属性的信息增益，并从中选择信息增益率最大的属性作为分类属性
		// System.out.println("Calculating max gainRate...");
		int maxIndex = -1;
		double maxGainRate = -1.0;
		for (int i : selattr) {
			double gain, splitInfo, gainRate;
			// System.out.println("----- " + attributes.get(i) + " ------");
			if (attributeTypes.get(i) == AttributeType.CONTINUOUS) {
				gain = calContinuousInfoGain(subset, i);
				splitInfo = splitInfoContinuous(subset, i, newSplitPoint[i]);
			}
			else {
				gain = calDiscreteInfoGain(subset, i);
				splitInfo = splitInfoDiscrete(subset, i);
			}
			gainRate = gain / splitInfo;
			// System.out.println("gainRate: " + gainRate);
			if (gainRate > maxGainRate) {
				maxIndex = i;
				maxGainRate = gainRate;
			}
		}

		//划分
		// System.out.println("Choose \"" + attributes.get(maxIndex)
						  // + "\" to decompose");
		node.decomposeAttribute = maxIndex;
		selattr.remove(new Integer(maxIndex));
		HashSet<Integer> unknownSet = new HashSet<Integer>();
		if (attributeTypes.get(maxIndex) == AttributeType.CONTINUOUS) {
			treeSize += 2;
			double splitPoint = newSplitPoint[maxIndex];
			HashSet<Integer> lowerSubset = new HashSet<Integer>();
			HashSet<Integer> upperSubset = new HashSet<Integer>();
			for (int j : subset) {
				int cmp = compareTo(splitPoint + "", maxIndex, traindata.get(j));
				if (cmp == UNKNOWN) {
					unknownSet.add(j);
					continue;
				}
				if (cmp >= 0)
					upperSubset.add(j);
				else
					lowerSubset.add(j);
			}
			// 将缺失属性的样本加入到多数的一边
			if (upperSubset.size() > lowerSubset.size())
				upperSubset.addAll(unknownSet);
			else
				lowerSubset.addAll(unknownSet);

			TreeNode lowerChild;
			if (lowerSubset.size() != 0)
				lowerChild = buildDecisionTree(new HashSet<Integer>(selattr), lowerSubset);
			else {
				lowerChild = new TreeNode();
				lowerChild.leaf = true;
				lowerChild.classLabel = MajorityVoting(subset);
			}
			lowerChild.pDecomposeValue = splitPoint + "";
			lowerChild.type = CompareType.LT;
			lowerChild.parent = node;
			node.children.add(lowerChild);

			TreeNode upperChild;
			if (upperSubset.size() != 0)
				upperChild = buildDecisionTree(new HashSet<Integer>(selattr), upperSubset);
			else {
				upperChild = new TreeNode();
				upperChild.leaf = true;
				upperChild.classLabel = MajorityVoting(subset);
			}
			upperChild.pDecomposeValue = splitPoint + "";
			upperChild.type = CompareType.GE;
			upperChild.parent = node;
			node.children.add(upperChild);
		} else {
			for (String attrval : attributeValues.get(maxIndex)) {
				treeSize ++;
				HashSet<Integer> subsubset = new HashSet<Integer>();
				for (int j : subset) {
					if (compareTo(attrval, maxIndex, traindata.get(j)) == 0) {
						subsubset.add(j);
					}
				}
				// System.out.println(attrval + ": " + subsubset.size());
				TreeNode child;
				if (subsubset.size() != 0)
					child = buildDecisionTree(new HashSet<Integer>(selattr), subsubset);
				else {
					child = new TreeNode();
					child.leaf = true;
					child.classLabel = MajorityVoting(subset);
				}
				child.pDecomposeValue = attrval;
				child.type = CompareType.EQ;
				child.parent = node;
				node.children.add(child);
			}
		}
		
		return node;
	}


	public void printTreeToFile(TreeNode node, String tab, BufferedWriter out) throws IOException{
		if (node.leaf) {
			out.write(tab + attributes.get(classAttributeIdx)
					+ " = \"" + node.classLabel+ "\";");
					// + "(" + node.data[0] + ", " + node.data[1] + ")");
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
			out.write(tab + "if( "
					+ attributes.get(node.decomposeAttribute) + classifier
					+ "\"" + child.pDecomposeValue + "\") {");
			out.newLine();
			printTreeToFile(child, tab + "\t", out);
			if (i != childsize - 1) {
				out.write(tab + "} else ");
			}
			else {
				out.write(tab + "}");
				out.newLine();
			}
		}

	}

	public void blackboxTest(BufferedWriter out) throws IOException {
		for (String[] test : testdata) {
			String predict = classifyOneRecord(test);
			out.write(predict + ".");
			out.newLine();
		}
	}

	public void task3() {
		readC45("./data/adult.names", "./data/adult.data", "./data/test.features");
		int validationSetSize = (int)(0.4 * traindata.size());
		validationSet = new HashSet<Integer>(validationSetSize);
		trainSet = new HashSet<Integer>(traindata.size() - validationSetSize);
		for (int i = 0; i < validationSetSize; i++)
			validationSet.add(i);
		for (int i = validationSetSize; i < traindata.size(); i++)
			trainSet.add(i);
		train();
		pruning();
		try {
			BufferedWriter testout = new BufferedWriter(new FileWriter(new File("./2012011335.test.result")));       
			blackboxTest(testout);
			testout.close();
		} catch(Exception e) {
			System.out.println(e);
		}
		return;
	}

	public double task1(double trainRate) {
		long startTime = System.currentTimeMillis();
		// 读取C4.5格式数据文件
		readC45("./data/adult.names", "./data/adult.data", "./data/adult.test");
		// printData();
		// 构建分类决策树
		generateTrainTestSet(trainRate, 0);
		train();

		try {
			BufferedWriter treeout = new BufferedWriter(new FileWriter(new File("./tree")));       
			printTreeToFile(treeRoot, "", treeout); 
			treeout.close();
		} catch(Exception e) {
			System.out.println(e);
		}
		
		long endTime = System.currentTimeMillis();
		System.out.println("========= Tree built, size:" + treeSize +
							" , running time: " + (endTime - startTime)/1000.0
							 + "s =========");
		
//		System.out.println("testing...");
		double accuracy = 0.0;
		try {
			BufferedWriter testout = new BufferedWriter(new FileWriter(new File("./test")));       
			accuracy = test(testout);
			System.out.println("accuracy : " + accuracy);
			testout.close();
		} catch(Exception e) {
			System.out.println(e);
		}
		return accuracy;
	}

	public double task2(double trainRate, double validationRate) {
		long startTime = System.currentTimeMillis();
		// 读取C4.5格式数据文件
		readC45("./data/adult.names", "./data/adult.data", "./data/adult.test");
		// printData();
		// 构建分类决策树
		generateTrainTestSet(trainRate, validationRate);
		train();

		try {
			BufferedWriter treeout = new BufferedWriter(new FileWriter(new File("./tree")));       
			printTreeToFile(treeRoot, "", treeout); 
			treeout.close();
		} catch(Exception e) {
			System.out.println(e);
		}
		
		long endTime = System.currentTimeMillis();
		System.out.println("========= Tree built, size:" + treeSize +
							" , running time: " + (endTime - startTime)/1000.0
							 + "s =========");
		
		System.out.println("pruning...");
		pruning();
		System.out.println("========= Pruned size:" + treeSize +" =========");
		System.out.println("testing...");
		double accuracy = 0.0;
		try {
			BufferedWriter testout = new BufferedWriter(new FileWriter(new File("./test")));       
			accuracy = test(testout);
			System.out.println("accuracy : " + accuracy);
			testout.close();
		} catch(Exception e) {
			System.out.println(e);
		}
		return accuracy;
	}

	public static void main(String[] args) {
		try {
			PrintStream myout = new PrintStream(new FileOutputStream(new File("./log3")));       
			System.setOut(myout);        
			System.setErr(myout);
		} catch (FileNotFoundException e) {
			System.out.println(e);
		}
		for (int i = 10; i < 20; i++) {
			double trainRate = i * 0.05;
			double tmp = 0.0;
			for (int j = 0; j < 5; j++) {
				ID3 id3 = new ID3();
				tmp += id3.task1(trainRate);				
			}
			System.out.println( "****** " + trainRate + ": " + tmp/5);
		}
		// id3.task2(0.05, 0.4);
		// id3.task3();
		
		
	}
}
