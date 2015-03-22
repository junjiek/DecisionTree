package DecisionTree;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.concurrent.ExecutionException;
import java.math.BigDecimal;
import static java.lang.Math.*;

/**
 * 选择最佳分裂属性
 * 
 * @author Pjq
 * @qq 378290226
 * @mail 378290226@qq.com
 * @date 2012.04.17
 */
public class Gain {
	private ArrayList<ArrayList<String>> D = null; // 训练元组
	private ArrayList<String> attrList = null; // 候选属性集

	public Gain(ArrayList<ArrayList<String>> datas, ArrayList<String> attrList) {
		this.D = datas;
		this.attrList = attrList;
	}

	/**
	 * 获取最佳侯选属性列上的值域（假定所有属性列上的值都是有限的名词或分类类型的）
	 * 
	 * @param attrIndex
	 *            指定的属性列的索引
	 * @return 值域集合
	 */
	public ArrayList<String> getValues(ArrayList<ArrayList<String>> datas,
			int attrIndex) {
		ArrayList<String> values = new ArrayList<String>();
		String r = "";
		for (int i = 0; i < datas.size(); i++) {
			r = datas.get(i).get(attrIndex);
			if (!values.contains(r)) {
				values.add(r);
			}
		}
		return values;
	}

	/**
	 * 获取指定数据集中指定属性列索引的域值及其计数
	 * 
	 * @param d
	 *            指定的数据集
	 * @param attrIndex
	 *            指定的属性列索引
	 * @return 类别及其计数的map
	 */
	public Map<String, Integer> valueCounts(ArrayList<ArrayList<String>> datas,
			int attrIndex) {
		Map<String, Integer> valueCount = new HashMap<String, Integer>();
		String c = "";
		ArrayList<String> tuple = null;
		for (int i = 0; i < datas.size(); i++) {
			tuple = datas.get(i);
			c = tuple.get(attrIndex);
			if (valueCount.containsKey(c)) {
				valueCount.put(c, valueCount.get(c) + 1);
			} else {
				valueCount.put(c, 1);
			}
		}
		return valueCount;
	}

	/**
	 * 获取指定属性列上指定值域的所有元组
	 * 
	 * @param attrIndex
	 *            指定属性列索引
	 * @param value
	 *            指定属性列的值域
	 * @return 指定属性列上指定值域的所有元组
	 */
	public ArrayList<ArrayList<String>> datasOfValue(int attrIndex, String value) {
		ArrayList<ArrayList<String>> Di = new ArrayList<ArrayList<String>>();
		ArrayList<String> t = null;
		for (int i = 0; i < D.size(); i++) {
			t = D.get(i);
			if (t.get(attrIndex).equals(value)) {
				Di.add(t);
			}
		}
		return Di;
	}

	/**
	 * 基于按指定属性划分对D的元组分类所需要的期望信息
	 * 
	 * @param attrIndex
	 *            指定属性的索引
	 * @return 按指定属性划分的期望信息值
	 */
	public double infoAttr(int attrIndex) {
		double info = 0.000;
		ArrayList<String> values = getValues(D, attrIndex);
		DecisionTree dt = new DecisionTree();
		Map<String, Integer> classes; // 获取候选属性中一个取值的（age-> youth-> yes:no）
		double n1 = 0.000;
		for (int i = 0; i < values.size(); i++) {
			double e = 0.0, f = 0.0;
			ArrayList<ArrayList<String>> dv = datasOfValue(attrIndex, values
					.get(i));
			classes = dt.classOfDatas(dv);
			n1 = ((double) dv.size()) / ((double) D.size());
			try {
				/*
				 * e = (double)classes.get("yes"); f =
				 * (double)classes.get("no");
				 */
				e = (double) classes.get("Play");
				f = (double) classes.get("Study");
			} catch (Exception exce) {

			}

			info += n1 * gerException(e, f);
		}
		return info;
	}

	/**
	 * 获取最佳分裂属性的索引
	 * 
	 * @return 最佳分裂属性的索引
	 */
	public int bestGainAttrIndex(double styWhoEx) {
		int index = -1;
		double gain = 0.000;
		double tempGain = 0.000;
		for (int i = 0; i < attrList.size(); i++) {
			tempGain = styWhoEx - infoAttr(i);
			if (tempGain > gain) {
				gain = tempGain;
				index = i;
			}
		}
		return index;
	}

	/**
	 * 获取样本整体期望值
	 * 
	 * @return 样本整体期望值
	 */
	public double getStylebookWholeExpection(Map<String, Integer> classes, int n) {
		double styWhoEx = 0.0;
		Iterator iter = classes.entrySet().iterator();
		for (int i = 0; iter.hasNext(); i++) {
			Map.Entry entry = (Map.Entry) iter.next();
			Integer val = (Integer) entry.getValue();
			double vn = (double) val / (double) n;
			styWhoEx += -(vn) * ((log((double) vn) / (log((double) 2))));
		}
		return styWhoEx;
	}

	/**
	 * 计算属性期望值
	 * 
	 * @return 最佳分裂属性的索引
	 */
	private double gerException(double e, double f) {
		double info = 0.0000;
		if (e == 0.0 || f == 0.0) {
			info = 0.0;
			return info;
		} else if (e == f) {
			info = 1.0;
			return info;
		} else {
			double sum = e + f;
			info = -(e / sum) * ((log((double) (e / sum)) / (log((double) 2))))
					- (f / sum)
					* ((log((double) (f / sum)) / (log((double) 2))));
		}
		return info;
	}

}
