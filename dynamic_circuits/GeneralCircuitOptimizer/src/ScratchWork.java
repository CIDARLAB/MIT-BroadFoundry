import java.util.ArrayList;

public class ScratchWork{
	public static void main(String[] args){
		/*
		ConstantProperties.setCostPerOp();
		ConstantProperties.setOpNames();
		
		Circuit a = new Circuit("a","0011");
		Circuit b = new Circuit("b","0101");
		/*System.out.println(a.getCost());
		System.out.println(a.getName());
		System.out.println(a.getNetlist());
		System.out.println(a.getOperator());
		System.out.println(a.getTruthValue());
		System.out.println("break");
		System.out.println(b.getCost());
		System.out.println(b.getName());
		System.out.println(b.getNetlist());
		System.out.println(b.getOperator());
		System.out.println(b.getTruthValue());
		System.out.println("break");
		ArrayList<Circuit> temp = new ArrayList<Circuit>();
		ArrayList<Circuit> temp2 = new ArrayList<Circuit>();
		temp.add(a);
		temp.add(b);
		temp2.add(b);
		temp2.add(a);
		Circuit c = new Circuit(".",temp);
		System.out.println(c.getCost());
		System.out.println(c.getName());
		System.out.println(c.getNetlist());
		System.out.println(c.getOperator());
		System.out.println(c.getTruthValue());
		System.out.println("break");
		
		Circuit d = new Circuit(".", temp2);

		System.out.println(d.getCost());
		System.out.println(d.getName());
		System.out.println(d.getNetlist());
		System.out.println(d.getOperator());
		System.out.println(d.getTruthValue());
		System.out.println(d.canBeUsed());
		System.out.println("break");
		System.out.println(c.getContainedCircuitNames());
		System.out.println(d.getContainedCircuitNames());
		System.out.println(c.getContainedCircuitNames().equals(d.getContainedCircuitNames()));
		System.out.println("break");
		boolean x = c.equals(d);
		ArrayList<Circuit> list = new ArrayList<Circuit>();
		list.add(c);
		System.out.println(list.contains(d));
		System.out.println(c.equals(d));
		System.out.println(x);
		/*ArrayList<String> a = new ArrayList<String> ();
		a.add("a");
		a.add("b");
		a.add("c");
		a.add("d");
		//a.add("e");
		//a.add("f");
		System.out.println(a);
		ArrayList<ArrayList> combinations = HelperFunctions.combinations(a, 3, true);
		System.out.println(a);
		System.out.println(combinations);*/

		//System.out.println(HelperFunctions.inputTruthValues(4));
		String a = "apple";
		String x = null;
		String y = "a";
		//System.out.println(a.startsWith(x));
		System.out.println(a.startsWith(y));
	}
}