import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

public class Circuit{
	
	//A quick way to check whether a circuit can be used to build other circuits or not.
	boolean isGood = true;
	
	//Other properties of a circuit
	String name;
	String truthValue;
	String operator = null;
	double circuitCost;
	int numInputs;
	//The circuit string in netlist form
	ArrayList<String> netlist = null;
	HashSet<String> containedCircuitNames = new HashSet<String>();
	
	//subcircuits is a HashMap that maps strings of subcircuits to their key operator. There will only be one copy of each unique 
	//subcircuit in HashMap. ((a.b).(a.b)) will only have (a.b) once.
	HashMap<String,String> subcircuits = new HashMap<String,String> ();
	
	//Keeps track of each of the truth values in its subcircuits. We do not want to save a circuit if it gives the same truth value
	//as one of its subcircuits.
	HashSet<String> containedTruthValues= new HashSet<String> ();
	
	//Constructors for the class

	//For building inputs
	public Circuit(String symbols, String values) {
		if (values == null|| symbols==null|| (!(values.contains("1")) && !(values.contains("0")))){
			isGood = false;
			return;
		}
		name = symbols;
		truthValue = values;
		containedTruthValues.add(truthValue);
		circuitCost = 0.0;
		numInputs = 0;
		netlist = makeNetlist();
	}
	
	//For building multiple input circuits
	public Circuit (String op, ArrayList<Circuit> circuits) {
		numInputs = circuits.size();
		for (Circuit circ : circuits){
			if (!(circ.canBeUsed()) || circ==null){
				isGood = false;
				return;
			}
		}
		if (circuits.size()==0){
			isGood = false;
			return;
		}
		if (op.equals("~") && circuits.size()>1){
			isGood = false;
			return;
		}
		
		ArrayList<String> inputTruthValues = new ArrayList<String>();
		for (Circuit circ : circuits){
			containedCircuitNames.add(circ.getName());
			inputTruthValues.add(circ.getTruthValue());
		}
		
		//Finds the truth value by using the operator to combine the two subcircuits' truth values.
		if (op.equals("~")){
			truthValue = NOT(inputTruthValues);
		}
		else if (op.equals("&")){
			truthValue = AND(inputTruthValues);
		}
		else if (op.equals("@")){
			truthValue = NAND(inputTruthValues);
		}
		else if (op.equals("+")){
			truthValue = OR(inputTruthValues);
		}
		else if (op.equals("^")){
			truthValue = XOR(inputTruthValues);
		}
		else if (op.equals(".")){
			truthValue = NOR(inputTruthValues);
		}
		else if (op.equals("=")){
			truthValue = XNOR(inputTruthValues);
		}
		/*
		else if (op.equals(">")){
			truthValue = IMPLIES(inputTruthValues);
		}
		else if (op.equals("$")){
			truthValue = NIMPLIES(inputTruthValues);
		}*/
		else{
			isGood = false;
			return;
		}
		
		for (Circuit circ : circuits){
			containedTruthValues.addAll(circ.getContainedTruthValues());
			subcircuits.putAll(circ.getSubcircuits());
			if (circ.getOperator() != null){
				subcircuits.put(circ.getName(), circ.getOperator());
			}
		}
		
		if(truthValue == null || equalsZero() ||equalsOne()||containedTruthValues.contains(truthValue)){
			isGood = false;
			return;
		}
		containedTruthValues.add(truthValue);
		operator = op;
		
		if (operator.equals("~")){
			name = "(" + operator +circuits.get(0).getName() + ")";
		}
		else{
			String tempName = "(";
			for (Circuit circ : circuits){
				tempName = tempName + circ.getName() + operator;
			}
			name = tempName.substring(0, tempName.length()-1) + ")";
		}

		
		circuitCost = calcCost();
		netlist = makeNetlist();
	}
	
/*	
	//When retrieving from a file.
	public Circuit(String circData){
		if (circData.length()==0){
			isGood = false;
			return;
		}
		ArrayList<String> info = new ArrayList<String> (Arrays.asList(circData.split("\t")));
		isGood = Boolean.valueOf(info.get(0));
		name = info.get(1);
		truthValue = info.get(2);
		if (!info.get(3).equals("null")){
			operator = info.get(3);
		}
		circuitCost = Double.valueOf(info.get(4));
		
		String subs = info.get(5);
		subs = subs.substring(1, subs.length()-1);
		if (!subs.equals("")){
			ArrayList<String> pairs = new ArrayList<String> (Arrays.asList(subs.split(", ")));
			for (String keyVal : pairs){
				ArrayList<String> items = new ArrayList<String> (Arrays.asList(keyVal.split("=")));
				subcircuits.put(items.get(0),items.get(1));
			}
		}
		
		subs = info.get(6);
		subs = subs.substring(1, subs.length()-1);
		containedTruthValues = new HashSet<String> (Arrays.asList(subs.split(", ")));
	}
	//Format for saving to a file.
	public String format(){
		String circRep = "";
		circRep = circRep + isGood + "\t";
		circRep = circRep + name + "\t";
		circRep = circRep + truthValue + "\t";
		circRep = circRep + operator + "\t";
		circRep = circRep + circuitCost + "\t";
		circRep = circRep + subcircuits + "\t";
		circRep = circRep + containedTruthValues+"\t";
		return circRep;
	}
*/
	
	//Getter functions
	public boolean canBeUsed(){
		return isGood;
	}
	public String getTruthValue(){
		return truthValue;
	}
	public boolean isTruthValue(String tv){
		//Checks if the truth value of this circuit is the same as the given truth value. Allows for don't cares in the truth value.
		if (tv.length()!=truthValue.length()){
			return false;
		}
		for (int i=0;i<tv.length();i++){
			if (tv.charAt(i)=='1' && truthValue.charAt(i)=='0' || tv.charAt(i)=='0' && truthValue.charAt(i)=='1'){
				return false;
			}
		}
		return true;
	}
	public String getName(){
		return name;
	}
	public HashSet<String> getContainedTruthValues(){
		return containedTruthValues;
	}
	public String getOperator(){
		return operator;
	}
	public HashMap<String, String> getSubcircuits(){
		return subcircuits;
	}
	@Override
	public String toString(){
		return name;
	}
	public ArrayList<String> getNetlist(){
		return netlist;
	}
	public double getCost(){
		return circuitCost;
	}
	public HashSet<String> getContainedCircuitNames(){
		return containedCircuitNames;
	}
	
	//Determines the cost of the circuit
	public double calcCost(){
		double totalCost = 0;
		Set<String> keys = subcircuits.keySet();
		for (String key:keys){
			totalCost += ConstantProperties.costPerOp.get(subcircuits.get(key));
		}
		totalCost += ConstantProperties.costPerOp.get(operator);
		totalCost = roundTwoDecimals(totalCost);
		return totalCost;
	}

	public ArrayList<String> makeNetlist(){
		ArrayList<String> netlist = new ArrayList<String>();
		//Assumes that parentheses are properly matched and converts the circuit string into a netlist.
		String copyName = name;
		int wireCount = 1;
		Integer[] indeces = HelperFunctions.getSmallestParenthesis(copyName);
		while (indeces!=null){
			int openParen = indeces[0];
			int closeParen = indeces[1];
			String wireName = "W"+wireCount;
			String subCircuit = copyName.substring(openParen, closeParen+1);
			copyName = copyName.replace(subCircuit, wireName);
			//Remove front and end parentheses
			subCircuit = subCircuit.substring(1, subCircuit.length()-1);
			String op = null;
			//Find the operator of the subcircuit
			for (String tempOp : ConstantProperties.approvedOperators){
				if (subCircuit.contains(tempOp)){
					op = tempOp;
					break;
				}
			}
			//Add output wire name
			String tempGate = ConstantProperties.opNames.get(op) + "(" + wireName + ",";
			//Add all input wire names
			String[] pieces = subCircuit.split("\\"+op);
			for (String sub : pieces){
				//To deal with NOT
				if (sub.length()==0){
					continue;
				}
				tempGate = tempGate + sub + ",";
			}
			//Remove final comma and close parentheses
			tempGate = tempGate.substring(0, tempGate.length()-1) + ")";
			netlist.add(tempGate);
			wireCount++;
			indeces = HelperFunctions.getSmallestParenthesis(copyName);
		}
		if (netlist.size()==0){
			netlist.add("BUF("+name+")");
		}
		return netlist;
	}
	
	//Used to round the cost of the circuits to two decimal places.
	double roundTwoDecimals(double d) {
        DecimalFormat twoDForm = new DecimalFormat("#.##");
    return Double.valueOf(twoDForm.format(d));
	}
	
	//Checks if two circuits are equal by comparing their names.
	@Override
	public boolean equals(Object obj){
		if (obj==null || getClass() != obj.getClass()){
			return false;
		}
	    final Circuit other = (Circuit) obj;
	    if (containedCircuitNames.size()!=other.getContainedCircuitNames().size()){
	    	return false;
	    }
	    if (containedCircuitNames.size()==0 && other.getContainedCircuitNames().size()==0){
	    	if (name.equals(other.getName())){
	    		return true;
	    	}
	    	else if(!name.equals(other.getName())){
	    		return false;
	    	}
	    	
	    }

		if (name.equals(other.getName())
				|| (other.getContainedCircuitNames().equals(containedCircuitNames)
						&& other.getOperator()==(operator))){
			return true;
		}
		return false;
	}
	
	//Check if it equals one or zero
	boolean equalsZero(){
		for (int i=0;i<truthValue.length();i++){
			if (truthValue.charAt(i)!='0'){
				return false;
			}
		}
		return true;
	}
	boolean equalsOne(){
		for (int i=0;i<truthValue.length();i++){
			if (truthValue.charAt(i)!='1'){
				return false;
			}
		}
		return true;
	}
	//Logic Operators
	String NOT(ArrayList<String> values){
		if (values==null || values.size()!=1){
			isGood = false;
			return null;
		}
		String y = "";
		for (char i:values.get(0).toCharArray()){
			if (i=='1'){
				y+="0";
			}
			else if (i=='0'){;
				y+="1";
			}
		}

		return y;
	}
	String AND(ArrayList<String> values){
		if (values.size()<=1){
			isGood = false;
			return null;
		}
		int numChars = values.get(0).length();
		for (String tv : values){
			if (tv==null || tv.length()!=numChars){
				isGood = false;
				return null;
			}
		}
		
		String y = "";
		
		for (int i=0; i<numChars;i++){
			boolean isTrue = true;
			for(String tv : values){
				if (tv.charAt(i)=='0'){
						isTrue = false;
						break;
				}
			}
			if (isTrue){
				y += "1";
			}
			else if (!isTrue){
				y += "0";
			}
		}
		return y;
		
	}
	String NAND(ArrayList<String> values){
		if (values.size()<=1){
			isGood = false;
			return null;
		}
		int numChars = values.get(0).length();
		for (String tv : values){
			if (tv==null || tv.length()!=numChars){
				isGood = false;
				return null;
			}
		}
		
		String y = "";
		
		for (int i=0; i<numChars;i++){
			boolean isTrue = false;
			for(String tv : values){
				if (tv.charAt(i)=='0'){
						isTrue = true;
						break;
				}
			}
			if (isTrue){
				y += "1";
			}
			else if (!isTrue){
				y += "0";
			}
		}
		return y;
		
	}
	String OR(ArrayList<String> values){
		if (values.size()<=1){
			isGood = false;
			return null;
		}
		int numChars = values.get(0).length();
		for (String tv : values){
			if (tv==null || tv.length()!=numChars){
				isGood = false;
				return null;
			}
		}
		
		String y = "";
		
		for (int i=0; i<numChars;i++){
			boolean isTrue = false;
			for(String tv : values){
				if (tv.charAt(i)=='1'){
						isTrue = true;
						break;
				}
			}
			if (isTrue){
				y += "1";
			}
			else if (!isTrue){
				y += "0";
			}
		}
		return y;
		
	}
	String XOR(ArrayList<String> values){
		if (values.size()<=1){
			isGood = false;
			return null;
		}
		int numChars = values.get(0).length();
		for (String tv : values){
			if (tv==null || tv.length()!=numChars){
				isGood = false;
				return null;
			}
		}
		
		String y = "";
		
		for (int i=0; i<numChars;i++){
			boolean isTrue = false;
			for(String tv : values){
				if (tv.charAt(i)=='1' && isTrue==false){
						isTrue = true;
				}
				else if (tv.charAt(i)=='1' && isTrue==true){
					isTrue = false;
					break;
				}
			}
			if (isTrue){
				y += "1";
			}
			else if (!isTrue){
				y += "0";
			}
		}
		return y;
	}
	String NOR(ArrayList<String> values){
		if (values.size()<=1){
			isGood = false;
			return null;
		}
		int numChars = values.get(0).length();
		for (String tv : values){
			if (tv==null || tv.length()!=numChars){
				isGood = false;
				return null;
			}
		}
		
		String y = "";
		
		for (int i=0; i<numChars;i++){
			boolean isTrue = true;
			for(String tv : values){
				if (tv.charAt(i)=='1'){
						isTrue = false;
						break;
				}
			}
			if (isTrue){
				y += "1";
			}
			else if (!isTrue){
				y += "0";
			}
		}
		return y;
		
	}
	String XNOR(ArrayList<String> values){
		if (values.size()<=1){
			isGood = false;
			return null;
		}
		int numChars = values.get(0).length();
		for (String tv : values){
			if (tv==null || tv.length()!=numChars){
				isGood = false;
				return null;
			}
		}
		
		String y = "";
		
		for (int i=0; i<numChars;i++){
			boolean isTrue = true;
			for(String tv : values){
				if (tv.charAt(i)=='1' && isTrue==true){
						isTrue = false;
				}
				else if (tv.charAt(i)=='1' && isTrue==false){
					isTrue = true;
					break;
				}
			}
			if (isTrue){
				y += "1";
			}
			else if (!isTrue){
				y += "0";
			}
		}
		return y;
	}
	
	/*
	String IMPLIES(String alphaVals, String betaVals){
		if (alphaVals==null || betaVals==null || alphaVals.length() != betaVals.length()){
			isGood = false;
			return null;
		}
		int x = alphaVals.length();
		String y = "";
		for (int i=0;i<x;i++){
			if (alphaVals.charAt(i)=='0' || betaVals.charAt(i)=='1'){
				y+="1";
			}
			else if (alphaVals.charAt(i)=='1' && betaVals.charAt(i)=='0'){
				y+="0";
			}
		}
		return y;
		
	}
	String NIMPLIES(String alphaVals, String betaVals){
		if (alphaVals==null || betaVals==null || alphaVals.length() != betaVals.length()){
			isGood = false;
			return null;
		}
		int x = alphaVals.length();
		String y = "";
		for (int i=0;i<x;i++){
			if (alphaVals.charAt(i)=='0' || betaVals.charAt(i)=='1'){
				y+="0";
			}
			else if (alphaVals.charAt(i)=='1' && betaVals.charAt(i)=='0'){
				y+="1";
			}
		}
		return y;
		
	}
	*/

}