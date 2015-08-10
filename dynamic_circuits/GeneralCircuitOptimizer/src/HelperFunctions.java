import java.util.ArrayList;

public class HelperFunctions{
	public static ArrayList<ArrayList> combinations(ArrayList objects, int groupSize, boolean allowRepeats){
		ArrayList<ArrayList> foundCombinations = new ArrayList<ArrayList>();
		ArrayList tempGroup = new ArrayList();
		addCombinations(foundCombinations,tempGroup,objects,0,groupSize,allowRepeats);
		return foundCombinations;
	}
	public static void addCombinations(ArrayList<ArrayList> foundCombinations, ArrayList tempGroup,ArrayList objects,int startingPoint,int groupSize,boolean allowRepeats){
		if (tempGroup.size()==groupSize){
			ArrayList newTemp = new ArrayList();
			newTemp.addAll(tempGroup);
			foundCombinations.add(newTemp);
			tempGroup.remove(tempGroup.size()-1);
			return;
		}
		for (int i=startingPoint;i<objects.size();i++){
			tempGroup.add(objects.get(i));
			int newStartPoint = i;
			if (!allowRepeats){
				newStartPoint++;
			}
			addCombinations(foundCombinations,tempGroup,objects,newStartPoint,groupSize, allowRepeats);
		}
		if (tempGroup.size()>0){
			tempGroup.remove(tempGroup.size()-1);
		}
		return;
	}
	public static ArrayList<String> inputTruthValues(int numInputs){
		//The first truth value is for the zero circuit, the last one is for the one circuit. The ones in between are for the inputs.
		ArrayList<String> listOfTV = new ArrayList<String> ();
		if (numInputs<0){
			System.out.println("This is not a valid number of inputs");
			return null;
		}
		int tvlength = 1;
		for (int i=0;i<numInputs;i++){
			tvlength = tvlength*2;
		}
		String zero = "0";
		while (zero.length()<tvlength){
			zero += "0";
		}
		listOfTV.add(zero);
		String base = "01";
		for (int i=0;i<numInputs;i++){
			String tempTV = base;
			while(tempTV.length()<tvlength){
				tempTV = tempTV + tempTV;
			}
			listOfTV.add(1, tempTV);
			base = base.substring(0,base.length()/2)+base+base.substring(base.length()/2,base.length());
		}
		
		String one = "1";
		while (one.length()<tvlength){
			one += "1";
		}
		listOfTV.add(one);
		return listOfTV;
	}
	public static String invert(String truthValue){
		String answer = "";
		for (int i=0;i<truthValue.length();i++){
			if(truthValue.charAt(i)=='1'){
				answer += "0";
			}
			else if (truthValue.charAt(i)=='0'){
				answer += "1";
			}
			else{
				answer += "x";
			}
		}
		return answer;
	}
	public static boolean truthValuesAreSame(String tv1, String tv2){
		if (tv1.length()!=tv2.length()){
			return false;
		}
		for (int i=0;i<tv1.length();i++){
			if (tv1.charAt(i)=='1' && tv2.charAt(i)=='0' || tv1.charAt(i)=='0' && tv2.charAt(i)=='1'){
				return false;
			}
		}
		return true;
	}
	
	public static Integer[] getSmallestParenthesis(String circuitString){
		int openParen = -1;
		for (int i=0;i<circuitString.length();i++){
			if(circuitString.charAt(i)=='('){
				openParen = i;
			}
			else if(circuitString.charAt(i)==')' && openParen==-1){
				return null;
			}
			else if(circuitString.charAt(i)==')' && openParen!=-1){
				return new Integer [] {openParen,i};
			}
		}
		return null;
	}
}